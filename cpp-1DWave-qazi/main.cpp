#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <fftw3.h>

using namespace std;
namespace fs = filesystem;

namespace CONST {
    const double& EPS = 1e-12;
    const double& STOCH_THRESHOLD = 10.0;

    const std::string OUTPUT = "out";
    const std::string DATA = "data";
}

std::string NAME = "test";

// ================== Model Parameters ==================
struct ModelParams {
    int seed = 42;               // Seed for rng
    int L = 25;                  // Grid size L x L
    int M = 2000;                // Number of time steps
    int tshow = M / 50;          // State saving interval
    int T0 = 0;                  // First step to save
    double R0 = 2.5;             // Basic reproduction number
    double Ub = 1e-3;            // Mutation rate
    double a = 7.0;              // Cross-immunity scale
    double dt = 0.5;             // Time step
    double init_infected = 1e-3; // Initial fraction of infected
    int64_t N = 1e8;             // Total population size
};

// ================== Main Simulation Class ==================
class EpidemicSimulator {
private:
    ModelParams p_;
    mt19937 rng_;
    uniform_real_distribution<double> uniform_{ 0.0, 1.0 };

    // States (fractions of the total population)
    std::vector<std::vector<double>> I_;   // Infected
    std::vector<std::vector<double>> R_;   // Recovered/Susceptible

    // Time series
    std::vector<double> norm_;   // Sum of I+R (should be ~1.0)
    std::vector<double> finf_;   // Fraction of infected in total population

    // Precomputed 1D kernel Kx[dx] = 1 - exp(-dx / a)
    std::vector<double> K2D_;

    // FFTW variables
    fftw_complex* fft_in_;
    fftw_complex* fft_out_;
    fftw_complex* fft_kernel_;
    fftw_plan plan_inv_;
    fftw_plan plan_fwd_;

public:
    EpidemicSimulator(const ModelParams& params, int seed = 42) : p_(params), rng_(seed) {
        PrecomputeKernel();
        InitFFT();
        InitializeState();
        norm_.resize(p_.M, 0.0);
        finf_.resize(p_.M, 0.0);
    }

    ~EpidemicSimulator() {
        CleanFFT();
    }

    void Run() {
        PrintParameters();
        auto start = chrono::high_resolution_clock::now();
        for (int step = 0; step < p_.M; ++step) {
            if (step >= p_.T0 && step % p_.tshow == 0) SaveState(step);

            CalculateStatistics(step);
            Step();

            if (step % max(1, p_.M / 10) == 0) PrintProgress(step);
        }
        auto end = chrono::high_resolution_clock::now();
        double elapsed = chrono::duration<double>(end - start).count();
        std::cout << "\nSimulation finished in " << elapsed << " s\n";

        SaveFinalResults();
    }

private:
    // Precompute 1D kernel K(dx) = 1 - exp(-dx / a)
    void PrecomputeKernel() {
        K2D_.resize(p_.L);
        for (int dx = 0; dx < p_.L; ++dx) {
            int d = std::min(dx, p_.L - dx); // periodic distance
            K2D_[dx] = 1.0 - exp(-static_cast<double>(d) / p_.a);
        }
    }

    void InitFFT() {
        int L = p_.L;

        fft_in_ = fftw_alloc_complex(L);
        fft_out_ = fftw_alloc_complex(L);
        fft_kernel_ = fftw_alloc_complex(L);

        plan_fwd_ = fftw_plan_dft_1d(L, fft_in_, fft_out_, FFTW_FORWARD, FFTW_MEASURE);
        plan_inv_ = fftw_plan_dft_1d(L, fft_in_, fft_out_, FFTW_BACKWARD, FFTW_MEASURE);

        for (int x = 0; x < L; ++x) {
            fft_in_[x][0] = K2D_[x];
            fft_in_[x][1] = 0.0;
        }
        fftw_execute(plan_fwd_);

        for (int i = 0; i < L; ++i) {
            fft_kernel_[i][0] = fft_out_[i][0];
            fft_kernel_[i][1] = fft_out_[i][1];
        }
    }

    void CleanFFT() {
        fftw_destroy_plan(plan_fwd_);
        fftw_destroy_plan(plan_inv_);
        fftw_free(fft_in_);
        fftw_free(fft_out_);
        fftw_free(fft_kernel_);
    }

    // Initialization
    void InitializeState() {
        int L = p_.L;
        I_.assign(L, std::vector<double>(L, 0.0));
        R_.assign(L, std::vector<double>(L, 0.0));

        int x0 = 15;
        double infected_total = p_.init_infected;
        int left_cells = x0 * L;

        for (int y = 0; y < L; ++y) {
            for (int x = 0; x < L; ++x) {
                if (x == x0) {
                    I_[x][y] = infected_total / L;
                    R_[x][y] = 0.0;
                }
                else if (x < x0) {
                    I_[x][y] = 0.0;
                    R_[x][y] = (1.0 - infected_total) / left_cells;
                }
                else {
                    I_[x][y] = 0.0;
                    R_[x][y] = 0.0;
                }
            }
        }

        double sumI = 0.0, sumR = 0.0;
        for (int y = 0; y < L; ++y)
            for (int x = 0; x < L; ++x) {
                sumI += I_[x][y];
                sumR += R_[x][y];
            }
        std::cout << "Initialization: sumI = " << sumI << ", sumR = " << sumR << ", total = " << sumI + sumR << "\n";
    }

    // Single time step
    void Step() {
        int L = p_.L;
        double dt = p_.dt;
        double R0 = p_.R0;
        double Ub = p_.Ub;

        std::vector<double> Qx(L, 0.0);
        std::vector<double> Px(L, 0.0);

        // 1. Analog of sum(R, 1) and sum(I, 1)
        std::vector<double> R_sum(L, 0.0);
        std::vector<double> I_sum(L, 0.0);
        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                R_sum[x] += R_[x][y];
                I_sum[x] += I_[x][y];
            }
        }

        // 2. FFT Convolution for R_sum -> Qx
        for (int x = 0; x < L; ++x) {
            fft_in_[x][0] = R_sum[x];
            fft_in_[x][1] = 0.0;
        }
        fftw_execute(plan_fwd_);

        for (int i = 0; i < L; ++i) {
            double r1 = fft_out_[i][0], i1 = fft_out_[i][1];
            double r2 = fft_kernel_[i][0], i2 = fft_kernel_[i][1];
            fft_in_[i][0] = r1 * r2 - i1 * i2;
            fft_in_[i][1] = r1 * i2 + i1 * r2;
        }
        fftw_execute(plan_inv_);

        for (int x = 0; x < L; ++x) {
            Qx[x] = fft_out_[x][0] / L;
        }

        // 3. FFT Convolution for I_sum -> Px
        for (int x = 0; x < L; ++x) {
            fft_in_[x][0] = I_sum[x];
            fft_in_[x][1] = 0.0;
        }
        fftw_execute(plan_fwd_);

        for (int i = 0; i < L; ++i) {
            double r1 = fft_out_[i][0], i1 = fft_out_[i][1];
            double r2 = fft_kernel_[i][0], i2 = fft_kernel_[i][1];
            fft_in_[i][0] = r1 * r2 - i1 * i2;
            fft_in_[i][1] = r1 * i2 + i1 * r2;
        }
        fftw_execute(plan_inv_);

        for (int x = 0; x < L; ++x) {
            Px[x] = fft_out_[x][0] / L;
        }

        // - Updates and Mutations -
        std::vector<std::vector<double>> I_new(L, std::vector<double>(L, 0.0));

        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                R_[x][y] = R_[x][y] * (1.0 - dt * R0 * Px[x]) + dt * I_[x][y];

                double newI = I_[x][y] * (1.0 + dt * (R0 * Qx[x] - 1.0));
                I_new[x][y] = GetStochasticCorrection(newI);

                int y_up = (y - 1 + L) % L;
                int y_down = (y + 1) % L;
                int x_left = (x - 1 + L) % L;
                int x_right = (x + 1) % L;

                double in_mut = dt * Ub * (I_[x][y_up] + I_[x][y_down] + I_[x_left][y] + I_[x_right][y]);
                double out_mut = dt * Ub * 4 * I_[x][y];

                I_new[x][y] += GetStochasticCorrection(in_mut) - out_mut;
                if (I_new[x][y] * p_.N < 1.0) I_new[x][y] = 0.0;
            }
        }

        I_ = move(I_new);
    }

    double GetStochasticCorrection(const double& val) {
        const int64_t& N = p_.N;
        const double& lambda = val * N;
        if (lambda >= 10.0) return val;

        double new_val = 0.0;
        const int& xm = static_cast<int>(round(6.0 * lambda));
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        int count = 0;
        for (int n = 0; n < xm; ++n) {
            if (xm * dist(rng_) < lambda) ++count;
        }

        return static_cast<double>(count) / N;
    }

    void CalculateStatistics(int step) {
        double sumI = 0.0, sumR = 0.0;
        for (int y = 0; y < p_.L; ++y)
            for (int x = 0; x < p_.L; ++x) {
                sumI += I_[x][y];
                sumR += R_[x][y];
            }
        double total = sumI + sumR;
        norm_[step] = total;
        finf_[step] = (total > 0.0) ? sumI / total : 0.0;
    }

    void SaveState(int step) {
        const std::string& output_dir = CONST::OUTPUT + "/" + NAME;
        const std::string& data_dir = output_dir + "/" + CONST::DATA;
        fs::create_directories(data_dir);

        auto save_matrix = [&](const std::vector<std::vector<double>>& mat, const std::string& fname) {
            std::ofstream fout(data_dir + "/" + fname);
            fout << scientific << setprecision(16);
            for (int y = 0; y < p_.L; ++y) {
                for (int x = 0; x < p_.L; ++x) {
                    fout << mat[x][y];
                    if (x < p_.L - 1) fout << ";";
                }
                fout << "\n";
            }
            };

        save_matrix(I_, "I_step_" + to_string(step) + ".csv");
        save_matrix(R_, "R_step_" + to_string(step) + ".csv");
    }

    void SaveFinalResults() {
        const std::string& output_dir = CONST::OUTPUT + "/" + NAME;
        const std::string& data_dir = output_dir + "/" + CONST::DATA;
        fs::create_directories(data_dir);

        std::ofstream norm_file(data_dir + "/norm.csv");
        norm_file << "step;norm\n";
        for (int i = 0; i < p_.M; ++i) norm_file << i << ";" << norm_[i] << "\n";

        std::ofstream finf_file(data_dir + "/finf.csv");
        finf_file << "step;finf\n";
        for (int i = 0; i < p_.M; ++i) finf_file << i << ";" << finf_[i] << "\n";

        std::ofstream param_file(data_dir + "/parameters.txt");
        param_file << "L = " << p_.L << " ";
        param_file << "M = " << p_.M << " ";
        param_file << "R0 = " << p_.R0 << " ";
        param_file << "Ub = " << p_.Ub << "\n";
        param_file << "a = " << p_.a << " ";
        param_file << "N = " << p_.N << " ";
        param_file << "Seed = " << p_.seed << "\n";
    }

    void PrintProgress(int step) {
        std::cout << "Step " << step << "/" << p_.M << " (" << (100 * step / p_.M) << "%)" << " | norm = " << norm_[step] << ", finf = " << finf_[step] << "\n";
    }

    void PrintParameters() {
        std::cout << "Grid size: " << p_.L << " x " << p_.L << "\n";
        std::cout << "Time steps: " << p_.M << "\n";
        std::cout << "R0 = " << p_.R0 << "\n";
        std::cout << "Ub = " << p_.Ub << "\n";
        std::cout << "a = " << p_.a << "\n";
        std::cout << "dt = " << p_.dt << "\n";
        std::cout << "init_infected = " << p_.init_infected << "\n";
        std::cout << "N = " << p_.N << "\n";
        std::cout << "Seed = " << p_.seed << "\n\n";
    }
};

void InitDirectory() {
    const std::string& output_dir = CONST::OUTPUT + "/" + NAME;
    const std::string& data_dir = output_dir + "/" + CONST::DATA;
    fs::create_directories(data_dir);

    for (const auto& entry : fs::directory_iterator(data_dir)) {
        fs::remove_all(entry.path());
    }
}

void TEST(const ModelParams& params, const std::string& name) {
    NAME = name;
    InitDirectory();
    EpidemicSimulator sim(params, params.seed);
    sim.Run();

    const std::string& output_dir = CONST::OUTPUT + "/" + NAME;
    const std::string& data_dir = output_dir + "/" + CONST::DATA;
    std::string cmd = "python 2D_visualize.py -d " + data_dir + " -s " +
        output_dir + "/" + NAME + "_snap.png -a " +
        output_dir + "/" + NAME + "_evo.gif -f" +
        output_dir + "/" + NAME + "_finf.png";
    system(cmd.c_str());
}

// ================== Entry Point ==================
int main() {
    ModelParams params;
    params.L = 50;
    params.T0 = 0;
    params.R0 = 1.8;
    params.Ub = 1e-3;
    params.a = 6.0;
    params.dt = 0.5;
    params.init_infected = 1e-3;
    params.N = 1e8;

    params.M = 500;
    params.tshow = params.M / 50;
    params.seed = 1;
    TEST(params, "test_short");

    params.M = 5000;
    params.tshow = params.M / 50;
    params.seed = 2;
    TEST(params, "test_long");

    return 0;
}
