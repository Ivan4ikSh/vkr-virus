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
    int shift = 2;
    double asym = 0.5;           // asymerty
    double R0 = 2.5;             // Basic reproduction number
    double Ub = 1e-3;            // Mutation rate
    double a = 7.0;              // Cross-immunity scale
    double dt = 0.5;             // Time step
    double init_infected = 1e-3; // Initial fraction of infected

    double r_i = 1; // Initial fraction of infected
    double r_r = 1; // Initial fraction of infected
    int64_t N = 1e8;             // Total population size
};

// ================== Main Simulation Class ==================
class EpidemicSimulator {
private:
    ModelParams p_;
    mt19937 rng_;
    double angle_;

    vector<vector<double>> I_;
    vector<vector<double>> R_;

    vector<double> norm_;
    vector<double> finf_;

    vector<vector<double>> K2D_;

    fftw_complex* fft_in_;
    fftw_complex* fft_out_;
    fftw_complex* fft_kernel_;
    fftw_plan plan_inv_;
    fftw_plan plan_fwd_;
public:
    EpidemicSimulator(const ModelParams& params, int seed = 42) : p_(params), rng_(seed) {
        std::uniform_real_distribution<double> dist_angle(0.0, 2.0 * acos(-1.0));
        angle_ = dist_angle(rng_);

        PrecomputeKernel();
        InitFFT();
        InitializeState();
        norm_.resize(p_.M, 0.0);
        finf_.resize(p_.M, 0.0);
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

    ~EpidemicSimulator() {
        CleanFFT(); // Не забываем освобождать память
    }

private:
    // Precompute 1D kernel K(dx) = 1 - exp(-dx / a)
    void PrecomputeKernel() {
        int L = p_.L;
        double asym = p_.asym;
        K2D_.assign(L, vector<double>(L, 0.0));

        double cos_a = cos(angle_);
        double sin_a = sin(angle_);

        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                double dx_signed = (x > L / 2) ? x - L : x;
                double dy_signed = (y > L / 2) ? y - L : y;

                double dx_rot = dx_signed * cos_a + dy_signed * sin_a;
                double dy_rot = -dx_signed * sin_a + dy_signed * cos_a;

                double dist = sqrt(dx_rot * dx_rot + asym * dy_rot * dy_rot);

                K2D_[x][y] = 1.0 - exp(-dist / p_.a);
            }
        }
    }

    // Initialization
    void InitializeState() {
        int L = p_.L;
        int shift = p_.shift;
        double asym = p_.asym;

        double r_i = p_.r_i;
        double r_r = p_.r_r;

        I_.assign(L, vector<double>(L, 0.0));
        R_.assign(L, vector<double>(L, 0.0));

        double cos_a = cos(angle_);
        double sin_a = sin(angle_);

        double center = L / 2.0;

        double cx_i = center + shift * cos_a;
        double cy_i = center + shift * sin_a;

        double cx_r = center - shift * cos_a;
        double cy_r = center - shift * sin_a;

        double sum_I = 0.0;
        double sum_R = 0.0;

        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                double dx = x - cx_i;
                double dy = y - cy_i;
                double dist2 = dx * dx + asym * dy * dy;
                I_[x][y] = exp(-dist2 / (2.0 * r_i * r_i));
                sum_I += I_[x][y];
            }
        }

        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                double dx = x - cx_r;
                double dy = y - cy_r;
                double dist2 = dx * dx + asym * dy * dy;
                R_[x][y] = exp(-dist2 / (2.0 * r_r * r_r));
                sum_R += R_[x][y];
            }
        }

        // Нормализация (остается без изменений)
        double target_R = 1.0 - p_.init_infected;
        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                R_[x][y] = (R_[x][y] / sum_R) * target_R;
            }
        }

        double target_I = p_.init_infected;
        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                I_[x][y] = (I_[x][y] / sum_I) * target_I;
            }
        }
    }

    void InitFFT() {
        int L = p_.L;
        int N = L * L;

        fft_in_ = fftw_alloc_complex(N);
        fft_out_ = fftw_alloc_complex(N);
        fft_kernel_ = fftw_alloc_complex(N);

        // Создаем планы прямого и обратного преобразования
        plan_fwd_ = fftw_plan_dft_2d(L, L, fft_in_, fft_out_, FFTW_FORWARD, FFTW_MEASURE);
        plan_inv_ = fftw_plan_dft_2d(L, L, fft_in_, fft_out_, FFTW_BACKWARD, FFTW_MEASURE);

        // Считаем FFT для ядра K2D_ (оно не меняется со временем, делаем это 1 раз)
        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                fft_in_[x * L + y][0] = K2D_[x][y];
                fft_in_[x * L + y][1] = 0.0;
            }
        }
        fftw_execute(plan_fwd_);

        // Сохраняем результат
        for (int i = 0; i < N; ++i) {
            fft_kernel_[i][0] = fft_out_[i][0];
            fft_kernel_[i][1] = fft_out_[i][1];
        }
    }

    // Single time step
    void Step() {
        int L = p_.L;
        int N_grid = L * L;
        double dt = p_.dt;
        double R0 = p_.R0;
        double Ub = p_.Ub;

        vector<vector<double>> Q(L, vector<double>(L, 0.0));
        vector<vector<double>> P(L, vector<double>(L, 0.0));

        // 1. Вычисляем Q = K * R
        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                fft_in_[x * L + y][0] = R_[x][y];
                fft_in_[x * L + y][1] = 0.0;
            }
        }
        fftw_execute(plan_fwd_); // Прямое FFT от R_

        for (int i = 0; i < N_grid; ++i) {
            // Комплексное умножение: (a+bi)*(c+di) = (ac-bd) + (ad+bc)i
            double real1 = fft_out_[i][0];
            double imag1 = fft_out_[i][1];
            double real2 = fft_kernel_[i][0];
            double imag2 = fft_kernel_[i][1];
            fft_in_[i][0] = real1 * real2 - imag1 * imag2;
            fft_in_[i][1] = real1 * imag2 + imag1 * real2;
        }
        fftw_execute(plan_inv_); // Обратное FFT

        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                Q[x][y] = fft_out_[x * L + y][0] / N_grid; // Нормализация
            }
        }

        // 2. Вычисляем P = K * I
        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                fft_in_[x * L + y][0] = I_[x][y];
                fft_in_[x * L + y][1] = 0.0;
            }
        }
        fftw_execute(plan_fwd_); // Прямое FFT от I_

        for (int i = 0; i < N_grid; ++i) {
            double real1 = fft_out_[i][0], imag1 = fft_out_[i][1];
            double real2 = fft_kernel_[i][0], imag2 = fft_kernel_[i][1];
            fft_in_[i][0] = real1 * real2 - imag1 * imag2;
            fft_in_[i][1] = real1 * imag2 + imag1 * real2;
        }
        fftw_execute(plan_inv_); // Обратное FFT

        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                P[x][y] = fft_out_[x * L + y][0] / N_grid; // Нормализация
            }
        }

        vector<vector<double>> I_new(L, vector<double>(L, 0.0));
        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                // Обновление R
                R_[x][y] = R_[x][y] * (1.0 - dt * R0 * P[x][y]) + dt * I_[x][y];

                // Обновление I
                double newI = I_[x][y] * (1.0 + dt * (R0 * Q[x][y] - 1.0));
                I_new[x][y] = GetStochasticCorrection(newI);

                // Мутации (диффузия)
                int y_up = (y - 1 + L) % L;
                int y_down = (y + 1) % L;
                int x_left = (x - 1 + L) % L;
                int x_right = (x + 1) % L;

                double in_mut = dt * Ub * (I_[x][y_up] + I_[x][y_down] + I_[x_left][y] + I_[x_right][y]);
                double out_mut = dt * Ub * 4.0 * I_[x][y];

                I_new[x][y] += GetStochasticCorrection(in_mut) - out_mut;
                if (I_new[x][y] * p_.N < 1.0) I_new[x][y] = 0.0;
                if (I_new[x][y] < 0) I_new[x][y] = 0;
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

    // Calculate statistics for the current step
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

    void CleanFFT() {
        fftw_destroy_plan(plan_fwd_);
        fftw_destroy_plan(plan_inv_);
        fftw_free(fft_in_);
        fftw_free(fft_out_);
        fftw_free(fft_kernel_);
    }

    // Save current state (I and R matrices)
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

    // Save time series and parameters
    void SaveFinalResults() {
        const std::string& output_dir = CONST::OUTPUT + "/" + NAME;
        const std::string& data_dir = output_dir + "/" + CONST::DATA;
        fs::create_directories(data_dir);

        // Norm file
        std::ofstream norm_file(data_dir + "/norm.csv");
        norm_file << "step;norm\n";
        for (int i = 0; i < p_.M; ++i)
            norm_file << i << ";" << norm_[i] << "\n";

        // Infection fraction file
        std::ofstream finf_file(data_dir + "/finf.csv");
        finf_file << "step;finf\n";
        for (int i = 0; i < p_.M; ++i)
            finf_file << i << ";" << finf_[i] << "\n";

        // Parameters file
        std::ofstream param_file(data_dir + "/parameters.txt");
        param_file << "L = " << p_.L << " ";
        param_file << "M = " << p_.M << " ";
        param_file << "R0 = " << p_.R0 << " ";
        param_file << "Ub = " << p_.Ub << "\n";
        param_file << "a = " << p_.a << " ";
        param_file << "dt = " << p_.dt << " ";
        param_file << "N = " << p_.N << " ";
        param_file << "Seed = " << p_.seed;
    }

    // Progress output
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
    // Parameters can be overridden here
    params.L = 100;
    params.M = 1000;
    params.tshow = params.M / 50;
    params.T0 = 0;
    params.asym = 1.0;
    params.R0 = 1.8;
    params.Ub = 1e-3;
    params.a = 14.0;
    params.dt = 0.5;
    params.init_infected = 1e-3;
    params.N = 1e8;
    params.shift = 0;

    params.r_i = 5;
    params.r_r = 5;

    params.seed = 1;
    TEST(params, "exp1");

    params.r_i = 3;
    params.r_r = 5;
    params.seed = 2;
    TEST(params, "exp2");

    params.asym = 0.3;
    params.shift = 10;
    params.seed = 3;
    TEST(params, "exp31");

    params.M = 5000;
    params.tshow = params.M / 50;
    TEST(params, "exp32");
    return 0;
}
