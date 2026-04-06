#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <numeric>
#include <string>
#include <filesystem>

struct Params {
    double R0;
    double Ub;
    double a;
    int L;
    double N;
    double Tmax;
    double T0;
    double Txy;
    double dt;
    unsigned seed;
};

namespace fs = std::filesystem;
const std::string OUTPUT_DIR = "out";

// --- Helper Functions ---
double GetStochasticCorrection(const double& val, const double& N, std::mt19937& gen) {
    const double& lambda = val * N;
    if (lambda >= 10.0) return val;

    double new_val = 0.0;
    int xm = static_cast<int>(round(6.0 * lambda));
    if (xm == 0) return 0.0;
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    int count = 0;
    for (int n = 0; n < xm; ++n) {
        if (xm * dist(gen) < lambda) ++count;
    }

    return static_cast<double>(count) / N;
}

std::vector<double> PrecomputeKernel(const int& L, const double& a) {
    std::vector<double> k(L);
    for (int u = 0; u < L; ++u) {
        k[u] = 1.0 - std::exp(-static_cast<double>(u) / a);
    }
    return k;
}

void SaveToFile(const std::string& filename, const std::vector<double>& data) {
    std::ofstream f(OUTPUT_DIR + "/" + filename);
    for (double v : data) f << v << "\n";
}

void WaveSimulation(const Params& params) {
    double R0 = params.R0;
    double Ub = params.Ub;
    double a = params.a;
    int L = params.L;
    double N = params.N;
    double Tmax = params.Tmax;
    double T0 = params.T0;
    double Txy = params.Txy;
    double dt = params.dt;
    unsigned seed = params.seed;

    fs::create_directories(OUTPUT_DIR);
    const int& M = static_cast<int>(round(Tmax / dt));
    std::mt19937 gen(seed);

    std::vector<double> I_buf1(L, 0.0);
    std::vector<double> I_buf2(L, 0.0);
    std::vector<double> R_buf1(L, 0.0);
    std::vector<double> R_buf2(L, 0.0);
    auto* I_cur = &I_buf1, * I_next = &I_buf2;
    auto* R_cur = &R_buf1, * R_next = &R_buf2;

    int start_node = 35;
    (*I_cur)[start_node] = 1e-3;
    const double& init_R = (1.0 - 1e-3) / start_node;
    for (int j = 0; j < start_node; ++j) (*R_cur)[j] = init_R;

    const std::vector<double>& K = PrecomputeKernel(L, a);

    std::vector<double> norm(M);
    std::vector<double> finf(M);
    std::vector<double> slice_times;
    std::vector<std::vector<double>> I_slices;
    std::vector<std::vector<double>> R_slices;
    const int& show_step = static_cast<int>(round((Tmax / 7.0) / dt));
    const int& t0_step = static_cast<int>(round(T0 / dt)) - 1;
    double meanI0 = 0.0;

    for (int i = 0; i < M; ++i) {
        const double& t = (i + 1) * dt;
        const auto& I_old = *I_cur;
        const auto& R_old = *R_cur;
        auto& Inew = *I_next;
        auto& Rnew = *R_next;

        const double& sumI_old = std::accumulate(I_old.begin(), I_old.end(), 0.0);
        const double& norm_i = sumI_old + std::accumulate(R_old.begin(), R_old.end(), 0.0);
        norm[i] = norm_i;
        finf[i] = sumI_old / norm_i;

        for (int j = 0; j < L; ++j) {
            double Qj = 0.0, Pj = 0.0;
            for (int k = 0; k <= j; ++k) Qj += R_old[k] * K[j - k];
            for (int k = j; k < L; ++k) Pj += I_old[k] * K[k - j];

            if (t > Txy) {
                Qj = 0.0; Pj = 0.0;
                for (int k = j; k < L; ++k) Qj += R_old[k] * K[k - j];
                for (int k = 0; k <= j; ++k) Pj += I_old[k] * K[j - k];
            }

            const double& newI = I_old[j] * (1.0 + dt * (R0 * Qj - 1.0));
            Inew[j] = GetStochasticCorrection(newI, N, gen);
            Rnew[j] = R_old[j] * (1.0 - dt * R0 * Pj) + dt * I_old[j];

            if (j > 0 && j < L - 1) {
                const double& in = dt * Ub * (I_old[j + 1] + I_old[j - 1]);
                const double& out = dt * Ub * 2.0 * I_old[j];
                Inew[j] += GetStochasticCorrection(in, N, gen) - out;
            }
        }

        Inew[0] = I_old[0] + dt * Ub * (I_old[1] - I_old[0]);
        Inew[L - 1] = I_old[L - 1] + dt * Ub * (I_old[L - 2] - I_old[L - 1]);

        if (norm_i > 0) for (double& val : Inew) val /= norm_i;

        if (i == t0_step && sumI_old > 0) {
            double sum_pos = 0.0, sI = 0.0;
            for (int j = 0; j < L; ++j) { sum_pos += Inew[j] * (j + 1); sI += Inew[j]; }
            meanI0 = sum_pos / sI;
        }

        if (t > T0 && (i % show_step == 0)) {
            slice_times.push_back(t);
            I_slices.push_back(Inew);
            R_slices.push_back(Rnew);
        }

        std::swap(I_cur, I_next);
        std::swap(R_cur, R_next);
        if (i % 500 == 0) std::cout << "t = " << t << " | norm = " << norm_i << " | finf = " << finf[i] << "\n";
    }

    const auto& I_final = *I_cur;
    const auto& R_final = *R_cur;
    const double& sumI = std::accumulate(I_final.begin(), I_final.end(), 0.0);
    const double& sumR = std::accumulate(R_final.begin(), R_final.end(), 0.0);

    double meanI = 0.0, meanR = 0.0, sdI = 0.0, sdR = 0.0, speed = 0.0;

    // ИСПРАВЛЕНИЕ 2: Защита от NaN
    if (sumI > 0) {
        for (int j = 0; j < L; ++j) {
            meanI += I_final[j] * (j + 1);
            meanR += R_final[j] * (j + 1);
        }
        meanI /= sumI;
        meanR /= sumR;

        speed = (meanI - meanI0) / (Tmax - T0);
        for (int j = 0; j < L; ++j) {
            sdI += I_final[j] * ((j + 1) - meanI) * ((j + 1) - meanI);
            sdR += R_final[j] * ((j + 1) - meanR) * ((j + 1) - meanR);
        }
        sdI = sqrt(sdI / sumI);
        sdR = sqrt(sdR / sumR);
    }

    SaveToFile("norm.txt", norm);
    SaveToFile("finf.txt", finf);
    SaveToFile("I_final.txt", I_final);
    SaveToFile("slice_times.txt", slice_times);

    for (size_t idx = 0; idx < I_slices.size(); ++idx) {
        SaveToFile("I_slice_" + std::to_string(idx) + ".txt", I_slices[idx]);
        SaveToFile("R_slice_" + std::to_string(idx) + ".txt", R_slices[idx]);
    }

    std::ofstream f_obs(OUTPUT_DIR + "/observables.txt");
    f_obs << "speed=" << speed << "\n";
    f_obs << "sdI=" << sdI << "\n";
    f_obs << "sdR=" << sdR << "\n";
    f_obs << "finf=" << finf.back() << "\n";
    std::cout << "speed=" << speed << "\n";
    std::cout << "sdI=" << sdI << "\n";
    std::cout << "sdR=" << sdR << "\n";
    std::cout << "finf=" << finf.back() << "\n";
    f_obs.close();

    std::cout << "\nSimulation Complete. Speed: " << speed << "\n";
}

int main() {
    Params params;
    params.R0 = 1.8;
    params.Ub = 1e-3;
    params.a = 14.0;
    params.L = 220;
    params.N = 1e8;
    params.Tmax = 5000;
    params.T0 = 500;
    params.Txy = 5000;
    params.dt = 0.5;
    params.seed = 1;

    WaveSimulation(params);

    if (std::system("python 1D_visualize.py") != 0) std::cerr << "Visualizer failed to start.\n";
    return 0;
}