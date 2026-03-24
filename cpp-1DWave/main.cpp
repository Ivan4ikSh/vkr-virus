#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <numeric>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;
const std::string OUTPUT_DIR = "out";

// --- Helper Functions ---
double GetStochasticCorrection(const double& lambda, const double& N, std::mt19937& gen) {
    if (lambda >= 10.0) return lambda;
    double new_val = 0.0;
    int xm = static_cast<int>(round(6.0 * lambda));
    if (xm > 0) {
        double prob = lambda / xm;
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        int count = 0;
        for (int n = 0; n < xm; ++n) {
            if (dist(gen) < prob) ++count;
        }
        new_val = static_cast<double>(count) / N;
    }
    return new_val;
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

void WaveSimulation(const double& R0, const double& Ub, const double& a, const int& L,
    const double& N, const double& Tmax, const double& T0,
    const double& Txy, const double& dt, const unsigned& seed = 42) {

    fs::create_directories(OUTPUT_DIR);
    const int& M = static_cast<int>(round(Tmax / dt));
    std::mt19937 gen(seed);

    // Buffers and Pointers
    std::vector<double> I_buf1(L, 0.0), I_buf2(L, 0.0);
    std::vector<double> R_buf1(L, 0.0), R_buf2(L, 0.0);
    auto* I_cur = &I_buf1, * I_next = &I_buf2;
    auto* R_cur = &R_buf1, * R_next = &R_buf2;

    // Initial Conditions
    (*I_cur)[19] = 1e-3;
    const double& init_R = (1.0 - 1e-3) / 19.0;
    for (int j = 0; j < 19; ++j) (*R_cur)[j] = init_R;

    const std::vector<double>& K = PrecomputeKernel(L, a);

    // Time series and Slices
    std::vector<double> norm(M);
    std::vector<double> finf(M);
    std::vector<double> QI(M);
    std::vector<double> PR(M);
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

        // Spatial update loop
        for (int j = 0; j < L; ++j) {
            double Qj = 0.0, Pj = 0.0;
            for (int k = 0; k <= j; ++k) Qj += R_old[k] * K[j - k];
            for (int k = j; k < L; ++k) Pj += I_old[k] * K[k - j];

            if (t > Txy) {
                for (int k = j; k < L; ++k) Qj += R_old[k] * K[k - j];
                for (int k = 0; k <= j; ++k) Pj += I_old[k] * K[j - k];
            }

            // Dynamics and Stochasticity
            const double& newI = I_old[j] * (1.0 + dt * (R0 * Qj - 1.0));
            Inew[j] = GetStochasticCorrection(newI * N, N, gen) / N;
            Rnew[j] = R_old[j] * (1.0 - dt * R0 * Pj) + dt * I_old[j];

            // Mutations
            if (j > 0 && j < L - 1) {
                const double& in = dt * Ub * (I_old[j + 1] + I_old[j - 1]);
                const double& out = dt * Ub * 2.0 * I_old[j];
                Inew[j] += (GetStochasticCorrection(in * N, N, gen) / N) - out;
            }
        }

        // Boundary conditions & Normalization
        Inew[0] = I_old[0] + dt * Ub * (I_old[1] - I_old[0]);
        Inew[L - 1] = I_old[L - 1] + dt * Ub * (I_old[L - 2] - I_old[L - 1]);
        for (double& val : Inew) val /= norm_i;

        // Statistics at T0
        if (i == t0_step) {
            double sum_pos = 0.0, sI = 0.0;
            for (int j = 0; j < L; ++j) { sum_pos += Inew[j] * (j + 1); sI += Inew[j]; }
            meanI0 = sum_pos / sI;
        }

        // Save snapshots
        if (t > T0 && (i % show_step == 0)) {
            slice_times.push_back(t);
            I_slices.push_back(Inew);
            R_slices.push_back(Rnew);
        }

        std::swap(I_cur, I_next);
        std::swap(R_cur, R_next);
        if (i % 500 == 0) std::cout << "t = " << t << " | norm = " << norm_i << "\n";
    }

    // Post-processing & Final Results
    const auto& I_final = *I_cur;
    const double& sumI = std::accumulate(I_final.begin(), I_final.end(), 0.0);
    double meanI = 0.0;
    for (int j = 0; j < L; ++j) meanI += I_final[j] * (j + 1);
    meanI /= sumI;

    const double& speed = (meanI - meanI0) / (Tmax - T0);

    // Save everything
    SaveToFile("norm.txt", norm);
    SaveToFile("finf.txt", finf);
    SaveToFile("I_final.txt", I_final);
    SaveToFile("slice_times.txt", slice_times);

    for (size_t idx = 0; idx < I_slices.size(); ++idx) {
        SaveToFile("I_slice_" + std::to_string(idx) + ".txt", I_slices[idx]);
    }

    std::cout << "\nSimulation Complete. Speed: " << speed << "\n";
}

int main() {
    WaveSimulation(1.8, 1e-3, 7.0, 100, 1e8, 2000.0, 500.0, 500.0, 0.5);

    if (std::system("python 1D_visualize.py") != 0) {
        std::cerr << "Visualizer failed to start.\n";
    }
    return 0;
}