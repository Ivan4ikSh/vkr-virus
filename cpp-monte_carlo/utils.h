#pragma once
#include "domain.h"

namespace FILEUTILS {
    void CreateDirectory(const string& path) { filesystem::create_directories(path); }
    bool FileExists(const string& filename) { return filesystem::exists(filename); }

    bool WriteCSV(const string& filename, const vector<vector<double>>& data, const vector<string>& headers) {
        ofstream file(filename);
        if (!file.is_open()) return false;

        // «аписываем заголовки
        if (!headers.empty()) {
            for (size_t i = 0; i < headers.size(); ++i) {
                file << headers[i];
                if (i < headers.size() - 1) file << ",";
            }
            file << "\n";
        }

        // «аписываем данные
        for (const auto& row : data) {
            for (size_t i = 0; i < row.size(); ++i) {
                file << scientific << setprecision(6) << row[i];
                if (i < row.size() - 1) file << ",";
            }
            file << "\n";
        }

        file.close();
        return true;
    }

    template<typename T>
    bool WriteVector(const string& filename, const vector<T>& data) {
        ofstream file(filename);
        if (!file.is_open()) return false;

        file << scientific << setprecision(6);
        for (const auto& val : data) {
            file << val << "\n";
        }

        file.close();
        return true;
    }

    template<typename T>
    bool WriteMatrix(const string& filename, const vector<vector<T>>& matrix) {
        ofstream file(filename);
        if (!file.is_open()) return false;

        file << scientific << setprecision(6);
        for (const auto& row : matrix) {
            for (size_t i = 0; i < row.size(); ++i) {
                file << row[i];
                if (i < row.size() - 1) file << ",";
            }
            file << "\n";
        }

        file.close();
        return true;
    }

    void SaveSubplot1Data(const string& dir, const StartParams& params, const SimulationParams& sim_params) {
        // —охран€ем данные дл€ каждого временного интервала
        ofstream file(dir);
        if (!file.is_open()) return;

        file << "time,bin_center,count" << endl;
        file << scientific << setprecision(6);

        for (size_t hist_idx = 0; hist_idx < sim_params.hist_times.size(); ++hist_idx) {
            int t = sim_params.hist_times[hist_idx];
            const auto& xx = sim_params.fitness_hist_xx[hist_idx];
            const auto& nn = sim_params.fitness_hist_nn[hist_idx];
            for (size_t i = 0; i < xx.size(); ++i) {
                // Ќормализуем на s0, как в MATLAB
                if (nn[i] > 0) file << t << "," << xx[i] / params.s0 << "," << nn[i] << endl;
            }
        }
        file.close();
    }

    void SaveSubplot2Data(const string& filename, const StartParams& params, const SimulationParams& sim_params) {
        ofstream file(filename);
        if (!file.is_open()) return;

        // «аголовки столбцов
        file << "Time";
        file << ",k_av_L";           // k_av/L - красный
        file << ",sqrt_Vark_L";      // sqrt(Vark)/L - синий  
        file << ",dist_L";           // dist/L - зеленый
        file << ",C";                // C - черный пунктир
        file << ",f_survive";        // f_survive - черный
        file << ",C_all";            // C_all - пурпурный
        file << ",f1_site";          // f1_site - черный пунктир с точками
        file << "\n";

        // ‘орматирование чисел
        file << fixed << setprecision(6);

        // «аписываем данные дл€ каждого временного шага
        for (size_t t = 0; t < sim_params.T.size(); ++t) {
            file << sim_params.T[t];

            // 1. k_av/L
            double k_av_L = 0.0;
            if (t < sim_params.k_av.size() && params.L > 0) {
                k_av_L = sim_params.k_av[t] / params.L;
            }
            file << "," << k_av_L;

            // 2. sqrt(Vark)/L
            double sqrt_Vark_L = 0.0;
            if (t < sim_params.V_ark.size() && params.L > 0) {
                // Ѕерем квадратный корень из Vark (котора€ уже нормирована на s0^2)
                sqrt_Vark_L = sqrt(sim_params.V_ark[t]) / params.L;
            }
            file << "," << sqrt_Vark_L;

            // 3. dist/L (генетическое разнообразие)
            double dist_L = 0.0;
            if (t < sim_params.dist_over_L.size()) {
                dist_L = sim_params.dist_over_L[t];
            }
            file << "," << dist_L;

            // 4. C (дол€ пар с общим предком)
            double C = 0.0;
            if (t < sim_params.C.size()) {
                C = sim_params.C[t];
            }
            file << "," << C;

            // 5. f_survive (дол€ локусов с хот€ бы одним благопри€тным аллелем)
            double f_survive = 0.0;
            if (t < sim_params.f_survive.size()) {
                f_survive = sim_params.f_survive[t];
            }
            file << "," << f_survive;

            // 6. C_all (дол€ локусов, где все особи имеют общего предка)
            double C_all = 0.0;
            if (t < sim_params.C_all.size()) {
                C_all = sim_params.C_all[t];
            }
            file << "," << C_all;

            // 7. f1_site (теоретическа€ частота дл€ 1-сайтовой модели)
            double f1_site = 0.0;
            if (t < sim_params.f1_site.size()) {
                f1_site = sim_params.f1_site[t];
            }
            file << "," << f1_site << "\n";
        }

        file.close();
    }

    void SaveSubplot3Data(const string& filename, const StartParams& params, const SimulationParams& sim_params) {
        ofstream file(filename);
        if (!file.is_open()) return;

        // «аголовок: Time, locus1, locus2, ..., locusL, f1_site_theoretical
        file << "Time";
        for (int locus = 0; locus < params.L; ++locus) {
            file << ",locus_" << locus + 1;
        }
        file << ",f1_site_theoretical\n";

        file << fixed << setprecision(6);

        for (size_t t = 0; t < sim_params.T.size(); ++t) {
            file << sim_params.T[t];

            // „астоты дл€ каждого локуса
            if (t < sim_params.f_site.size()) {
                for (int locus = 0; locus < params.L; ++locus) {
                    file << "," << sim_params.f_site[t][locus];
                }
            }
            else {
                for (int locus = 0; locus < params.L; ++locus) {
                    file << "," << 0.0;
                }
            }

            // “еоретическа€ частота
            if (t < sim_params.f1_site.size()) file << "," << sim_params.f1_site[t];
            else file << "," << 0.0;
            
            file << "\n";
        }

        file.close();
    }

    void SaveSubplot4Data(const string& filename, const StartParams& params, const SimulationParams& sim_params) {
        ofstream file(filename);
        if (!file.is_open()) return;

        file << "time,bin_center,count" << endl;
        file << fixed << setprecision(6);

        // —охран€ем гистограммы частот аллелей (без нормализации на s0, так как это частоты)
        for (size_t hist_idx = 0; hist_idx < sim_params.freq_hist_times.size(); ++hist_idx) {
            int t = sim_params.freq_hist_times[hist_idx];
            const auto& xx = sim_params.freq_hist_xx[hist_idx];
            const auto& nn = sim_params.freq_hist_nn[hist_idx];

            for (size_t i = 0; i < xx.size(); ++i) {
                if (nn[i] > 0) {
                    file << t << "," << xx[i] << "," << nn[i] << endl;
                }
            }
        }

        file.close();
    }

    void SaveSimulationParams(const string& filename, const StartParams& params, const SimulationParams& sim) {
        ofstream file(filename);
        if (!file.is_open()) return;

        file << "Simulation Parameters:\n";
        file << "N = " << params.N << "\n";
        file << "L = " << params.L << "\n";
        file << "s0 = " << params.s0 << "\n";
        file << "r = " << params.r << "\n";
        file << "f0 = " << params.f0 << "\n";
        file << "muL = " << params.muL << "\n";
        file << "tf = " << params.tf << "\n";
        file << "M = " << params.M << "\n";
        file << "distribution = " << params.distribution << "\n";
        file << "exp_name = " << params.exp_name << "\n";
        file << "V_an = " << sim.V_an << "\n";
        file << "V_num = " << sim.V_num << "\n";

        file.close();
    }

    void SaveFigureData(const StartParams& params, const SimulationParams& sim_params) {
        string dir = params.output_dir + "/" + params.exp_name;

        SaveSimulationParams(dir + "/simulation_params.txt", params, sim_params);
        if (!sim_params.fitness_hist_xx.empty()) {
            SaveSubplot1Data(dir + "/fig1_221.csv", params, sim_params);
            SaveSubplot2Data(dir + "/fig1_222.csv", params, sim_params);
            SaveSubplot3Data(dir + "/fig1_223.csv", params, sim_params);
            SaveSubplot4Data(dir + "/fig1_224.csv", params, sim_params);
        }
        cout << "Data saved in directory: " << dir << endl;
    }
}
