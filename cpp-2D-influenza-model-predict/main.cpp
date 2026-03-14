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

using namespace std;
namespace fs = filesystem;

namespace CONST {
    const double EPS = 1e-16;
    const string OUTPUT = "out";
    const string DATA = "data";
}

string NAME = "test";

// ================== Параметры модели ==================
struct ModelParams {
    int seed = 42;
    int L = 25;                 // размер сетки L x L
    int M = 2000;               // число шагов по времени
    int tshow = M / 50;         // интервал сохранения состояний
    int T0 = 0;                 // первый шаг сохранения
    double R0 = 2.5;            // базовое репродуктивное число
    double Ub = 1e-3;           // частота мутаций
    double a = 7.0;             // масштаб кросс-иммунитета
    double dt = 0.5;            // шаг по времени (в единицах t_rec)
    double init_infected = 1e-3; // начальная доля инфицированных
    int64_t N = 1e8;            // общая численность популяции (не используется напрямую)
    double init_radius = 2.0;   // радиус начального очага
    double mutation_step = 1.0; // не используется (оставлено для совместимости)
};

class EpidemicSimulator {
private:
    ModelParams p_;
    mt19937 rng_;
    uniform_real_distribution<double> uniform_{ 0.0, 1.0 };
    normal_distribution<double> normal_{ 0.0, 1.0 };

    vector<vector<double>> I_; // инфицированные штаммом (x,y)
    vector<vector<double>> R_; // выздоровевшие от штамма (x,y)

    vector<double> norm_;
    vector<double> finf_;
    vector<double> frec_;
    vector<double> total_I_;
    vector<double> wave_center_;
    vector<double> wave_r_;      // расстояние центра масс от геометрического центра
    vector<double> max_I_;       // максимальная доля инфицированных
    vector<double> diversity_;   // число активных штаммов

public:
    EpidemicSimulator(const ModelParams& params, int seed = time(NULL)) : p_(params), rng_(seed) {
        InitializeState();
        PrintParameters();

        norm_.resize(p_.M, 0.0);
        finf_.resize(p_.M, 0.0);
        frec_.resize(p_.M, 0.0);
        total_I_.resize(p_.M, 0.0);
        wave_center_.resize(p_.M, 0.0);
        wave_r_.resize(p_.M, 0.0);
        max_I_.resize(p_.M, 0.0);
        diversity_.resize(p_.M, 0.0);
    }

    void Run() {
        PrintTableHeader();
        for (int step = 0; step < p_.M; ++step) {
            if (step >= p_.T0 && step % p_.tshow == 0) SaveState(step);
            CalculateStatistics(step);
            Step();
            if (step % max(1, p_.M / 10) == 0) PrintProgress(step);
        }
        SaveFinalResults();
        cout << "\n" << string(60, '=') << "\n";
    }

private:
    // Кросс-иммунитет: K(d) = 1 - exp(-d / a)
    double K(double dist) const {
        return 1.0 - exp(-dist / p_.a);
    }

    void InitializeState() {
        int L = p_.L;
        I_.assign(L, vector<double>(L, 0.0));
        R_.assign(L, vector<double>(L, 0.0));

        int cx = L / 2;          // центр начального очага инфекции
        int cy = L / 2;
        double r = p_.init_radius;

        // ---- Задаём начальное распределение инфицированных I ----
        vector<pair<int, int>> infected_cells;
        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                double dist = sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
                if (dist <= r) infected_cells.emplace_back(x, y);
            }
        }

        double I_total = p_.init_infected;   // желаемая суммарная доля инфицированных
        int cloud_size = infected_cells.size();
        for (auto& [x, y] : infected_cells) {
            I_[x][y] = 1.0 / cloud_size;
        }

        // Проверка сумм
        double sumI = 0.0;
        double sumR = 0.0;
        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                sumI += I_[x][y];
                sumR += R_[x][y];
            }
        }
        cout << "Initialization: I = " << sumI << ", R = " << sumR << ", total = " << sumI + sumR << "\n";
    }

    // Один шаг по времени (явный Эйлер)
    void Step() {
        int L = p_.L;
        double dt = p_.dt;
        double R0 = p_.R0;

        // Вычисляем свёртки с ядром K
        vector<vector<double>> sumR(L, vector<double>(L, 0.0));
        vector<vector<double>> sumI(L, vector<double>(L, 0.0));

        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                double sR = 0.0, sI = 0.0;
                for (int x1 = 0; x1 < L; ++x1) {
                    for (int y1 = 0; y1 < L; ++y1) {
                        // периодические границы
                        int dx = abs(x - x1);
                        int dy = abs(y - y1);
                        dx = min(dx, L - dx);
                        dy = min(dy, L - dy);
                        double dist = sqrt(dx * dx + dy * dy);
                        double kval = K(dist);
                        sR += kval * R_[x1][y1];
                        sI += kval * I_[x1][y1];
                    }
                }
                sumR[x][y] = sR;
                sumI[x][y] = sI;
            }
        }

        // Сохраняем старые значения для согласованного обновления
        auto I_old = I_;
        auto R_old = R_;

        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                double I = I_old[x][y];
                double R = R_old[x][y];

                // dI/dt = I*(R0*sumR - 1) — нетто рост/спад инфицированных
                double growth = I * (R0 * sumR[x][y] - 1.0) * dt;
                I_[x][y] = I + growth;

                // dR/dt = I - R*R0*sumI — приход выздоровевших минус реинфекция
                double recovery = I * dt;               // выздоровевшие из I (используем I_old!)
                double reinfection = R * R0 * sumI[x][y] * dt; // реинфицированные
                R_[x][y] = R + recovery - reinfection;
            }
        }
        // Мутации (диффузия в пространстве штаммов)
        ApplyMutation();
    }

    void ApplyMutation() {
        int L = p_.L;
        double dt = p_.dt;
        double Ub = p_.Ub;
        vector<vector<double>> newI = I_;
        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                double out = 4 * Ub * I_[x][y];
                double in = 0.0;

                if (x > 0) in += I_[x - 1][y];
                if (x < L - 1) in += I_[x + 1][y];
                if (y > 0) in += I_[x][y - 1];
                if (y < L - 1)in += I_[x][y + 1];

                in *= Ub;
                newI[x][y] += dt * (in - out);
                if (newI[x][y] < 0) newI[x][y] = 0.0;
            }
        }
        I_.swap(newI);

        // Небольшой стохастический шум (можно убрать, если не нужен)
        double noise = 1e-6;
        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                double n = 1.0 + noise * normal_(rng_);
                I_[x][y] *= max(0.0, n);
            }
        }
    }

    void CalculateStatistics(int step) {
        double sumI = 0.0;
        double sumR = 0.0;
        double weighted_x = 0.0;
        double weighted_y = 0.0;
        double max_val = 0.0;
        int div_count = 0;
        double threshold = 1e-8;
        for (int x = 0; x < p_.L; ++x) {
            for (int y = 0; y < p_.L; ++y) {
                double i_val = I_[x][y];
                sumI += i_val;
                sumR += R_[x][y];
                weighted_x += x * i_val;
                weighted_y += y * i_val;

                if (i_val > max_val) max_val = i_val;
                if (i_val > threshold) ++div_count;
            }
        }

        double total = sumI + sumR;
        norm_[step] = total;
        finf_[step] = (total > 0.0) ? sumI / total : 0.0;
        frec_[step] = (total > 0.0) ? sumR / total : 0.0;
        total_I_[step] = sumI;
        max_I_[step] = max_val;
        diversity_[step] = static_cast<double>(div_count);

        if (sumI > 0) {
            double wx = weighted_x / sumI;
            double wy = weighted_y / sumI;
            double cx = p_.L / 2.0;
            double cy = p_.L / 2.0;
            wave_r_[step] = sqrt((wx - cx) * (wx - cx) + (wy - cy) * (wy - cy));
            wave_center_[step] = sqrt(wx * wx + wy * wy); // от (0,0), для совместимости
        }
        else {
            wave_r_[step] = 0.0;
            wave_center_[step] = 0.0;
        }
    }

    void SaveState(int step) {
        const string& output_dir = CONST::OUTPUT + "/" + NAME;
        const string& data_dir = output_dir + "/" + CONST::DATA;
        fs::create_directories(data_dir);

        auto save_matrix = [&](const vector<vector<double>>& mat, const string& fname) {
            ofstream fout(data_dir + "/" + fname);
            for (int x = 0; x < p_.L; ++x) {
                for (int y = 0; y < p_.L; ++y) {
                    fout << mat[x][y];
                    if (y < p_.L - 1) fout << ";";
                }
                fout << "\n";
            }
            };

        save_matrix(I_, "I_step_" + to_string(step) + ".csv");
        save_matrix(R_, "R_step_" + to_string(step) + ".csv");
    }

    void SaveFinalResults() {
        const string& output_dir = CONST::OUTPUT + "/" + NAME;
        const string& data_dir = output_dir + "/" + CONST::DATA;
        fs::create_directories(data_dir);

        ofstream param_file(data_dir + "/parameters.txt");
        param_file << "L = " << p_.L << "\n";
        param_file << "M = " << p_.M << "\n";
        param_file << "R0 = " << p_.R0 << "\n";
        param_file << "Ub = " << p_.Ub << "\n";
        param_file << "a = " << p_.a << "\n";
        param_file << "dt = " << p_.dt << "\n";
        param_file << "init_infected = " << p_.init_infected << "\n";
        param_file << "N = " << p_.N << "\n";
        param_file << "Seed = " << p_.seed << "\n";
        param_file.close();
    }

    void PrintTableHeader() {
        cout << left << setw(5) << "step" << " | "
            << left << setw(8) << "f_inf" << " | "
            << left << setw(8) << "f_rec" << " | "
            << left << setw(6) << "norm" << " | "
            << left << setw(7) << "wave_r" << " | "
            << left << setw(10) << "max_I" << " | "
            << left << setw(9) << "diversity" << "\n";
        cout << string(70, '-') << "\n";
    }

    void PrintProgress(int step) {
        cout << left << setw(5) << step << " | "
            << fixed << setprecision(6) << setw(8) << finf_[step] << " | "
            << setprecision(6) << setw(8) << frec_[step] << " | "
            << setprecision(2) << setw(6) << norm_[step] << " | "
            << setprecision(3) << setw(7) << wave_r_[step] << " | "
            << scientific << setprecision(2) << setw(10) << max_I_[step] << " | "
            << fixed << setw(9) << static_cast<int>(diversity_[step]) << "\n";
    }

    void PrintParameters() {
        cout << "Grid size: " << p_.L << " x " << p_.L << "\n";
        cout << "Time steps: " << p_.M << "\n";
        cout << "R0 = " << p_.R0 << "\n";
        cout << "Ub = " << p_.Ub << "\n";
        cout << "a = " << p_.a << "\n";
        cout << "dt = " << p_.dt << "\n";
        cout << "init_infected = " << p_.init_infected << "\n";
        cout << "N = " << p_.N << "\n";
        cout << "Seed = " << p_.seed << "\n\n";
    }
};

void InitDirectory() {
    const string& output_dir = CONST::OUTPUT + "/" + NAME;
    const string& data_dir = output_dir + "/" + CONST::DATA;
    fs::create_directories(data_dir);
    for (const auto& entry : fs::directory_iterator(data_dir)) fs::remove_all(entry.path());
}

void TEST(ModelParams& params, const string& name) {
    NAME = name;
    ++params.M;
    InitDirectory();
    EpidemicSimulator sim(params, params.seed);
    sim.Run();

    const string& output_dir = CONST::OUTPUT + "/" + NAME;
    const string& data_dir = output_dir + "/" + CONST::DATA;
    system(("python 2D_visualize.py -d " + data_dir + " -s " + output_dir + "/" + NAME + ".png -a " + output_dir + "/" + NAME + ".gif").c_str());
}

int main() {
    ModelParams params;
    params.L = 30;
    params.M = 1000;
    params.tshow = params.M / 50;
    params.T0 = 0;
    params.R0 = 1.8;
    params.Ub = 1e-4;
    params.a = 0.5;
    params.dt = 1.0 / (params.R0 * 5);   // ~0.111
    params.N = 1e5;
    params.init_infected = double(params.L * params.L) / params.N;
    params.seed = 2;
    params.init_radius = 1.0;

    TEST(params, "test");
    return 0;
}