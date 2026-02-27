#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

namespace CONST {
    const string VISUALIZE_PYTHON = "visualize.py";
    const string OUTPUT_DIR_NAME = "output";
    const string DATA_DIR = "data";
    const double EPS = 1e-8;
    const double PI = 3.1415926535897932;
    const double STOCH_THRESHOLD = 10.0;
}

string OUTPUT_NAME = "output";
double RADIUS = 2.0;

struct Point {
    int x;
    int y;

    bool operator==(const Point& other) const {
        return x == other.x && y == other.y;
    }
};

struct PointHash {
    size_t operator()(const Point& p) const {
        return hash<int>()(p.x) ^ (hash<int>()(p.y) << 1);
    }
};

struct ModelParameters {
    int L;                // размер решётки LxL
    int M;                // число шагов по времени
    int tshow;            // интервал отображения
    int T0;               // начало отображения
    double beta;
    double R0;            // базовый репродуктивный номер
    double Ub;            // частота мутаций на геном за шаг
    double a;             // полурасстояние иммунитета
    double Varx;          // разброс координат
    double asym;          // коэффициент асимметрии
    double dt;            // шаг по времени
    double init_infected; // начальная доля инфицированных
    int64_t N;            // общая численность популяции
};

class EpidemicSimulator {
private:
    ModelParameters params_;
    // Доли инфицированных и переболевших (сумма всех I+R = 1)
    vector<vector<double>> I_;
    vector<vector<double>> R_;
    vector<vector<double>> X_;   // случайные x-координаты каждой клетки
    vector<vector<double>> Y_;   // случайные y-координаты каждой клетки
    // Предвычисленное ядро кросс-иммунитета: K_matrix_[i][j][k][l] = K(расстояние между (i,j) и (k,l))
    vector<vector<vector<vector<double>>>> K_matrix_;

    vector<double> norm_;
    vector<double> finf_;

    double Txy_;            // не используется в 2D

    mt19937 generator_;
    uniform_real_distribution<double> uniform_dist_;

public:
    EpidemicSimulator(ModelParameters params, int seed = time(nullptr))
        : params_(params), generator_(seed), uniform_dist_(0.0, 1.0) {
        GenerateCoordinates();
        PrecomputeKernel();
        InitializeState();
    }

    void Run() {
        cout << "========================================\n";
        cout << "2D MODEL OF ANTIGENIC EVOLUTION WITH COMETS\n";
        cout << "========================================\n\n";

        PrintParameters();
        cout << "\nStarting simulation...\n";

        auto start_time = chrono::high_resolution_clock::now();

        for (int step = 0; step < params_.M; ++step) {
            if (step >= params_.T0 && step % params_.tshow == 0) SaveCurrentState(step);
            CalculateStatistics(step);
            StepSimulation();
            if (params_.M >= 10 && step % (params_.M / 10) == 0) PrintStepInfo(step);
        }

        auto end_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);

        cout << "\nSimulation completed in " << duration.count() / 1000.0 << " seconds\n";

        SaveFinalResults();
        PrintFinalStatistics();
    }

private:
    void GenerateCoordinates() {
        int L = params_.L;
        X_.resize(L, vector<double>(L));
        Y_.resize(L, vector<double>(L));

        uniform_real_distribution<> dis(-params_.Varx, params_.Varx);

        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                X_[i][j] = j + dis(generator_);
                Y_[i][j] = i + dis(generator_);
            }
        }
    }

    void PrecomputeKernel() {
        int L = params_.L;
        double a = params_.a;
        double asym = params_.asym;

        K_matrix_.resize(L, vector<vector<vector<double>>>(L, vector<vector<double>>(L, vector<double>(L, 0.0))));

        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                for (int k = 0; k < L; ++k) {
                    for (int l = 0; l < L; ++l) {
                        double dx = X_[k][l] - X_[i][j];
                        double dy = Y_[k][l] - Y_[i][j];
                        double dist = sqrt(dx * dx + asym * dy * dy) / a;
                        K_matrix_[i][j][k][l] = dist / (1.0 + dist);
                    }
                }
            }
        }
    }

    void InitializeState() {
        int L = params_.L;
        int center_x = L / 2;
        int center_y = L / 2;
        int64_t N_total = params_.N;

        // Равномерное распределение населения по клеткам
        int64_t people_per_cell = N_total / (L * L);
        int64_t remainder = N_total % (L * L);

        // Определяем, какие клетки попадают в начальный круг
        vector<vector<bool>> in_circle(L, vector<bool>(L, false));
        int cells_in_circle = 0;
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                double dx = j - center_x;
                double dy = i - center_y;
                if (dx * dx + dy * dy <= RADIUS * RADIUS) {
                    in_circle[i][j] = true;
                    ++cells_in_circle;
                }
            }
        }

        // Общее число инфицированных
        int64_t infected_total = static_cast<int64_t>(params_.init_infected * N_total);
        int64_t infected_per_cell = cells_in_circle ? infected_total / cells_in_circle : 0;
        int64_t infected_rem = cells_in_circle ? infected_total % cells_in_circle : 0;

        // Инициализация долей
        I_.resize(L, vector<double>(L, 0.0));
        R_.resize(L, vector<double>(L, 0.0));

        int infected_idx = 0;
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                int64_t tot = people_per_cell + (i * L + j < remainder ? 1 : 0); // население клетки
                int64_t I_val = 0;
                if (in_circle[i][j] && infected_idx < cells_in_circle) {
                    I_val = infected_per_cell + (infected_idx < infected_rem ? 1 : 0);
                    ++infected_idx;
                }
                int64_t R_val = tot - I_val;
                // Доли от общей популяции
                I_[i][j] = static_cast<double>(I_val) / N_total;
                R_[i][j] = static_cast<double>(R_val) / N_total;
            }
        }

        norm_.resize(params_.M, 0.0);
        finf_.resize(params_.M, 0.0);
    }

    void StepSimulation() {
        int L = params_.L;
        double dt = params_.dt;
        double R0 = params_.R0;
        double Ub = params_.Ub;
        int64_t N_total = params_.N;

        // Сохраняем старые состояния (доли)
        vector<vector<double>> I_old = I_;
        vector<vector<double>> R_old = R_;

        // Вычисление Q и P (свёртки с ядром) — доли, поэтому Q и P ≤ 1
        vector<vector<double>> Q(L, vector<double>(L, 0.0));
        vector<vector<double>> P(L, vector<double>(L, 0.0));

        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                double sumQ = 0.0, sumP = 0.0;
                for (int k = 0; k < L; ++k) {
                    for (int l = 0; l < L; ++l) {
                        double K = K_matrix_[i][j][k][l];
                        sumQ += R_old[k][l] * K;
                        sumP += I_old[k][l] * K;
                    }
                }
                Q[i][j] = sumQ;
                P[i][j] = sumP;
            }
        }

        // Временные массивы для новых долей
        vector<vector<double>> I_new(L, vector<double>(L, 0.0));
        vector<vector<double>> R_new(L, vector<double>(L, 0.0));

        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                // --- Новое I без мутации (доля) ---
                double newI = I_old[i][j] * (1.0 + dt * (R0 * Q[i][j] - 1.0));

                // --- Стохастическая коррекция для newI (работаем с абсолютными числами) ---
                double lambda = newI * N_total;          // ожидаемое абсолютное число
                if (lambda < CONST::STOCH_THRESHOLD) {
                    int xm = static_cast<int>(round(6.0 * lambda));
                    if (xm > 0) {
                        double prob = lambda / xm;
                        int count = 0;
                        for (int n = 0; n < xm; ++n) {
                            if (uniform_dist_(generator_) < prob) ++count;
                        }
                        newI = static_cast<double>(count) / N_total;  // доля
                    }
                    else {
                        newI = 0.0;
                    }
                }
                I_new[i][j] = newI;

                // --- Обновление R (доля) ---
                double newR = R_old[i][j] * (1.0 - dt * R0 * P[i][j]) + dt * I_old[i][j];
                if (newR < 0) newR = 0.0; // защита от отрицательных
                R_new[i][j] = newR;
            }
        }

        // --- Мутационный член (диффузия) — добавляется к I_new ---
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                double in = 0.0;
                int neighbor_count = 0;
                if (i > 0) { in += I_old[i - 1][j]; neighbor_count++; }
                if (i < L - 1) { in += I_old[i + 1][j]; neighbor_count++; }
                if (j > 0) { in += I_old[i][j - 1]; neighbor_count++; }
                if (j < L - 1) { in += I_old[i][j + 1]; neighbor_count++; }

                double in_mut = dt * Ub * in;          // доля
                double out_mut = dt * Ub * neighbor_count * I_old[i][j]; // доля

                // Стохастика для входного потока (абсолютное число)
                double lambda_in = in_mut * N_total;
                if (lambda_in < CONST::STOCH_THRESHOLD) {
                    int xm = static_cast<int>(round(6.0 * lambda_in));
                    if (xm > 0) {
                        double prob = lambda_in / xm;
                        int count = 0;
                        for (int n = 0; n < xm; ++n) {
                            if (uniform_dist_(generator_) < prob) ++count;
                        }
                        in_mut = static_cast<double>(count) / N_total; // доля
                    }
                    else {
                        in_mut = 0.0;
                    }
                }

                I_new[i][j] += in_mut - out_mut;
                if (I_new[i][j] < 0) I_new[i][j] = 0.0;
            }
        }

        // --- Нормировка: сохраняем общую сумму I+R = 1 ---
        double sum_total = 0.0;
        for (int i = 0; i < L; ++i)
            for (int j = 0; j < L; ++j)
                sum_total += I_new[i][j] + R_new[i][j];

        if (sum_total > 0) {
            for (int i = 0; i < L; ++i) {
                for (int j = 0; j < L; ++j) {
                    I_[i][j] = I_new[i][j] / sum_total;
                    R_[i][j] = R_new[i][j] / sum_total;
                }
            }
        }
        else {
            // В极端 случае всё население вымерло (не должно происходить)
            cerr << "Warning: total population zero at step!" << endl;
            I_ = I_new; // сохраняем как есть (всё равно нули)
            R_ = R_new;
        }
    }

    void CalculateStatistics(int step) {
        double total_I = 0.0;
        double total_R = 0.0;

        for (int i = 0; i < params_.L; ++i) {
            for (int j = 0; j < params_.L; ++j) {
                total_I += I_[i][j];
                total_R += R_[i][j];
            }
        }

        double total = total_I + total_R;
        norm_[step] = total;                // должно быть близко к 1
        finf_[step] = (total > 0) ? total_I / total : 0.0;
    }

    int CountInfectedCells() const {
        int count = 0;
        for (int i = 0; i < params_.L; ++i) {
            for (int j = 0; j < params_.L; ++j) {
                if (I_[i][j] > CONST::EPS) ++count;
            }
        }
        return count;
    }

    void SaveCurrentState(int step) {
        int L = params_.L;
        // Сохраняем доли
        string dir_path = CONST::OUTPUT_DIR_NAME + "/" + OUTPUT_NAME + "/" + CONST::DATA_DIR;
        filesystem::create_directories(dir_path);

        ofstream I_file(dir_path + "/state_I_step_" + to_string(step) + ".csv");
        ofstream R_file(dir_path + "/state_R_step_" + to_string(step) + ".csv");

        if (I_file.is_open() && R_file.is_open()) {
            for (int i = 0; i < L; ++i) {
                for (int j = 0; j < L; ++j) {
                    I_file << fixed << setprecision(16) << I_[i][j];
                    R_file << fixed << setprecision(16) << R_[i][j];
                    if (j < L - 1) {
                        I_file << ";";
                        R_file << ";";
                    }
                }
                I_file << "\n";
                R_file << "\n";
            }
        }
    }

    void SaveFinalResults() {
        string dir_path = CONST::OUTPUT_DIR_NAME + "/" + OUTPUT_NAME + "/" + CONST::DATA_DIR;

        ofstream norm_file(dir_path + "/norm_time_series.csv");
        ofstream finf_file(dir_path + "/finf_time_series.csv");

        if (norm_file.is_open() && finf_file.is_open()) {
            norm_file << "Step;Norm\n";
            finf_file << "Step;Finf\n";

            for (int i = 0; i < params_.M; ++i) {
                norm_file << i << ";" << norm_[i] << "\n";
                finf_file << i << ";" << finf_[i] << "\n";
            }
        }
        norm_file.close();
        finf_file.close();

        ofstream param_file(dir_path + "/parameters.txt");
        if (param_file.is_open()) {
            param_file << "MODEL PARAMETERS:\n";
            param_file << "L: " << params_.L << "\n";
            param_file << "M: " << params_.M << "\n";
            param_file << "R0: " << params_.R0 << "\n";
            param_file << "Ub: " << params_.Ub << "\n";
            param_file << "a: " << params_.a << "\n";
            param_file << "Varx: " << params_.Varx << "\n";
            param_file << "asym: " << params_.asym << "\n";
            param_file << "dt: " << params_.dt << "\n";
            param_file << "init_infected: " << params_.init_infected << "\n";
            param_file << "N: " << params_.N << "\n";
        }
        param_file.close();
    }

    void PrintStepInfo(const int& step) {
        cout << "Step " << step << "/" << params_.M << " (" << (100 * step / params_.M) << "%)" << " - Infected cells: " << CountInfectedCells() << " | ";
        cout << "norm = " << norm_[step] << ", finf = " << finf_[step] << "\n";
    }

    void PrintParameters() {
        cout << "Model Parameters:\n";
        cout << "  Grid size: " << params_.L << "x" << params_.L << "\n";
        cout << "  Time steps: " << params_.M << "\n";
        cout << "  R0: " << params_.R0 << "\n";
        cout << "  Mutation rate (Ub): " << params_.Ub << "\n";
        cout << "  Cross-immunity distance (a): " << params_.a << "\n";
        cout << "  Coordinate variance (Varx): " << params_.Varx << "\n";
        cout << "  Asymmetry: " << params_.asym << "\n";
        cout << "  Time step (dt): " << params_.dt << "\n";
        cout << "  Initial infected fraction: " << params_.init_infected << "\n";
        cout << "  Total population size: " << params_.N << "\n";
    }

    void PrintFinalStatistics() {
        double total_I = 0.0;
        double total_R = 0.0;
        int infected_cells = 0;

        for (int i = 0; i < params_.L; ++i) {
            for (int j = 0; j < params_.L; ++j) {
                total_I += I_[i][j];
                total_R += R_[i][j];
                if (I_[i][j] > CONST::EPS) ++infected_cells;
            }
        }

        cout << "\nFinal Statistics:\n";
        cout << "  Total infected fraction: " << total_I << "\n";
        cout << "  Total recovered fraction: " << total_R << "\n";
        cout << "  Total fraction: " << total_I + total_R << "\n";
        cout << "  Infected cells: " << infected_cells << "\n";
    }
};

void InitDirectory(const string& data_dir) {
    filesystem::create_directories(data_dir + "/" + CONST::DATA_DIR);
    for (const auto& entry : filesystem::directory_iterator(data_dir + "/" + CONST::DATA_DIR)) {
        filesystem::remove_all(entry.path());
    }
}

void RunExperiment(const string& name) {
    cout << "========================================\n";
    cout << "Experiment: " << name << endl;
    cout << "========================================\n\n";

    ModelParameters params;
    params.L = 50;
    params.M = 800;
    params.tshow = params.M / 50;
    params.T0 = 0;
    params.R0 = 2.5;
    params.Ub = 1e-3;
    params.a = 7.0;
    params.Varx = 0.1;
    params.asym = 1.0;
    params.dt = 0.5;
    params.init_infected = 0.01;
    params.N = 1e8;
    params.beta = 1.0;

    OUTPUT_NAME = name;
    RADIUS = 2.0;

    string data_dir = CONST::OUTPUT_DIR_NAME + "/" + OUTPUT_NAME;

    InitDirectory(data_dir);
    EpidemicSimulator simulator(params, 42);
    simulator.Run();

    // Вызов Python для визуализации (раскомментировать при наличии скрипта)
    string cmd = "python " + CONST::VISUALIZE_PYTHON +
        " --data-dir " + data_dir + "/" + CONST::DATA_DIR +
        " --animation " + data_dir + "/" + OUTPUT_NAME + ".gif" +
        " --snapshot " + data_dir + "/" + OUTPUT_NAME + ".png";
    system(cmd.c_str());
    system(("python stats.py --data-dir " + data_dir + "/" + CONST::DATA_DIR + " --output-dir " + data_dir).c_str());
}

int main() {
    RunExperiment("test");
    return 0;
}