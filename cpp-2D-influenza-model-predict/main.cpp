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

struct Strain {
    int64_t I; // абсолютное число инфицированных
    int64_t R; // абсолютное число восприимчивых
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
    double R0;            // базовый репродуктивный номер
    double Ub;            // частота мутаций на геном за шаг
    double a;             // полурасстояние иммунитета
    double Varx;          // разброс координат
    double asym;          // коэффициент асимметрии
    double dt;            // шаг по времени
    double init_infected; // начальная доля инфицированных
    int64_t N_total;      // общая численность популяции
};

class EpidemicSimulator {
private:
    ModelParameters params_;
    unordered_map<Point, Strain, PointHash> strains_;
    vector<vector<double>> X_;
    vector<vector<double>> Y_;
    vector<double> norm_;
    vector<double> finf_;
    mt19937 generator_;

    uniform_real_distribution<double> uniform_dist_;

public:
    EpidemicSimulator(ModelParameters params, int seed = time(nullptr)) : params_(params), generator_(seed), uniform_dist_(0.0, 1.0) {
        GenerateCoordinates();
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
            if (params_.M >= 10 && step % (params_.M / 10) == 0)
                cout << "Step " << step << "/" << params_.M << " (" << (100 * step / params_.M) << "%)"
                << " - Infected cells: " << CountInfectedCells() << endl;
            StepSimulation();
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

    void InitializeState() {
        int L = params_.L;
        int center_x = L / 2;
        int center_y = L / 2;

        // Общее число клеток
        int64_t total_cells = L * L;
        int64_t people_per_cell = params_.N_total / total_cells;
        int64_t remainder = params_.N_total % total_cells;

        // Равномерное распределение населения по клеткам
        vector<vector<int64_t>> population(L, vector<int64_t>(L, people_per_cell));
        for (int i = 0; i < L && remainder > 0; ++i) {
            for (int j = 0; j < L && remainder > 0; ++j) {
                ++population[i][j];
                --remainder;
            }
        }

        // Подсчёт клеток внутри начального круга
        int cells_in_circle = 0;
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                double dx = j - center_x;
                double dy = i - center_y;
                if (dx * dx + dy * dy <= RADIUS * RADIUS) ++cells_in_circle;
            }
        }

        int64_t infected_total = static_cast<int64_t>(params_.init_infected * params_.N_total);
        int64_t infected_per_cell = cells_in_circle ? infected_total / cells_in_circle : 0;
        int64_t infected_rem = cells_in_circle ? infected_total % cells_in_circle : 0;

        strains_.clear();
        int infected_idx = 0;
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                int64_t tot = population[i][j];
                int64_t I = 0;
                double dx = j - center_x;
                double dy = i - center_y;
                if (dx * dx + dy * dy <= RADIUS * RADIUS && infected_idx < cells_in_circle) {
                    I = infected_per_cell + (infected_idx < infected_rem ? 1 : 0);
                    ++infected_idx;
                }
                int64_t S = tot - I;  // остальные восприимчивы
                strains_[{i, j}] = { I, S };
            }
        }

        norm_.resize(params_.M, 0.0);
        finf_.resize(params_.M, 0.0);
    }

    // Асимметричная кросс-иммунитетная функция: заражение возможно только если антигенная координата цели больше, чем источника
    double KFunc(int i1, int j1, int i2, int j2) {
        double dx = X_[i2][j2] - X_[i1][j1];
        //if (dx >= 0) return 0.0;  // запрет заражения «назад»
        double dy = Y_[i2][j2] - Y_[i1][j1];
        //if (dy >= 0) return 0.0;  // запрет заражения «назад»
        double dist = sqrt(dx * dx + params_.asym * dy * dy);
        dist /= params_.a;
        return dist / (1.0 + dist);
    }

    void StepSimulation() {
        int L = params_.L;
        double dt = params_.dt;
        double R0 = params_.R0;
        double Ub = params_.Ub;
        int64_t Ntot = params_.N_total;

        unordered_map<Point, Strain, PointHash> new_strains = strains_; // Новое состояние после детерминированного шага (копируем текущее)
        unordered_map<Point, int64_t, PointHash> mutation_flow; // Поток мутаций (изменение абсолютного числа инфицированных в клетках)

        // Детерминированное обновление I и R и подготовка мутаций
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                Point p{ i, j };

                auto it = strains_.find(p);
                if (it == strains_.end()) continue;

                int64_t I_abs = it->second.I;
                int64_t S_abs = it->second.R;

                if (I_abs == 0 && S_abs == 0) continue;

                // Переводим в double для вычислений
                double I_frac = static_cast<double>(I_abs) / Ntot;
                double S_frac = static_cast<double>(S_abs) / Ntot;

                // Вычисляем Q (сумма по восприимчивым, взвешенная на K) и P (сумма по инфицированным)
                double Q = 0.0, P = 0.0;
                for (const auto& [point, strain] : strains_) {
                    int m = point.x;
                    int n = point.y;
                    Q += static_cast<double>(strain.R) / Ntot * KFunc(i, j, m, n); // Q: влияние восприимчивых на заражение от (i,j)
                    P += static_cast<double>(strain.I) / Ntot * KFunc(m, n, i, j); // P: влияние инфицированных на заражение (i,j)
                }

                double dI_dt = I_frac * (R0 * Q - 1.0);
                double dS_dt = I_frac - R0 * S_frac * P;

                double new_I_frac = I_frac + dt * dI_dt;
                double new_S_frac = S_frac + dt * dS_dt;

                // Предотвращаем отрицательные значения
                new_I_frac = max(0.0, new_I_frac);
                new_S_frac = max(0.0, new_S_frac);

                // Переводим обратно в абсолютные числа
                int64_t new_I_abs = llround(new_I_frac * Ntot);
                int64_t new_S_abs = llround(new_S_frac * Ntot);

                // Сохраняем детерминированные значения
                new_strains[p] = { new_I_abs, new_S_abs };
                if (new_I_abs <= 0) continue;

                double mean_mut = Ub * new_I_abs;
                if (mean_mut <= 0.0) continue;

                // Генерируем пуассоновское число мутантов с оптимизацией
                int64_t n_mut;
                if (mean_mut > 10.0) n_mut = static_cast<int64_t>(std::round(mean_mut));
                else if (mean_mut >= 0.01) {
                    poisson_distribution<int64_t> poiss(mean_mut);
                    n_mut = poiss(generator_);
                }
                else n_mut = 0;
                if (n_mut == 0) continue;

                // Распределяем мутантов по 4 соседним клеткам
                const int dirs[4][2] = { {1,0}, {-1,0}, {0,1}, {0,-1} };
                for (int m = 0; m < n_mut; ++m) {
                    int dir = static_cast<int>(uniform_dist_(generator_) * 4); // случайное направление 0..3
                    int ni = i + dirs[dir][0];
                    int nj = j + dirs[dir][1];

                    mutation_flow[p] -= 1; // Мутант покидает исходную клетку
                    if (ni >= 0 && ni < L && nj >= 0 && nj < L) mutation_flow[{ni, nj}] += 1; // Мутант попадает в целевую клетку
                    // иначе мутант теряется (за границей)
                }
            }
        }

        // Применяем мутационный поток к new_strains
        for (const auto& [point, delta_abs] : mutation_flow) {
            if (delta_abs == 0) continue;
            new_strains[point].I += delta_abs;
            if (new_strains[point].I < 0) new_strains[point].I = 0; // защита от отрицательных
        }

        // Обновляем глобальное состояние
        strains_ = move(new_strains);
    }

    void CalculateStatistics(int step) {
        double total_I = 0.0;
        double total_R = 0.0;

        for (const auto& [point, strain] : strains_) {
            total_I += strain.I;
            total_R += strain.R;
        }

        double total = total_I + total_R;
        norm_[step] = total / params_.N_total; // нормированная полная численность
        finf_[step] = (total > 0) ? total_I / total : 0.0;
    }

    int CountInfectedCells() const {
        int count = 0;
        for (const auto& [point, strain] : strains_) {
            if (strain.I > 0) ++count;
        }
        return count;
    }

    void SaveCurrentState(int step) {
        int L = params_.L;
        vector<vector<double>> I_dense(L, vector<double>(L, 0.0));
        vector<vector<double>> R_dense(L, vector<double>(L, 0.0));

        for (const auto& [point, strain] : strains_) {
            int i = point.x;
            int j = point.y;
            if (i >= 0 && i < L && j >= 0 && j < L) {
                I_dense[i][j] = static_cast<double>(strain.I) / params_.N_total;
                R_dense[i][j] = static_cast<double>(strain.R) / params_.N_total;
            }
        }

        string dir_path = CONST::OUTPUT_DIR_NAME + "/" + OUTPUT_NAME + "/" + CONST::DATA_DIR;
        filesystem::create_directories(dir_path);

        ofstream I_file(dir_path + "/state_I_step_" + to_string(step) + ".csv");
        ofstream R_file(dir_path + "/state_R_step_" + to_string(step) + ".csv");

        if (I_file.is_open() && R_file.is_open()) {
            for (int i = 0; i < L; ++i) {
                for (int j = 0; j < L; ++j) {
                    I_file << fixed << setprecision(16) << I_dense[i][j];
                    R_file << fixed << setprecision(16) << R_dense[i][j];
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
            param_file << "N_total: " << params_.N_total << "\n";
        }
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
        cout << "  Total population size: " << params_.N_total << "\n";
    }

    void PrintFinalStatistics() {
        double total_I = 0.0;
        double total_R = 0.0;
        int infected_cells = 0;
        int susceptible_cells = 0;

        for (const auto& [point, strain] : strains_) {
            total_I += strain.I;
            total_R += strain.R;

            if (strain.I > 0) ++infected_cells;
            if (strain.R > 0) ++susceptible_cells;
        }

        cout << "\nFinal Statistics:\n";
        cout << "  Total infected fraction: " << total_I / params_.N_total << "\n";
        cout << "  Total susceptible fraction: " << total_R / params_.N_total << "\n";
        cout << "  Total fraction: " << (total_I + total_R) / params_.N_total << "\n";
        cout << "  Infected cells: " << infected_cells << "\n";
        cout << "  Susceptible cells: " << susceptible_cells << "\n";
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
    params.L = 40;
    params.M = 2'000; // увеличено для достаточного времени развития волны
    params.tshow = params.M / 50;
    params.T0 = 0;
    params.R0 = 2;
    params.Ub = 1e-3;
    params.a = 2.0;
    params.Varx = 0.1;
    params.asym = 1.0;
    params.dt = 1e-1;
    params.init_infected = 0.1;
    params.N_total = 1e8;

    OUTPUT_NAME = name;
    RADIUS = 2.0;

    string data_dir = CONST::OUTPUT_DIR_NAME + "/" + OUTPUT_NAME;

    InitDirectory(data_dir);
    EpidemicSimulator simulator(params);
    simulator.Run();

    // Вызов Python для визуализации
    string cmd = "python " + CONST::VISUALIZE_PYTHON +
        " --data-dir " + data_dir + "/" + CONST::DATA_DIR +
        " --animation " + data_dir + "/" + OUTPUT_NAME + ".gif" +
        " --snapshot " + data_dir + "/" + OUTPUT_NAME + ".png";
    system(cmd.c_str());
}

int main() {
    RunExperiment("test");
    return 0;
}