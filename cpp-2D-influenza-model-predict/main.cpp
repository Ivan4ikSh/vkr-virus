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
    double i_fitness;  // infected
    double s_fitness;  // susceptible/recovered
};

struct PointHash {
    size_t operator()(const Point& p) const {
        return hash<int>()(p.x) ^ (hash<int>()(p.y) << 1);
    }
};

struct ModelParameters {
    int L;                    // Размер решётки LxL
    int M;                    // Число шагов по времени
    int tshow;                // Интервал отображения
    int T0;                   // Начало отображения
    double R0;                // Базовый репродуктивный номер
    double Ub;                // Частота мутаций на геном за шаг
    double a;                 // Полурасстояние иммунитета
    double Varx;              // Разброс координат
    double asym;              // Коэффициент асимметрии
    double beta;              // Параметр формы распределения скачков
    double mutation_fraction; // Доля переносимых инфицированных при мутации
    double dt;                // Шаг по времени
    double init_infected;     // Начальная доля инфицированных
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
    exponential_distribution<double> exp_dist_;
    normal_distribution<double> normal_dist_;

public:
    EpidemicSimulator(ModelParameters params, int seed = time(nullptr)) : params_(params), generator_(seed),
        uniform_dist_(0.0, 1.0),
        exp_dist_(1.0),
        normal_dist_(0.0, 1.0) {

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
            // Сохранение состояния
            if (step >= params_.T0 && step % params_.tshow == 0) SaveCurrentState(step);
            // Расчет статистик
            CalculateStatistics(step);
            // Отображение прогресса
            if (params_.M >= 10 && step % (params_.M / 10) == 0) cout << "Step " << step << "/" << params_.M << " (" << (100 * step / params_.M) << "%)" << " - Infected cells: " << CountInfectedCells() << endl;
            // Один шаг симуляции
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

        strains_.clear();

        // Создаем начальное облако инфицированных
        double infected_per_cell = params_.init_infected / (CONST::PI * RADIUS * RADIUS);

        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                double dx = j - center_x;
                double dy = i - center_y;
                double dist2 = dx * dx + dy * dy;
                if (dist2 <= RADIUS * RADIUS) strains_[{ i, j }] = { infected_per_cell, 0.0 };
            }
        }

        // Распределяем восприимчивых по всем клеткам
        double susceptible_total = 1.0 - params_.init_infected;
        double susceptible_per_cell = susceptible_total / (L * L);

        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                if (strains_.find({ i, j }) == strains_.end()) strains_[{ i, j }] = { 0.0, susceptible_per_cell };
                else strains_[{ i, j }].s_fitness += susceptible_per_cell;
            }
        }

        norm_.resize(params_.M, 0.0);
        finf_.resize(params_.M, 0.0);
    }

    double KFunc(int i1, int j1, int i2, int j2) {
        double dx = X_[i2][j2] - X_[i1][j1];
        double dy = Y_[i2][j2] - Y_[i1][j1];
        double dist = sqrt(dx * dx + params_.asym * dy * dy);
        dist /= params_.a;
        return dist / (1.0 + dist);
    }

    pair<double, double> GenerateMutationJump() {
        double dx;
        double dy;

        if (params_.beta == 1) {
            // Экспоненциальное распределение
            double r = exp_dist_(generator_);
            double angle = uniform_dist_(generator_) * 2.0 * CONST::PI;
            dx = r * cos(angle);
            dy = r * sin(angle);
        }
        else if (params_.beta == 2) {
            // Гауссово распределение
            dx = normal_dist_(generator_);
            dy = normal_dist_(generator_);
        }
        else {
            // Общий случай (Вейбулла)
            double u = uniform_dist_(generator_);
            double r = pow(-log(u), 1.0 / params_.beta);
            double angle = uniform_dist_(generator_) * 2.0 * CONST::PI;
            dx = r * cos(angle);
            dy = r * sin(angle);
        }

        return { dx, dy };
    }

    void StepSimulation() {
        int L = params_.L;
        double dt = params_.dt;
        double R0 = params_.R0;

        unordered_map<Point, Strain, PointHash> new_strains;

        // Копируем текущее состояние
        for (const auto& [point, strain] : strains_) {
            new_strains[point] = strain;
        }

        // Обновляем каждую клетку
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                Point p{ i, j };

                double I = strains_[p].i_fitness;
                double S = strains_[p].s_fitness;

                if (I < CONST::EPS && S < CONST::EPS) continue;
                // Вычисляем Q (сумма по восприимчивым) и P (сумма по инфицированным)
                double Q = 0.0;
                double P = 0.0;
                for (const auto& [point, strain] : strains_) {
                    int m = point.x;
                    int n = point.y;
                    Q += strain.s_fitness * KFunc(i, j, m, n); // Q: сумма S * K для текущей клетки как источника инфекции
                    P += strain.i_fitness * KFunc(m, n, i, j); // P: сумма I * K для текущей клетки как цели заражения
                }
                double dI_dt = I * (R0 * Q - 1.0); // dI/dt = I*(R0*Q - 1)
                double new_I = I + dt * dI_dt;
                double dS_dt = I - R0 * S * P; // dS/dt = -R0*S*P + I
                double new_S = S + dt * dS_dt;
                // Мутации
                if (new_I > CONST::EPS && uniform_dist_(generator_) < params_.Ub * dt) {
                    auto [dx, dy] = GenerateMutationJump();
                    int new_i = i + (int)round(dx);
                    int new_j = j + (int)round(dy);

                    // Проверяем границы
                    if (new_i >= 0 && new_i < L && new_j >= 0 && new_j < L) {
                        Point new_p{ new_i, new_j };
                        double mutated = params_.mutation_fraction * new_I;
                        // Уменьшаем родительскую популяцию
                        new_I -= mutated;
                        // Добавляем мутанта
                        if (new_strains.find(new_p) == new_strains.end()) new_strains[new_p] = { mutated, 0.0 };
                        else new_strains[new_p].i_fitness += mutated;
                    }
                }
                // Обновляем значения
                if (new_I > CONST::EPS || new_S > CONST::EPS) new_strains[p] = { new_I, new_S };
            }
        }

        // Обновляем глобальное состояние
        strains_ = move(new_strains);
    }

    void CalculateStatistics(int step) {
        double total_I = 0.0, total_S = 0.0;

        for (const auto& [point, strain] : strains_) {
            total_I += strain.i_fitness;
            total_S += strain.s_fitness;
        }

        double total = total_I + total_S;
        norm_[step] = total;
        finf_[step] = (total > 0) ? total_I / total : 0.0;
    }

    int CountInfectedCells() const {
        int count = 0;
        for (const auto& [point, strain] : strains_) {
            if (strain.i_fitness > CONST::EPS) ++count;
        }
        return count;
    }

    void SaveCurrentState(int step) {
        int L = params_.L;

        // Создаем плотные матрицы
        vector<vector<double>> I_dense(L, vector<double>(L, 0.0));
        vector<vector<double>> S_dense(L, vector<double>(L, 0.0));

        for (const auto& [point, strain] : strains_) {
            int i = point.x;
            int j = point.y;
            if (i >= 0 && i < L && j >= 0 && j < L) {
                I_dense[i][j] = strain.i_fitness;
                S_dense[i][j] = strain.s_fitness;
            }
        }

        // Сохраняем
        string dir_path = CONST::OUTPUT_DIR_NAME + "/" + OUTPUT_NAME + "/" + CONST::DATA_DIR;
        filesystem::create_directories(dir_path);

        ofstream I_file(dir_path + "/state_I_step_" + to_string(step) + ".csv");
        ofstream S_file(dir_path + "/state_S_step_" + to_string(step) + ".csv");

        if (I_file.is_open() && S_file.is_open()) {
            for (int i = 0; i < L; ++i) {
                for (int j = 0; j < L; ++j) {
                    I_file << fixed << setprecision(6) << I_dense[i][j];
                    S_file << fixed << setprecision(6) << S_dense[i][j];
                    if (j < L - 1) {
                        I_file << ";";
                        S_file << ";";
                    }
                }
                I_file << "\n";
                S_file << "\n";
            }
        }
    }

    void SaveFinalResults() {
        string dir_path = CONST::OUTPUT_DIR_NAME + "/" + OUTPUT_NAME + "/" + CONST::DATA_DIR;

        // Сохраняем временные ряды
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

        // Сохраняем параметры
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
            param_file << "beta: " << params_.beta << "\n";
            param_file << "mutation_fraction: " << params_.mutation_fraction << "\n";
            param_file << "dt: " << params_.dt << "\n";
            param_file << "init_infected: " << params_.init_infected << "\n";
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
        cout << "  Mutation distribution (beta): " << params_.beta << "\n";
        cout << "  Initial infected: " << params_.init_infected << "\n";
    }

    void PrintFinalStatistics() {
        double total_I = 0.0;
        double total_S = 0.0;
        int infected_cells = 0;
        int susceptible_cells = 0;

        for (const auto& [point, strain] : strains_) {
            total_I += strain.i_fitness;
            total_S += strain.s_fitness;

            if (strain.i_fitness > CONST::EPS) ++infected_cells;
            if (strain.s_fitness > CONST::EPS) ++susceptible_cells;
        }

        cout << "\nFinal Statistics:\n";
        cout << "  Total infected: " << total_I << "\n";
        cout << "  Total susceptible: " << total_S << "\n";
        cout << "  Total population: " << total_I + total_S << "\n";
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
    params.M = 300;
    params.tshow = 10;
    params.T0 = 0;
    params.R0 = 3.0;
    params.Ub = 0.5;
    params.a = 15.0;
    params.Varx = 0.1;
    params.asym = 1.0;
    params.beta = 2;
    params.mutation_fraction = 0.5;
    params.dt = 0.1;
    params.init_infected = 0.01;

    OUTPUT_NAME = name;
    RADIUS = 1.0;

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
    RunExperiment("exp3");
    return 0;
}