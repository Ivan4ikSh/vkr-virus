#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <unordered_map>
#include <sstream>
#include <string>
#include <vector>
#include <set>

using namespace std;

namespace CONST {
    const string VISUALIZE_PYTHON = "visualize.py";
    const double EPS = 1e-10;
}

string OUTPUT_NAME = "output";

// Тип асимметрии иммунитета
enum AsymmetryType {
    ANI_ASY, // Анизотропная, асимметричная (основной режим статьи)
    ANI_SYM, // Анизотропная, симметричная
    ISO      // Изотропная (две антигенные координаты)
};

struct BoundingBox {
    int i_min, i_max; // Границы по Y
    int j_min, j_max; // Границы по X
    int margin;       // Запас вокруг волны

    BoundingBox() : i_min(0), i_max(0), j_min(0), j_max(0), margin(5) {}

    void reset() {
        i_min = INT_MAX;
        i_max = INT_MIN;
        j_min = INT_MAX;
        j_max = INT_MIN;
    }

    bool is_valid() const {
        return i_min <= i_max && j_min <= j_max;
    }

    void expand(int L) {
        i_min = max(0, i_min - margin);
        i_max = min(L - 1, i_max + margin);
        j_min = max(0, j_min - margin);
        j_max = min(L - 1, j_max + margin);
    }

    void update(int i, int j) {
        i_min = min(i_min, i);
        i_max = max(i_max, i);
        j_min = min(j_min, j);
        j_max = max(j_max, j);
    }
};

struct ModelParameters {
    int L; // Размер решётки L×L (50)
    double R0;    // Базовый репродуктивный номер (2.6)
    double D;     // Мутационная вероятность (0.0001)
    double a;     // Полурасстояние иммунитета (7)
    double Varx;  // Разброс координат (0.01)
    double N; // Общая численность популяции (1e10)
    double init_infected; // Начальная доля инфицированных (1e-2)
    int Tmax;  // Максимальное время (600)
    int tshow; // Интервал отображения (100)
    int T0;    // Начало отображения (0)
    AsymmetryType asymmetry; // Тип асимметрии иммунитета
    int M; // Число шагов по времени = Tmax/stept
};

ModelParameters InitParameters() {
    ModelParameters params_;
    params_.R0 = 2.6;
    params_.D = 0.0001;
    params_.a = 7;
    params_.Varx = 0.01;
    params_.N = 1e10;
    params_.init_infected = 1e-2;
    params_.L = 50;
    params_.Tmax = 1200;
    params_.tshow = 50;
    params_.T0 = 0;
    params_.M = params_.Tmax;
    params_.asymmetry = ANI_ASY;
    return params_;
}

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ (h2 << 1);
    }
};

using SparseMatrix = unordered_map<pair<double, double>, double, pair_hash>;

class EpidemicSimulator {
public:
    EpidemicSimulator() {
        params_= InitParameters();
        X = GenerateX();
        Y = GenerateY();
        K = GenerateImmunityMatrix();
        InitializeState(I_sparse, R_sparse);
        norm.resize(params_.M, 0.0);
        finf.resize(params_.M, 0.0);
    }

    void Run() {
        cout << "========================================\n";
        cout << "TWO-DIMENSIONAL MODEL OF VIRUS ANTIGENIC EVOLUTION\n";
        cout << "========================================\n\n";

        cout << "1. Initializing model parameters...\n";
        PrintParameters();
        cout << "2. Initializing sparse state matrices.\n";
        CheckInitialConditions();
        cout << "\n========================================\n";
        cout << "INITIALIZATION COMPLETE\n";
        cout << "Ready to start simulation of " << params_.M << " time steps\n";
        cout << "========================================\n";

        RunSimulation();

        cout << "3. Saving final results...\n";
        SaveFinalMatrices();
        SaveTimeSeriesData();

        PrintFinalStatistics();
    }

private:
    ModelParameters params_;
    vector<vector<vector<vector<double>>>> K;
    vector<vector<double>> X;
    vector<vector<double>> Y;
    SparseMatrix I_sparse; // Разреженная матрица инфицированных
    SparseMatrix R_sparse; // Разреженная матрица выздоровевших
    vector<double> norm;
    vector<double> finf;

    // Сохранение sparse матрицы в файл
    void SaveSparseMatrixToFile(const SparseMatrix& matrix, const string& filename) {
        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "Error: Could not open file " << filename << endl;
            return;
        }

        // Преобразуем sparse в плотную матрицу для сохранения
        vector<vector<double>> dense(params_.L, vector<double>(params_.L, 0.0));
        for (const auto& entry : matrix) {
            int i = entry.first.first;
            int j = entry.first.second;
            dense[i][j] = entry.second;
        }

        file << fixed << setprecision(6);
        for (int i = 0; i < params_.L; ++i) {
            for (int j = 0; j < params_.L; ++j) {
                file << dense[i][j];
                if (j < params_.L - 1) file << ",";
            }
            file << endl;
        }
        file.close();
    }

    vector<vector<double>> GenerateX() {
        int L = params_.L;
        double Varx = params_.Varx;
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(-1.0, 1.0);

        vector<vector<double>> X(L, vector<double>(L));
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                X[i][j] = j + Varx * dis(gen);
            }
        }
        return X;
    }

    vector<vector<double>> GenerateY() {
        int L = params_.L;
        double Varx = params_.Varx;
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(-1.0, 1.0);

        vector<vector<double>> Y(L, vector<double>(L));
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                Y[i][j] = i + Varx * dis(gen);
            }
        }
        return Y;
    }

    vector<vector<vector<vector<double>>>> GenerateImmunityMatrix() {
        int L = params_.L;
        double a = params_.a;
        AsymmetryType asymmetry = params_.asymmetry;
        vector<vector<vector<vector<double>>>> K(L, vector<vector<vector<double>>>(L, vector<vector<double>>(L, vector<double>(L, 0.0))));

        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                double X_ij = X[i][j];
                double Y_ij = Y[i][j];

                for (int m = 0; m < L; ++m) {
                    for (int n = 0; n < L; ++n) {
                        double dist = 0.0;
                        switch (asymmetry) {
                        case ANI_SYM:
                        case ANI_ASY:
                            dist = fabs(X[m][n] - X_ij) / a;
                            break;
                        case ISO:
                            double dx = X[m][n] - X_ij;
                            double dy = Y[m][n] - Y_ij;
                            dist = sqrt(dx * dx + dy * dy) / a;
                            break;
                        }

                        if (asymmetry == ANI_ASY) {
                            if (X[m][n] < X_ij) K[i][j][m][n] = dist / (1 + dist);
                            else K[i][j][m][n] = 0.0;
                        }
                        else K[i][j][m][n] = dist / (1 + dist);
                    }
                }
            }
        }
        return K;
    }

    void InitializeState(SparseMatrix& I, SparseMatrix& R) {
        int L = params_.L;
        I.clear();
        R.clear();

        // Инициализация инфицированных: вертикальная линия при x=16 (индекс 15)
        int infected_line = 15;
        double infected_per_cell = params_.init_infected / L;
        for (int i = 0; i < L; ++i) {
            I[{i, infected_line}] = infected_per_cell;
        }

        // Инициализация выздоровевших: полоса слева от инфицированных
        int susceptible_end = 15;
        double susceptible_total = (1.0 - params_.init_infected);
        double susceptible_per_cell = susceptible_total / (susceptible_end * L);

        for (int j = 0; j < susceptible_end; ++j) {
            for (int i = 0; i < L; ++i) {
                R[{i, j}] = susceptible_per_cell;
            }
        }
    }

    // Получить множество всех клеток, которые могут быть затронуты на текущем шаге
    set<pair<double, double>> GetAffectedCells() {
        set<pair<double, double>> affected;

        // Все текущие ненулевые клетки I и R
        for (const auto& entry : I_sparse) {
            affected.insert(entry.first);
        }
        for (const auto& entry : R_sparse) {
            affected.insert(entry.first);
        }

        // Добавляем соседей для учёта мутаций
        set<pair<double, double>> neighbors;
        for (const auto& cell : affected) {
            int i = cell.first;
            int j = cell.second;
            if (i > 0) neighbors.insert({ i - 1, j });
            if (i < params_.L - 1) neighbors.insert({ i + 1, j });
            if (j > 0) neighbors.insert({ i, j - 1 });
            if (j < params_.L - 1) neighbors.insert({ i, j + 1 });
        }

        affected.insert(neighbors.begin(), neighbors.end());
        return affected;
    }

    void RunSimulation() {
        int L = params_.L;
        int M = params_.M;
        double R0 = params_.R0;
        double D = params_.D;
        double N = params_.N;

        auto sim_start = chrono::high_resolution_clock::now();

        for (int k = 0; k < M; ++k) {
            // Вычисление нормировки и доли инфицированных
            double total_I = 0.0;
            for (const auto& entry : I_sparse) total_I += entry.second;

            double total_R = 0.0;
            for (const auto& entry : R_sparse) total_R += entry.second;

            norm[k] = total_I + total_R;
            finf[k] = total_I / norm[k];

            SparseMatrix I_old = I_sparse;
            SparseMatrix R_old = R_sparse;

            SparseMatrix Inew;
            SparseMatrix Rnew = R_old;

            // Получаем все клетки, которые нужно обновить
            set<pair<double, double>> affected_cells = GetAffectedCells();

            // Предварительно вычисляем суммы для Q и P для каждой затронутой клетки
            SparseMatrix Q_values, P_values;

            for (const auto& cell : affected_cells) {
                int i = cell.first;
                int j = cell.second;

                // Вычисляем Q: сумма по всем R[m][n] * K[i][j][m][n]
                double Q = 0.0;
                for (const auto& R_entry : R_old) {
                    int m = R_entry.first.first;
                    int n = R_entry.first.second;
                    Q += R_entry.second * K[i][j][m][n];
                }
                Q_values[{i, j}] = Q;

                // Вычисляем P: сумма по всем I_old[m][n] * K[m][n][i][j]
                double P = 0.0;
                for (const auto& I_entry : I_old) {
                    int m = I_entry.first.first;
                    int n = I_entry.first.second;
                    P += I_entry.second * K[m][n][i][j];
                }
                P_values[{i, j}] = P;
            }

            // Обновляем значения для затронутых клеток
            for (const auto& cell : affected_cells) {
                int i = cell.first;
                int j = cell.second;

                double I_old_val = (I_old.find(cell) != I_old.end()) ? I_old[cell] : 0.0;
                double R_old_val = (R_old.find(cell) != R_old.end()) ? R_old[cell] : 0.0;

                double Q = Q_values[{i, j}];
                double P = P_values[{i, j}];

                // Обновление инфицированных
                double I_new_val = I_old_val * R0 * Q;

                // Мутации
                if (i > 0 && i < L - 1 && j > 0 && j < L - 1) {
                    double mutation_term = 0.0;
                    mutation_term += (I_old.find({ i, j + 1 }) != I_old.end()) ? I_old[{i, j + 1}] : 0.0;
                    mutation_term += (I_old.find({ i, j - 1 }) != I_old.end()) ? I_old[{i, j - 1}] : 0.0;
                    mutation_term += (I_old.find({ i - 1, j }) != I_old.end()) ? I_old[{i - 1, j}] : 0.0;
                    mutation_term += (I_old.find({ i + 1, j }) != I_old.end()) ? I_old[{i + 1, j}] : 0.0;
                    mutation_term -= 4.0 * I_old_val;
                    I_new_val += D * mutation_term;
                }

                // Обновление выздоровевших
                double R_new_val = R_old_val * (1.0 - R0 * P) + I_old_val;

                // Применяем порог отсечения
                if (I_new_val / norm[k] > 1.0 / N) Inew[cell] = I_new_val;

                if (R_new_val > CONST::EPS) Rnew[cell] = R_new_val;
                else Rnew.erase(cell);
            }

            // Обновляем глобальные матрицы
            I_sparse = Inew;
            R_sparse = Rnew;

            // Сохранение состояния
            if (k >= params_.T0 && (static_cast<int>(k) % params_.tshow == 0)) {
                stringstream filename_I, filename_R;
                filename_I << OUTPUT_NAME << "/state_I_step_" << k << ".csv";
                filename_R << OUTPUT_NAME << "/state_R_step_" << k << ".csv";

                SaveSparseMatrixToFile(I_sparse, filename_I.str());
                SaveSparseMatrixToFile(R_sparse, filename_R.str());
            }

            // Прогресс
            if (M >= 10 && k % (M / 10) == 0 && k > 0) {
                double progress = 100.0 * k / M;
                cout << "\tProgress: " << fixed << setprecision(1) << progress << "%" << endl;
            }
        }

        auto sim_end = chrono::high_resolution_clock::now();
        auto sim_duration = chrono::duration_cast<chrono::seconds>(sim_end - sim_start);
        cout << "\nSimulation completed in " << sim_duration.count() << " seconds" << endl;
    }

    void PrintFinalStatistics() {
        double final_total_I = 0.0;
        for (const auto& entry : I_sparse) final_total_I += entry.second;

        double final_total_R = 0.0;
        for (const auto& entry : R_sparse) final_total_R += entry.second;

        cout << "\n========================================\n";
        cout << "FINAL SIMULATION STATISTICS\n";
        cout << "========================================\n";
        cout << "Final infected fraction: " << final_total_I << endl;
        cout << "Final recovered fraction: " << final_total_R << endl;
        cout << "Total (should be ~1.0): " << final_total_I + final_total_R << endl;
        cout << "Non-zero infected cells: " << I_sparse.size() << endl;
        cout << "Non-zero recovered cells: " << R_sparse.size() << endl;
    }

    void SaveTimeSeriesData() {
        int M = params_.M;

        ofstream norm_file(OUTPUT_NAME + "/norm_time_series.csv");
        ofstream finf_file(OUTPUT_NAME + "/finf_time_series.csv");

        if (norm_file.is_open() && finf_file.is_open()) {
            norm_file << "Step,Time,Norm\n";
            finf_file << "Step,Time,Finf\n";

            for (int k = 0; k < M; ++k) {
                double time_val = k;
                norm_file << k << "," << time_val << "," << norm[k] << "\n";
                finf_file << k << "," << time_val << "," << finf[k] << "\n";
            }

            norm_file.close();
            finf_file.close();
            cout << "\tSaved time series data\n";
        }
    }

    void SaveFinalMatrices() {
        SaveSparseMatrixToFile(I_sparse, OUTPUT_NAME + "/final_I.csv");
        SaveSparseMatrixToFile(R_sparse, OUTPUT_NAME + "/final_R.csv");
    }

    void PrintParameters() {
        cout << fixed << setprecision(6);
        cout << "\tL = " << params_.L << "x" << params_.L << endl;
        cout << "\tR0 = " << params_.R0 << endl;
        cout << "\tD = " << params_.D << endl;
        cout << "\ta = " << params_.a << endl;
        cout << "\tVarx = " << params_.Varx << endl;
        cout << "\tN = " << scientific << params_.N << fixed << endl;
        cout << "\tTmax = " << params_.Tmax << endl;
        cout << "\tSteps = " << params_.M << endl;
        cout << "\tSave interval = " << params_.tshow << endl;
    }

    void CheckInitialConditions() {
        double total_I = 0.0;
        double total_R = 0.0;
        for (const auto& entry : I_sparse) total_I += entry.second;
        for (const auto& entry : R_sparse) total_R += entry.second;

        cout << "\tTotal infected: " << total_I << " (expected: " << params_.init_infected << ")\n";
        cout << "\tTotal recovered: " << total_R << " (expected: " << 1.0 - params_.init_infected << ")\n";
        cout << "\tTotal population: " << total_I + total_R << " (should be 1.0)\n";
        cout << "\tInitial non-zero infected cells: " << I_sparse.size() << endl;
        cout << "\tInitial non-zero recovered cells: " << R_sparse.size() << endl;
    }
};

void TEST(const string& experiment_name) {
    OUTPUT_NAME = experiment_name;
    filesystem::create_directory(OUTPUT_NAME);
 
    EpidemicSimulator simulator;
    simulator.Run();

    system(("python " + CONST::VISUALIZE_PYTHON + " --data-dir " + OUTPUT_NAME + " --animation " + OUTPUT_NAME + ".gif --snapshot " + OUTPUT_NAME + ".png").c_str());
}

int main() {
    srand(time(NULL));
    TEST("exp1");
    TEST("exp2");
    TEST("exp3");
    TEST("exp4");
    TEST("exp5");
    return 0;
}