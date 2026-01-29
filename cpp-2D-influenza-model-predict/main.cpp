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
    const string OUTPUT_DIR_NAME = "output";
    const string DATA_DIR = "data";
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
    int L;                   // Размер решётки L×L (50)
    double R0;               // Базовый репродуктивный номер (2.6)
    double D;                // Мутационная вероятность (0.0001)
    double a;                // Полурасстояние иммунитета (7)
    double Varx;             // Разброс координат (0.01)
    double N;                // Общая численность популяции (1e10)
    double init_infected;    // Начальная доля инфицированных (1e-2)
    int tshow;               // Интервал отображения (100)
    int T0;                  // Начало отображения (0)
    AsymmetryType asymmetry; // Тип асимметрии иммунитета
    int M;                   // Число шагов по времени = Tmax/stept
};

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
    EpidemicSimulator(ModelParameters params) : params_(params){
        GenerateX();
        GenerateY();
        K_ = GenerateImmunityMatrix();
        InitializeState();
        norm_.resize(params_.M, 0.0);
        finf_.resize(params_.M, 0.0);
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

        string param_file = CONST::OUTPUT_DIR_NAME + "/" + OUTPUT_NAME + "/" + CONST::DATA_DIR + "/model_parameters.txt";
        SaveParametersToFile(param_file);

        PrintFinalStatistics();
    }

private:
    ModelParameters params_;
    vector<vector<vector<vector<double>>>> K_;
    vector<vector<double>> X_;
    vector<vector<double>> Y_;
    SparseMatrix I_;
    SparseMatrix R_;
    vector<double> norm_;
    vector<double> finf_;

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
            for (const auto& entry : I_) total_I += entry.second;

            double total_R = 0.0;
            for (const auto& entry : R_) total_R += entry.second;

            norm_[k] = total_I + total_R;
            finf_[k] = total_I / norm_[k];

            SparseMatrix I_old = I_;
            SparseMatrix R_old = R_;

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
                    Q += R_entry.second * K_[i][j][m][n];
                }
                Q_values[{i, j}] = Q;

                // Вычисляем P: сумма по всем I_old[m][n] * K[m][n][i][j]
                double P = 0.0;
                for (const auto& I_entry : I_old) {
                    int m = I_entry.first.first;
                    int n = I_entry.first.second;
                    P += I_entry.second * K_[m][n][i][j];
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
                if (I_new_val / norm_[k] > 1.0 / N) Inew[cell] = I_new_val;

                if (R_new_val > CONST::EPS) Rnew[cell] = R_new_val;
                else Rnew.erase(cell);
            }

            // Обновляем глобальные матрицы
            I_ = Inew;
            R_ = Rnew;

            // Сохранение состояния
            if (k >= params_.T0 && (static_cast<int>(k) % params_.tshow == 0)) {
                stringstream filename_I, filename_R;
                filename_I << CONST::OUTPUT_DIR_NAME + "/" + OUTPUT_NAME + "/" + CONST::DATA_DIR + "/state_I_step_" << k << ".csv";
                filename_R << CONST::OUTPUT_DIR_NAME + "/" + OUTPUT_NAME + "/" + CONST::DATA_DIR + "/state_R_step_" << k << ".csv";

                SaveSparseMatrixToFile(I_, filename_I.str());
                SaveSparseMatrixToFile(R_, filename_R.str());
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

    void InitializeState() {
        int L = params_.L;
        I_.clear();
        R_.clear();

        // Инициализация инфицированных: вертикальная линия при x=16 (индекс 15)
        int infected_line = 15;
        double infected_per_cell = params_.init_infected / L;
        for (int i = 0; i < L; ++i) {
            I_[{i, infected_line}] = infected_per_cell;
        }

        // Инициализация выздоровевших: полоса слева от инфицированных
        int susceptible_end = infected_line;
        double susceptible_total = (1.0 - params_.init_infected);
        double susceptible_per_cell = susceptible_total / (susceptible_end * L);

        for (int j = 0; j < susceptible_end; ++j) {
            for (int i = 0; i < L; ++i) {
                R_[{i, j}] = susceptible_per_cell;
            }
        }
    }


    void GenerateX() {
        int L = params_.L;
        double Varx = params_.Varx;
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(-1.0, 1.0);

        for (int i = 0; i < L; ++i) {
            vector<double> row;
            for (int j = 0; j < L; ++j) {
                row.push_back(j + Varx * dis(gen));
            }
            X_.push_back(row);
        }
    }

    void GenerateY() {
        int L = params_.L;
        double Varx = params_.Varx;
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(-1.0, 1.0);

        for (int i = 0; i < L; ++i) {
            vector<double> row;
            for (int j = 0; j < L; ++j) {
                row.push_back(i + Varx * dis(gen));
            }
            Y_.push_back(row);
        }
    }

    vector<vector<vector<vector<double>>>> GenerateImmunityMatrix() {
        int L = params_.L;
        double a = params_.a;
        AsymmetryType asymmetry = params_.asymmetry;
        vector<vector<vector<vector<double>>>> K(L, vector<vector<vector<double>>>(L, vector<vector<double>>(L, vector<double>(L, 0.0))));

        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                double X_ij = X_[i][j];
                double Y_ij = Y_[i][j];

                for (int m = 0; m < L; ++m) {
                    for (int n = 0; n < L; ++n) {
                        double dist = 0.0;
                        switch (asymmetry) {
                        case ANI_SYM:
                        case ANI_ASY:
                            dist = fabs(X_[m][n] - X_ij) / a;
                            break;
                        case ISO:
                            double dx = X_[m][n] - X_ij;
                            double dy = Y_[m][n] - Y_ij;
                            dist = sqrt(dx * dx + dy * dy) / a;
                            break;
                        }

                        if (asymmetry == ANI_ASY) {
                            if (X_[m][n] < X_ij) K[i][j][m][n] = dist / (1 + dist);
                            else K[i][j][m][n] = 0.0;
                        }
                        else K[i][j][m][n] = dist / (1 + dist);
                    }
                }
            }
        }
        return K;
    }

    // Получить множество всех клеток, которые могут быть затронуты на текущем шаге
    set<pair<double, double>> GetAffectedCells() {
        set<pair<double, double>> affected;

        // Все текущие ненулевые клетки I и R
        for (const auto& entry : I_) {
            affected.insert(entry.first);
        }
        for (const auto& entry : R_) {
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

    void PrintFinalStatistics() {
        double final_total_I = 0.0;
        for (const auto& entry : I_) final_total_I += entry.second;

        double final_total_R = 0.0;
        for (const auto& entry : R_) final_total_R += entry.second;

        cout << "\n========================================\n";
        cout << "FINAL SIMULATION STATISTICS\n";
        cout << "========================================\n";
        cout << "Final infected fraction: " << final_total_I << endl;
        cout << "Final recovered fraction: " << final_total_R << endl;
        cout << "Total (should be ~1.0): " << final_total_I + final_total_R << endl;
        cout << "Non-zero infected cells: " << I_.size() << endl;
        cout << "Non-zero recovered cells: " << R_.size() << endl;
    }

    void SaveTimeSeriesData() {
        int M = params_.M;

        ofstream norm__file(OUTPUT_NAME + "/norm__time_series.csv");
        ofstream finf__file(OUTPUT_NAME + "/finf__time_series.csv");

        if (norm__file.is_open() && finf__file.is_open()) {
            norm__file << "Step,Time,norm_\n";
            finf__file << "Step,Time,finf_\n";

            for (int k = 0; k < M; ++k) {
                double time_val = k;
                norm__file << k << "," << time_val << "," << norm_[k] << "\n";
                finf__file << k << "," << time_val << "," << finf_[k] << "\n";
            }

            norm__file.close();
            finf__file.close();
            cout << "\tSaved time series data\n";
        }
    }

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

    void SaveFinalMatrices() {
        SaveSparseMatrixToFile(I_, CONST::OUTPUT_DIR_NAME + "/" + OUTPUT_NAME + "/" + CONST::DATA_DIR + "/final_I.csv");
        SaveSparseMatrixToFile(R_, CONST::OUTPUT_DIR_NAME + "/" + OUTPUT_NAME + "/" + CONST::DATA_DIR + "/final_R.csv");
    }
    
    void SaveParametersToFile(const string& filename) {
        ofstream file(filename);
        if (file.is_open()) {
            file << fixed << setprecision(6);
            file << "MODEL_PARAMETERS:\n";
            file << "===============\n";
            file << "Grid size L: " << params_.L << "x" << params_.L << "\n";
            file << "Reproductive number R0: " << params_.R0 << "\n";
            file << "Mutation probability D: " << params_.D << "\n";
            file << "Immunity distance a: " << params_.a << "\n";
            file << "Coordinate spread Varx: " << params_.Varx << "\n";
            file << "Total population N: " << scientific << params_.N << fixed << "\n";
            file << "Initial infected fraction: " << params_.init_infected << "\n";
            file << "Show interval tshow: " << params_.tshow << "\n";
            file << "Start time T0: " << params_.T0 << "\n";
            file << "Total steps M: " << params_.M << "\n";

            // Тип асимметрии
            string asymmetry_str;
            switch (params_.asymmetry) {
            case ANI_ASY: asymmetry_str = "Anisotropic Asymmetric (Main)"; break;
            case ANI_SYM: asymmetry_str = "Anisotropic Symmetric"; break;
            case ISO: asymmetry_str = "Isotropic (2 antigenic coordinates)"; break;
            default: asymmetry_str = "Unknown";
            }
            file << "Asymmetry type: " << asymmetry_str << "\n";

            file << "\nSIMULATION_INFO:\n";
            file << "===============\n";
            file << "Output directory: " << CONST::OUTPUT_DIR_NAME + "/" + OUTPUT_NAME << "\n";
            file << "Data subdirectory: " << CONST::DATA_DIR << "\n";
            file << "Total steps saved: " << min(params_.M, (int)(params_.M / params_.tshow) + 1) << "\n";

            file.close();
            cout << "\tModel parameters saved to " << filename << endl;
        }
        else {
            cerr << "Error: Could not open parameters file " << filename << endl;
        }
    }

    void PrintParameters() {
        cout << fixed << setprecision(6);
        cout << "\tL = " << params_.L << "x" << params_.L << endl;
        cout << "\tR0 = " << params_.R0 << endl;
        cout << "\tD = " << params_.D << endl;
        cout << "\ta = " << params_.a << endl;
        cout << "\tVarx = " << params_.Varx << endl;
        cout << "\tN = " << scientific << params_.N << fixed << endl;
        cout << "\tSteps = " << params_.M << endl;
        cout << "\tSave interval = " << params_.tshow << endl;
    }

    void CheckInitialConditions() {
        double total_I = 0.0;
        double total_R = 0.0;
        for (const auto& entry : I_) total_I += entry.second;
        for (const auto& entry : R_) total_R += entry.second;

        cout << "\tTotal infected: " << total_I << " (expected: " << params_.init_infected << ")\n";
        cout << "\tTotal recovered: " << total_R << " (expected: " << 1.0 - params_.init_infected << ")\n";
        cout << "\tTotal population: " << total_I + total_R << " (should be 1.0)\n";
        cout << "\tInitial non-zero infected cells: " << I_.size() << endl;
        cout << "\tInitial non-zero recovered cells: " << R_.size() << endl;
    }
};

void ClearDirectory(const filesystem::path& directory_path) {
    for (const auto& entry : filesystem::recursive_directory_iterator(directory_path)) {
        if (filesystem::is_regular_file(entry.path())) filesystem::remove(entry.path());
    }
}


void InitDirectory(const string& data_dir) {
    filesystem::create_directory(CONST::OUTPUT_DIR_NAME);
    filesystem::create_directory(data_dir);
    filesystem::create_directory(data_dir + "/" + CONST::DATA_DIR);
    ClearDirectory(data_dir);
}

void TEST(const string& experiment_name) {
    ModelParameters params;
    params.R0 = 3.1;
    params.D = 0.0001;
    params.a = 7;
    params.Varx = 0.01;
    params.N = 1e10;
    params.init_infected = 1e-2;
    params.L = 50;
    params.tshow = 20;
    params.T0 = 0;
    params.M = 1200;
    params.asymmetry = ISO;

    OUTPUT_NAME = experiment_name;
    string data_dir = CONST::OUTPUT_DIR_NAME + "/" + OUTPUT_NAME;
    InitDirectory(data_dir);
    
    EpidemicSimulator simulator(params);
    simulator.Run();
        
    system(("python " + CONST::VISUALIZE_PYTHON + " --data-dir " + (data_dir + "/" + CONST::DATA_DIR) + " --animation " + data_dir + "/" + OUTPUT_NAME + ".gif --snapshot " + data_dir + "/" + OUTPUT_NAME + ".png").c_str());
}

int main() {
    srand(time(NULL));
    TEST("exp2");
    //TEST("exp2");
    //TEST("exp3");
    return 0;
}