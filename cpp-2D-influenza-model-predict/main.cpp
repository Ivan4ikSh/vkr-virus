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

using namespace std;

namespace CONST {
    const string OUTPUT_DIR_NAME = "output";
    const string VISUALIZE_PYTHON = "visualize.py";
    const double EPS = 1e-5;
}

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

    // Конструктор с инициализацией по умолчанию
    BoundingBox() : i_min(0), i_max(0), j_min(0), j_max(0), margin(1) {}

    // Метод для сброса границ
    void reset() {
        i_min = INT_MAX;
        i_max = INT_MIN;
        j_min = INT_MAX;
        j_max = INT_MIN;
    }

    // Метод для проверки, инициализирован ли bounding box
    bool is_valid() const {
        return i_min <= i_max && j_min <= j_max;
    }

    // Метод для расширения границ с учетом margin
    void expand(int L) {
        i_min = max(0, i_min - margin);
        i_max = min(L - 1, i_max + margin);
        j_min = max(0, j_min - margin);
        j_max = min(L - 1, j_max + margin);
    }

    // Метод для обновления границ с учетом новой точки
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
    int stept; // Шаг по времени (1)
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
    params_.Tmax = 600;
    params_.stept = 1;
    params_.tshow = 50;
    params_.T0 = 0;
    params_.M = params_.Tmax / params_.stept;
    params_.asymmetry = ANI_ASY;
    return params_;
}

// Функция для вычисления суммы всех элементов матрицы
template <typename T>
T SumMatrix(const vector<vector<T>>& matrix) {
    T sum = 0;
    for (const auto& row : matrix) {
        sum += accumulate(row.begin(), row.end(), 0.0);
    }
    return sum;
}

class EpidemicSimulator {
public:
    EpidemicSimulator() {
        params_ = InitParameters();
        X = GenerateX();
        Y = GenerateY();
        K = GenerateImmunityMatrix(X, Y);
        InitializeState(I, R);
        norm.resize(params_.M, 0.0);
        finf.resize(params_.M, 0.0);

        box.margin = 3;
        UpdateBoundingBox();
    }

    void Run() {
        cout << "========================================\n";
        cout << "TWO-DIMENSIONAL MODEL OF VIRUS ANTIGENIC EVOLUTION\n";
        cout << "========================================\n\n";

        cout << "1. Initializing model parameters...\n";
        PrintParameters();
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
    vector<vector<double>> I;
    vector<vector<double>> R;
    vector<double> norm;
    vector<double> finf;
    BoundingBox box;

    void UpdateBoundingBox() {
        bool found = false;
        int L = params_.L;

        // Проверяем обе матрицы
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                if (I[i][j] > CONST::EPS || R[i][j] > CONST::EPS) {
                    if (!found) {
                        box.reset();
                        found = true;
                    }
                    box.update(i, j);
                }
            }
        }

        // Если ничего не нашли, используем весь диапазон
        if (!found) {
            box.i_min = 0;
            box.i_max = L - 1;
            box.j_min = 0;
            box.j_max = L - 1;
        }

        // Расширяем с учетом margin
        box.expand(L);
    }

    void SaveMatrixToFile(const vector<vector<double>>& matrix, const string& filename) {
        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "Error: Could not open file " << filename << endl;
            return;
        }

        int rows = matrix.size();
        int cols = matrix[0].size();

        file << fixed << setprecision(6);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                file << matrix[i][j];
                if (j < cols - 1) file << ",";
            }
            file << endl;
        }
        file.close();
    }

    vector<vector<vector<vector<double>>>> GenerateImmunityMatrix(const vector<vector<double>>& X, const vector<vector<double>>& Y) {
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

    void InitializeState(vector<vector<double>>& I, vector<vector<double>>& R) {
        int L = params_.L;
        I = vector<vector<double>>(L, vector<double>(L, 0.0));
        R = vector<vector<double>>(L, vector<double>(L, 0.0));

        // Инициализация инфицированных: вертикальная линия при x=16 (индекс 15)
        int infected_line = 15; // соответствует столбцу 16 в MATLAB
        for (int i = 0; i < L; ++i) {
            I[i][infected_line] = params_.init_infected / L;
        }

        // Инициализация выздоровевших: полоса слева от инфицированных
        int susceptible_end = 15; // столбцы 0-14 (в MATLAB 1:15)
        double susceptible_total = (1.0 - params_.init_infected);
        double susceptible_per_cell = susceptible_total / (susceptible_end * L);

        for (int j = 0; j < susceptible_end; ++j) {
            for (int i = 0; i < L; ++i) {
                R[i][j] = susceptible_per_cell;
            }
        }
    }

    void RunSimulation() {
        int L = params_.L;
        int M = params_.M;
        double stept = params_.stept;
        double R0 = params_.R0;
        double D = params_.D;
        double N = params_.N;

        vector<vector<double>> Inew(L, vector<double>(L, 0.0));
        auto sim_start = chrono::high_resolution_clock::now();

        cout << "Initial bounding box: [" << box.i_min << ":" << box.i_max << ", " << box.j_min << ":" << box.j_max << "]" << endl;

        for (int k = 0; k < M; ++k) {
            // Вычисление нормировки и доли инфицированных
            double total_I = SumMatrix(I);
            double total_R = SumMatrix(R);
            norm[k] = total_I + total_R;
            finf[k] = total_I / norm[k];

            vector<vector<double>> I_old = I;

            // Обнуляем Inew только в bounding box
            for (int i = box.i_min; i <= box.i_max; ++i) {
                for (int j = box.j_min; j <= box.j_max; ++j) {
                    Inew[i][j] = 0.0;
                }
            }

            for (int i = box.i_min; i <= box.i_max; ++i) {
                for (int j = box.j_min; j <= box.j_max; ++j) {
                    // Q: сумма по всем R[m][n] * K[i][j][m][n]
                    double Q = 0.0;
                    // Ограничиваем и внутренний цикл bounding box'ом
                    for (int m = box.i_min; m <= box.i_max; ++m) {
                        for (int n = box.j_min; n <= box.j_max; ++n) {
                            Q += R[m][n] * K[i][j][m][n];
                        }
                    }

                    // P: сумма по всем I_old[m][n] * K[m][n][i][j]
                    double P = 0.0;
                    for (int m = box.i_min; m <= box.i_max; ++m) {
                        for (int n = box.j_min; n <= box.j_max; ++n) {
                            P += I_old[m][n] * K[m][n][i][j];
                        }
                    }

                    // Обновление инфицированных
                    Inew[i][j] = I_old[i][j] * (1.0 + stept * (R0 * Q - 1.0));
                    R[i][j] = R[i][j] * (1.0 - stept * R0 * P) + stept * I_old[i][j];
                }
            }

            // Добавление мутационного члена (диффузия) ТОЛЬКО внутри bounding box
            for (int i = max(1, box.i_min); i <= min(L - 2, box.i_max); ++i) {
                for (int j = max(1, box.j_min); j <= min(L - 2, box.j_max); ++j) {
                    double mutation_term = I_old[i][j + 1] + I_old[i][j - 1] + I_old[i - 1][j] + I_old[i + 1][j] - 4.0 * I_old[i][j];
                    Inew[i][j] += stept * D * mutation_term;
                }
            }

            // Отсечка по популяционному размеру ТОЛЬКО внутри bounding box
            for (int i = box.i_min; i <= box.i_max; ++i) {
                for (int j = box.j_min; j <= box.j_max; ++j) {
                    if (Inew[i][j] / norm[k] <= 1.0 / N) Inew[i][j] = 0.0;
                }
            }

            // Обновляем I значениями из Inew
            for (int i = box.i_min; i <= box.i_max; ++i) {
                for (int j = box.j_min; j <= box.j_max; ++j) {
                    I[i][j] = Inew[i][j];
                }
            }

            // Сохранение состояния
            if (k * stept >= params_.T0 && (static_cast<int>(k * stept) % params_.tshow == 0)) {
                stringstream filename_I, filename_R;
                filename_I << CONST::OUTPUT_DIR_NAME << "/state_I_step_" << k << ".csv";
                filename_R << CONST::OUTPUT_DIR_NAME << "/state_R_step_" << k << ".csv";

                SaveMatrixToFile(I, filename_I.str());
                SaveMatrixToFile(R, filename_R.str());
            }
            // Прогресс
            if (M >= 10 && k % (M / 10) == 0 && k > 0) {
                double progress = 100.0 * k / M;
                cout << "\tProgress: " << fixed << setprecision(1) << progress << "%" << endl;
                cout << "\tBounding box size: " << (box.i_max - box.i_min + 1) << "x" << (box.j_max - box.j_min + 1) << endl;
                cout << "\tBounding box: [" << box.i_min << ":" << box.i_max << ", " << box.j_min << ":" << box.j_max << "]" << endl;

            }

            UpdateBoundingBox();
        }

        auto sim_end = chrono::high_resolution_clock::now();
        auto sim_duration = chrono::duration_cast<chrono::seconds>(sim_end - sim_start);
        cout << "\nSimulation completed in " << sim_duration.count() << " seconds" << endl;
        cout << "Final bounding box: [" << box.i_min << ":" << box.i_max << ", " << box.j_min << ":" << box.j_max << "]" << endl;
    }

    void PrintFinalStatistics() {
        int L = params_.L;
        double final_total_I = SumMatrix(I);
        double final_total_R = SumMatrix(R);

        cout << "\n========================================\n";
        cout << "FINAL SIMULATION STATISTICS\n";
        cout << "========================================\n";
        cout << "Final infected fraction: " << final_total_I << endl;
        cout << "Final recovered fraction: " << final_total_R << endl;
        cout << "Total (should be ~1.0): " << final_total_I + final_total_R << endl;

        // Поиск максимума инфицированных
        double max_I = 0.0;
        int max_i = 0, max_j = 0;
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                if (I[i][j] > max_I) {
                    max_I = I[i][j];
                    max_i = i;
                    max_j = j;
                }
            }
        }
        cout << "Maximum infection at cell (" << max_i << "," << max_j << "): " << max_I << endl;
        cout << "Corresponding X coordinate: " << X[max_i][max_j] << endl;
        cout << "Corresponding Y coordinate: " << Y[max_i][max_j] << endl;
    }

    void SaveTimeSeriesData() {
        int M = params_.M;
        double stept = params_.stept;

        ofstream norm_file(CONST::OUTPUT_DIR_NAME + "/norm_time_series.csv");
        ofstream finf_file(CONST::OUTPUT_DIR_NAME + "/finf_time_series.csv");

        if (norm_file.is_open() && finf_file.is_open()) {
            norm_file << "Step,Time,Norm\n";
            finf_file << "Step,Time,Finf\n";

            for (int k = 0; k < M; ++k) {
                double time_val = k * stept;
                norm_file << k << "," << time_val << "," << norm[k] << "\n";
                finf_file << k << "," << time_val << "," << finf[k] << "\n";
            }

            norm_file.close();
            finf_file.close();
            cout << "\tSaved time series data\n";
        }
    }

    void SaveFinalMatrices() {
        SaveMatrixToFile(I, CONST::OUTPUT_DIR_NAME + "/final_I.csv");
        SaveMatrixToFile(R, CONST::OUTPUT_DIR_NAME + "/final_R.csv");
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
        cout << "\tTime step = " << params_.stept << endl;
        cout << "\tSave interval = " << params_.tshow << endl;
        cout << "\tBounding box margin: " << box.margin << endl;
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
};

int main() {
    srand(time(NULL));

    filesystem::create_directory(CONST::OUTPUT_DIR_NAME);
    EpidemicSimulator simulator;
    simulator.Run();

    system(("python " + CONST::VISUALIZE_PYTHON).c_str());
    return 0;
}