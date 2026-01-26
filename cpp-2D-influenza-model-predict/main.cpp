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
}

// Тип асимметрии иммунитета
enum AsymmetryType {
    ANI_ASY, // Анизотропная, асимметричная (основной режим статьи)
    ANI_SYM, // Анизотропная, симметричная
    ISO      // Изотропная (две антигенные координаты)
};

struct ModelParameters {
    // Размеры решётки
    int L; // Размер решётки L×L (в MATLAB коде L=50)
    // Параметры модели
    double R0;    // Базовый репродуктивный номер (2.6)
    double D;     // Мутационная вероятность (0.0001) - в MATLAB D
    double a;     // Полурасстояние иммунитета (7)
    double Varx;  // Разброс координат (0.01) - в MATLAB Varx
    // Параметры популяции
    double N; // Общая численность популяции (1e10)
    double init_infected; // Начальная доля инфицированных (1e-2)
    // Временные параметры
    int Tmax;  // Максимальное время (600)
    int stept; // Шаг по времени (1)
    int tshow; // Интервал отображения (100)
    int T0;    // Начало отображения (0)
    // Тип асимметрии
    AsymmetryType asymmetry; // Тип асимметрии иммунитета
    // Вычисляемые параметры
    int M; // Число шагов по времени = Tmax/stept
};

ModelParameters InitParameters() {
    ModelParameters params_;
    // Основные параметры из MATLAB кода
    params_.R0 = 2.6;
    params_.D = 0.0001;
    params_.a = 7;
    params_.Varx = 0.01;
    params_.N = 1e10;
    params_.init_infected = 1e-2;
    // Параметры решётки и времени
    params_.L = 50;
    params_.Tmax = 400;
    params_.stept = 1;
    params_.tshow = 100;
    params_.T0 = 0;
    // Вычисляемые параметры
    params_.M = params_.Tmax / params_.stept;
    // Тип асимметрии
    params_.asymmetry = ANI_ASY;
    return params_;
}

// Функция для вычисления суммы всех элементов матрицы
template <typename T>
T SumMatrix(const vector<vector<T>>& matrix) {
    return accumulate(matrix.begin(), matrix.end(), 0.0,
        [](T sum, const vector<T>& row) { return sum + accumulate(row.begin(), row.end(), 0.0); });
}

struct Point {
    double x;
    double y;
};

struct Stain {
    vector<vector<double>> dist_matrix;
};

class EpidemicSimulator {
public:
    EpidemicSimulator() {
        params_ = InitParameters();
        X = GenerateX(params_);
        Y = GenerateY(params_);
        K = GenerateImmunityMatrix(X, Y, params_);
        InitializeState(I, S, R, params_);
        norm.resize(params_.M, 0.0);
        finf.resize(params_.M, 0.0);
    }

    // Метод для запуска всей симуляции
    void Run() {
        cout << "========================================\n";
        cout << "TWO-DIMENSIONAL MODEL OF VIRUS ANTIGENIC EVOLUTION\n";
        cout << "========================================\n\n";

        // Шаг 1: Инициализация параметров
        cout << "1. Initializing model parameters...\n";
        PrintParameters();
        cout << "2. Initializing state matrices.\n";
        CheckInitialConditions();
        cout << "\n========================================\n";
        cout << "INITIALIZATION COMPLETE\n";
        cout << "Ready to start simulation of " << params_.M << " time steps\n";
        cout << "========================================\n";
        // Шаг 6: Запуск симуляции
        RunSimulation();
        // Шаг 7: Сохранение итоговых результатов
        cout << "3. Saving final results...\n";
        SaveFinalMatrices();
        SaveTimeSeriesData();

        PrintFinalStatistics();

        system(("python " + CONST::VISUALIZE_PYTHON).c_str());
    }
private:
    ModelParameters params_;
    vector<vector<vector<vector<double>>>> K;
    vector<vector<double>> X;
    vector<vector<double>> Y;
    vector<vector<double>> I;
    vector<vector<double>> S;
    vector<vector<double>> R;
    vector<double> norm;
    vector<double> finf;

    // Функция для сохранения матрицы в файл
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
    // Для X: зависит от столбца j
    vector<vector<double>> GenerateX(const ModelParameters& params_) {
        int L = params_.L;
        double Varx = params_.Varx;
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(-Varx, Varx);

        vector<vector<double>> X;
        for (int i = 0; i < L; ++i) {
            vector<double> row;
            for (int j = 0; j < L; ++j) {
                row.push_back((j + 1) + dis(gen));
            }
            X.push_back(row);
        }
        
        return X;
    }

    // Для Y: зависит от строки i
    vector<vector<double>> GenerateY(const ModelParameters& params_) {
        int L = params_.L;
        double Varx = params_.Varx;
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(-Varx, Varx);

        vector<vector<double>> Y;
        for (int i = 0; i < L; ++i) {
            vector<double> row;
            for (int j = 0; j < L; ++j) {
                row.push_back((i + 1) + dis(gen));
            }
            Y.push_back(row);
        }

        return Y;
    }
    // Генерация матрицы иммунитета K (4D матрица)
    vector<vector<vector<vector<double>>>> GenerateImmunityMatrix(const vector<vector<double>>& X, const vector<vector<double>>& Y, const ModelParameters& params_) {
        int L = params_.L;
        double a = params_.a;
        AsymmetryType asymmetry = params_.asymmetry;
        // Инициализация 4D матрицы: K[L][L][L][L]
        vector<vector<vector<vector<double>>>> K(L, vector<vector<vector<double>>>(L, vector<vector<double>>(L, vector<double>(L, 0.0))));
        // Оптимизированная версия с предвычислением X_ij
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                double X_ij = X[i][j];
                double Y_ij = Y[i][j]; // Для ISO режима
                vector<vector<double>> dist_matrix(L, vector<double>(L));
                // Сначала вычисляем матрицу расстояний для данной паре (i,j)
                for (int m = 0; m < L; ++m) {
                    for (int n = 0; n < L; ++n) {
                        switch (asymmetry) {
                        case ANI_SYM:
                        case ANI_ASY:
                            dist_matrix[m][n] = fabs(X[m][n] - X_ij) / a;
                            break;
                        case ISO:
                            double dx = X[m][n] - X_ij;
                            double dy = Y[m][n] - Y_ij;
                            dist_matrix[m][n] = sqrt(dx * dx + dy * dy) / a;
                            break;
                        }
                    }
                }

                // Затем заполняем K
                for (int m = 0; m < L; ++m) {
                    for (int n = 0; n < L; ++n) {
                        double dist = dist_matrix[m][n];
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

    // Инициализация начальных условий
    void InitializeState(vector<vector<double>>& I, vector<vector<double>>& S, vector<vector<double>>& R, const ModelParameters& params_) {
        int L = params_.L;
        // Инициализация нулями
        I = vector<vector<double>>(L, vector<double>(L, 0.0));
        S = vector<vector<double>>(L, vector<double>(L, 0.0));
        R = vector<vector<double>>(L, vector<double>(L, 0.0));
        // MATLAB: I(:,16)=1e-2/L; (индексирование с 1, у нас с 0)
        // Вертикальная линия инфицированных при x=15 (поскольку индексация с 0)
        int infected_line = 15; // соответствует 16 в MATLAB
        // MATLAB: S(:,1:15)=(1-1e-2)/15/L;
        // Полоса восприимчивых слева от инфицированной линии
        int susceptible_end = 15; // столбцы 0-14 (в MATLAB 1:15)
        // Распределение инфицированных
        for (int i = 0; i < L; ++i) {
            I[i][infected_line] = params_.init_infected / L;
        }

        // Распределение восприимчивых
        double susceptible_total = (1.0 - params_.init_infected);
        double susceptible_per_cell = susceptible_total / (susceptible_end * L);

        for (int j = 0; j < susceptible_end; ++j) {
            for (int i = 0; i < L; ++i) {
                S[i][j] = susceptible_per_cell;
            }
        }

        // В начальный момент выздоровевших нет
        // R уже инициализирован нулями
    }

    void RunSimulation() {
        int L = params_.L;
        int M = params_.M;
        double stept = params_.stept;
        double R0 = params_.R0;
        double D = params_.D;
        double N = params_.N;
        // Матрица для новых значений инфицированных
        vector<vector<double>> Inew(L, vector<double>(L, 0.0));
        auto sim_start = chrono::high_resolution_clock::now();
        // Основной цикл по времени
        for (int k = 0; k < M; ++k) {
            // Вычисление нормировки и доли инфицированных
            double total_I = SumMatrix(I);
            double total_S = SumMatrix(S);
            norm[k] = total_I + total_S;
            finf[k] = total_I / norm[k];
            // Создаём копию текущего I для мутационного члена
            vector<vector<double>> I_old = I;
            // Обнуляем матрицу Inew
            for (int i = 0; i < L; ++i) {
                for (int j = 0; j < L; ++j) {
                    Inew[i][j] = 0.0;
                }
            }

            // Основной цикл по всем клеткам решётки
            for (int i = 0; i < L; ++i) {
                for (int j = 0; j < L; ++j) {
                    // Вычисление Q: сумма по всем S[m][n] * K[i][j][m][n]
                    double Q = 0.0;
                    for (int m = 0; m < L; ++m) {
                        for (int n = 0; n < L; ++n) {
                            Q += S[m][n] * K[i][j][m][n];
                        }
                    }

                    // Вычисление P: сумма по всем I[m][n] * K[m][n][i][j]
                    double P = 0.0;
                    for (int m = 0; m < L; ++m) {
                        for (int n = 0; n < L; ++n) {
                            P += I_old[m][n] * K[m][n][i][j];
                        }
                    }

                    // Обновление инфицированных (уравнение для I)
                    // MATLAB: Inew(i,j) = I(i,j)*(1+stept*(R0*Q-1));
                    Inew[i][j] = I_old[i][j] * (1.0 + stept * (R0 * Q - 1.0));
                    // Обновление восприимчивых (уравнение для S)
                    // MATLAB: S(i,j) = S(i,j)*(1-stept*R0*P) + stept*I(i,j);
                    S[i][j] = S[i][j] * (1.0 - stept * R0 * P) + stept * I_old[i][j];
                    // Обновление выздоровевших (уравнение для R)
                    // MATLAB: R(i,j) = R(i,j) + stept*I(i,j);
                    R[i][j] = R[i][j] + stept * I_old[i][j];
                }
            }

            // Добавление мутационного члена (диффузия)
            // Только для внутренних клеток (не на границах)
            for (int i = 1; i < L - 1; ++i) {
                for (int j = 1; j < L - 1; ++j) {
                    // MATLAB: Inew(i,j) = Inew(i,j) + stept*D*(I(i,j+1)+I(i,j-1)+I(i-1,j)+I(i+1,j)-4*I(i,j));
                    double mutation_term = I_old[i][j + 1] + I_old[i][j - 1] + I_old[i - 1][j] + I_old[i + 1][j] - 4.0 * I_old[i][j];
                    Inew[i][j] += stept * D * mutation_term;
                }
            }

            // Отсечка по популяционному размеру
            // MATLAB: Inew = Inew.*(Inew/norm(k) > 1/N);
            for (int i = 0; i < L; ++i) {
                for (int j = 0; j < L; ++j) {
                    if (Inew[i][j] / norm[k] <= 1.0 / N) Inew[i][j] = 0.0;
                }
            }
            // Обновление матрицы инфицированных
            I = Inew;
            // Сохранение состояния на каждом tshow шаге (после T0)
            if (k * stept >= params_.T0 && (static_cast<int>(k * stept) % params_.tshow == 0)) {
                stringstream filename_I, filename_S, filename_R;
                filename_I << CONST::OUTPUT_DIR_NAME + "/state_I_step_" << k << ".csv";
                filename_S << CONST::OUTPUT_DIR_NAME + "/state_S_step_" << k << ".csv";
                filename_R << CONST::OUTPUT_DIR_NAME + "/state_R_step_" << k << ".csv";

                SaveMatrixToFile(I, filename_I.str());
                SaveMatrixToFile(S, filename_S.str());
                SaveMatrixToFile(R, filename_R.str());
            }

            // Прогресс каждые 10% симуляции
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
        int L = params_.L;
        // Итоговая статистика
        double final_total_I = SumMatrix(I);
        double final_total_S = SumMatrix(S);
        double final_total_R = SumMatrix(R);

        cout << "\n========================================\n";
        cout << "FINAL SIMULATION STATISTICS\n";
        cout << "========================================\n";
        cout << "Final infected fraction: " << final_total_I << endl;
        cout << "Final susceptible fraction: " << final_total_S << endl;
        cout << "Final recovered fraction: " << final_total_R << endl;
        cout << "Total (should be ~1.0): " << final_total_I + final_total_S + final_total_R << endl;

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
        // Сохранение временных рядов
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
        SaveMatrixToFile(S, CONST::OUTPUT_DIR_NAME + "/final_S.csv");
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
    }

    void CheckInitialConditions() {
        double total_I = 0.0;
        double total_S = 0.0;
        for (int i = 0; i < params_.L; ++i) {
            for (int j = 0; j < params_.L; ++j) {
                total_I += I[i][j];
                total_S += S[i][j];
            }
        }
        cout << "\tTotal infected: " << total_I << " (expected: " << params_.init_infected << ")\n";
        cout << "\tTotal susceptible: " << total_S << " (expected: " << 1.0 - params_.init_infected << ")\n";
        cout << "\tTotal population: " << total_I + total_S << " (should be 1.0)\n";
    }

};

int main() {
    srand(time(NULL));
    EpidemicSimulator simulator;
    simulator.Run();
    return 0;
}