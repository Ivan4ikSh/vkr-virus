#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <numeric>
#include <string>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

const string OUTPUT_DIR = "out";

void wave_simulation(double R0, double Ub, double a, int L, double N,
    double Tmax, double T0, double Txy, double stept,
    unsigned seed = random_device{}()) {
    fs::create_directories(OUTPUT_DIR);

    int M = static_cast<int>(round(Tmax / stept));

    // Два буфера для I и R (текущее и следующее состояние)
    vector<double> I_buf1(L, 0.0), I_buf2(L, 0.0);
    vector<double> R_buf1(L, 0.0), R_buf2(L, 0.0);

    // Указатели на текущие и следующие буферы (меняются местами на каждом шаге)
    vector<double>* I_cur = &I_buf1;
    vector<double>* I_next = &I_buf2;
    vector<double>* R_cur = &R_buf1;
    vector<double>* R_next = &R_buf2;

    // Начальные условия (MATLAB: I(20)=1e-3, остальное в R)
    (*I_cur)[19] = 1e-3;
    double initR = (1.0 - 1e-3) / 19.0;
    for (int j = 0; j < 19; ++j) {
        (*R_cur)[j] = initR;
    }

    // Ядро перекрёстного иммунитета K(d) = (d/a)/(1 + d/a)
    const vector<double> K([&] {
        vector<double> k(L);
        for (int d = 0; d < L; ++d) {
            double x = static_cast<double>(d) / a;
            k[d] = x / (1.0 + x);
        }
        return k;
        }());

    // Массивы для временных рядов
    vector<double> norm(M), finf(M), QI(M), PR(M);

    // Генератор случайных чисел
    mt19937 gen(seed);

    // Переменные для срезов
    vector<double> slice_times;
    vector<vector<double>> I_slices, R_slices;
    double stepshow = Tmax / 7.0;
    int show_step = static_cast<int>(round(stepshow / stept));

    double meanI0 = 0.0;
    int t0_step = static_cast<int>(round(T0 / stept)) - 1; // индекс шага, на котором t≈T0 (индекс с 0)

    // Лямбда для стохастической коррекции (broken stick)
    auto stochastic_correct = [&](double lambda) -> double {
        if (lambda >= 10.0) return lambda;
        double new_val = 0.0;
        int xm = static_cast<int>(round(6.0 * lambda));
        if (xm > 0) {
            double prob = lambda / xm;
            uniform_real_distribution<double> dist(0.0, 1.0);
            int count = 0;
            for (int n = 0; n < xm; ++n) {
                if (dist(gen) < prob) ++count;
            }
            new_val = static_cast<double>(count) / N;
        }
        return new_val;
        };

    // Основной цикл
    for (int i = 0; i < M; ++i) {
        double t = (i + 1) * stept;

        // Константные ссылки на текущие векторы (для удобства и скорости)
        const auto& I_old = *I_cur;
        const auto& R_old = *R_cur;

        double sumI_old = accumulate(I_old.begin(), I_old.end(), 0.0);
        double sumR_old = accumulate(R_old.begin(), R_old.end(), 0.0);
        double norm_i = sumI_old + sumR_old;

        norm[i] = norm_i;
        finf[i] = sumI_old / norm_i;

        vector<double> Q(L, 0.0), P(L, 0.0);
        auto& Inew = *I_next;  // ссылка на буфер следующего шага
        auto& Rnew = *R_next;

        // Расчёт Q и P для всех j
        for (int j = 0; j < L; ++j) {
            double Qj = 0.0, Pj = 0.0;

            // Левая часть для Q (R[0..j] * K[j], K[j-1], ...)
            for (int k = 0; k <= j; ++k) {
                Qj += R_old[k] * K[j - k];
            }
            // Правая часть для P (I[j..L-1] * K[0..L-1-j])
            for (int k = j; k < L; ++k) {
                Pj += I_old[k] * K[k - j];
            }

            if (t > Txy) {
                // Добавляем правую часть для Q (R[j..L-1] * K[0..L-1-j])
                for (int k = j; k < L; ++k) {
                    Qj += R_old[k] * K[k - j];
                }
                // Добавляем левую часть для P (I[0..j] * K[j], K[j-1], ...)
                for (int k = 0; k <= j; ++k) {
                    Pj += I_old[k] * K[j - k];
                }
            }

            Q[j] = Qj;
            P[j] = Pj;

            // Новое I без мутаций
            double newI = I_old[j] * (1.0 + stept * (R0 * Qj - 1.0));
            // Стохастическая коррекция
            newI = stochastic_correct(newI * N) / N;  // делим на N, т.к. stochastic_correct возвращает число особей, делённое на N

            Inew[j] = newI;

            // Обновление R
            Rnew[j] = R_old[j] * (1.0 - stept * R0 * Pj) + stept * I_old[j];

            // Мутационный член (кроме границ)
            if (j > 0 && j < L - 1) {
                double in = stept * Ub * (I_old[j + 1] + I_old[j - 1]);
                double out = stept * Ub * 2.0 * I_old[j];

                in = stochastic_correct(in * N) / N;  // коррекция только для in (как в оригинале)
                Inew[j] += in - out;
            }
        }

        // Граничные условия для мутаций (без стохастики, как в оригинале)
        Inew[0] = I_old[0] + stept * Ub * (I_old[1] - I_old[0]);
        Inew[L - 1] = I_old[L - 1] + stept * Ub * (I_old[L - 2] - I_old[L - 1]);

        // Нормировка I
        for (int j = 0; j < L; ++j) {
            Inew[j] /= norm_i;
        }

        // Вычисляем QI и PR
        double qival = 0.0, prval = 0.0;
        for (int j = 0; j < L; ++j) {
            qival += Q[j] * Inew[j];
            prval += P[j] * Rnew[j];
        }
        QI[i] = qival;
        PR[i] = prval;

        // Сохраняем среднюю позицию I в момент T0
        if (i == t0_step) {
            double sum_pos = 0.0, sum_I = 0.0;
            for (int j = 0; j < L; ++j) {
                sum_pos += Inew[j] * (j + 1);
                sum_I += Inew[j];
            }
            meanI0 = sum_pos / sum_I;
        }

        // Сохраняем срезы (каждые stepshow, начиная с T0)
        if (t > T0 && (i % show_step == 0)) {
            slice_times.push_back(t);
            I_slices.push_back(Inew);   // сохраняем новое состояние
            R_slices.push_back(Rnew);
        }

        // Меняем местами буферы для следующего шага
        swap(I_cur, I_next);
        swap(R_cur, R_next);

        // Прогресс
        if (i % 100 == 0) {
            cout << "t = " << t << ", norm = " << norm_i
                << ", finf = " << finf[i] << endl;
        }
    }

    // После цикла I_cur и R_cur указывают на финальное состояние
    const auto& I_final = *I_cur;
    const auto& R_final = *R_cur;

    double sumI_final = accumulate(I_final.begin(), I_final.end(), 0.0);
    double sumR_final = accumulate(R_final.begin(), R_final.end(), 0.0);

    double meanI = 0.0, meanR = 0.0;
    for (int j = 0; j < L; ++j) {
        meanI += I_final[j] * (j + 1);
        meanR += R_final[j] * (j + 1);
    }
    meanI /= sumI_final;
    meanR /= sumR_final;

    double sdI = 0.0, sdR = 0.0;
    for (int j = 0; j < L; ++j) {
        sdI += I_final[j] * ((j + 1) - meanI) * ((j + 1) - meanI);
        sdR += R_final[j] * ((j + 1) - meanR) * ((j + 1) - meanR);
    }
    sdI = sqrt(sdI / sumI_final);
    sdR = sqrt(sdR / sumR_final);

    double speed = (meanI - meanI0) / (Tmax - T0);
    double finf_final = finf.back();

    // Сохранение результатов (как в оригинале)
    ofstream f_norm(OUTPUT_DIR + "/norm.txt");
    for (double v : norm) f_norm << v << "\n";

    ofstream f_finf(OUTPUT_DIR + "/finf.txt");
    for (double v : finf) f_finf << v << "\n";

    ofstream f_QI(OUTPUT_DIR + "/QI.txt");
    for (double v : QI) f_QI << v << "\n";

    ofstream f_PR(OUTPUT_DIR + "/PR.txt");
    for (double v : PR) f_PR << v << "\n";

    ofstream f_I(OUTPUT_DIR + "/I_final.txt");
    for (double v : I_final) f_I << v << "\n";

    ofstream f_R(OUTPUT_DIR + "/R_final.txt");
    for (double v : R_final) f_R << v << "\n";

    ofstream f_slice_times(OUTPUT_DIR + "/slice_times.txt");
    for (double t : slice_times) f_slice_times << t << "\n";

    for (size_t idx = 0; idx < I_slices.size(); ++idx) {
        ofstream f_I_slice(OUTPUT_DIR + "/I_slice_" + to_string(idx) + ".txt");
        for (double v : I_slices[idx]) f_I_slice << v << "\n";

        ofstream f_R_slice(OUTPUT_DIR + "/R_slice_" + to_string(idx) + ".txt");
        for (double v : R_slices[idx]) f_R_slice << v << "\n";
    }

    cout << "\n=== Observables ===\n";
    cout << "speed = " << speed << "\n";
    cout << "sdI   = " << sdI << "\n";
    cout << "sdR   = " << sdR << "\n";
    cout << "finf  = " << finf_final << endl;
}

int main() {
    double R0 = 2.0;
    double Ub = 1e-3;
    double a = 7.0;
    int L = 120;
    double N = 1e8;
    double Tmax = 2000.0;
    double T0 = 400.0;
    double Txy = T0;
    double stept = 0.5;

    wave_simulation(R0, Ub, a, L, N, Tmax, T0, Txy, stept, 42);

    int ret = system("python plot_results.py");
    if (ret != 0) cerr << "Cannot load Python script" << endl;

    return 0;
}