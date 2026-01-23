#pragma once
#include "random.h"
#include "utils.h"
#include "framework.h"
#include "luksza.h"

class MonteCarlo {
public:
    MonteCarlo(const StartParams& start_params, LukszaParams luksza_params) : start_params_(start_params), luksza_params_(luksza_params), luksza_(start_params, luksza_params) {
        InitDirectory(start_params);
        InitStartParams();
    }

    void RunSimulation() {
        cout << "Monte-Carlo simulation starting...\n";
        int total = simulation_params_.T.size();
        int current = 0;

        PrintDebugHeader();
        for (const int& t : simulation_params_.T) {
            SimulationStep(t);
            ++current;
            PrintDebugStep(t);
            //PrintProgress(current, total);
        }
        cout << "\n";
        ComputeNumericalVelocity(start_params_, simulation_params_);
        FILEUTILS::SaveFigureData(start_params_, simulation_params_);
        cout << "Monte-Carlo simulation ended successfully!\n\n";
    }
private:
    StartParams start_params_;
    SimulationParams simulation_params_;
    LukszaParams luksza_params_;
    LukszaModel luksza_;

    void InitStartParams() {
        simulation_params_.mu = start_params_.muL / start_params_.L;
        simulation_params_.T.resize(start_params_.tf);
        iota(simulation_params_.T.begin(), simulation_params_.T.end(), 0);

        simulation_params_.s = InitAdaptiveLandscape(start_params_);
        simulation_params_.V_an = ComputeAnalyticalVelocity(start_params_);

        // Инициализация матрицы предков A: (1:N)'*ones(1,L)
        for (int i = 0; i < start_params_.N; ++i) {
            vector<int> loci;
            for (int j = 0; j < start_params_.L; ++j)
                loci.push_back(i);

            simulation_params_.A.push_back(loci);
        }

        simulation_params_.t_int = round(start_params_.tf / 10.0);
        simulation_params_.f_sample = 0.1;

        // Инициализация популяции K: (rand(N,L) < f0) или zeros(N,L)
        if (start_params_.f0 != 0.0) {
            for (int i = 0; i < start_params_.N; ++i) {
                vector<int> loci;
                for (int j = 0; j < start_params_.L; ++j)
                    loci.push_back((RNG.Random() < start_params_.f0) ? 1 : 0);

                simulation_params_.K.push_back(loci);
            }
        }
        else simulation_params_.K.resize(start_params_.N, vector<int>(start_params_.L, 0));

        int time_steps = start_params_.tf;

        simulation_params_.W.resize(start_params_.N, vector<double>(time_steps, 0.0));
        simulation_params_.P1.resize(start_params_.N, vector<int>(time_steps, 0));
        simulation_params_.PL.resize(start_params_.N, vector<int>(time_steps, 0));

        simulation_params_.f_site.resize(time_steps, vector<double>(start_params_.L, 0.0));
        simulation_params_.k_av.resize(time_steps, 0.0);
        simulation_params_.V_ark.resize(time_steps, 0.0);
        simulation_params_.f_survive.resize(time_steps, 0.0);
        simulation_params_.C.resize(time_steps, 0.0);
        simulation_params_.C_all.resize(time_steps, 0.0);
        simulation_params_.mean_W.resize(time_steps, 0.0);

        // Инициализация дополнительных векторов
        simulation_params_.dist_over_L.resize(start_params_.tf, 0.0);
        simulation_params_.f1_site.resize(start_params_.tf, 0.0);
    }

    void SimulationStep(const int& t) {
        const vector<vector<int>> K_prev = simulation_params_.K;
        // 1. Мутации
        Mutation(start_params_, simulation_params_);
        // 2. Начальные метки родителей: P = (1:N)'*ones(1,L)
        // Каждый локус имеет родителя - саму особь
        for (int i = 0; i < start_params_.N; ++i) {
            vector<int> loci;
            for (int j = 0; j < start_params_.L; ++j)
                loci.push_back(i);

            simulation_params_.P.push_back(loci);
        }
        // 3. Расчет приспособленностей: w = K * s'
        vector<double> w;
        if (CONST::USE_LUKSZA) w = luksza_.CalculateFitness(simulation_params_.K, simulation_params_.A, K_prev);
        else w = ComputeFitness(start_params_, simulation_params_);
        // 4. Расчет среднего числа потомков: nprogav = exp(w) / mean(exp(w))
        vector<double> n_prog_av = ComputeDescNumAv(start_params_, simulation_params_, w);
        // 5. Метод "сломанной палки": расчет кумулятивных сумм
        vector<int> n_prog = BrokenStick(n_prog_av);
        // 6. Обновление популяции (отбор)
        UpdatePopulation(n_prog);
        // 7. Рекомбинация
        if (start_params_.r > 0) Recombination();
        // 8. Запись наблюдаемых величин
        RecordObservables(t, w);
        if (simulation_params_.t_int * round(t / simulation_params_.t_int) == t) {
            auto hist = ComputeHistogram(w);
            simulation_params_.fitness_hist_xx.push_back(hist.first);
            simulation_params_.fitness_hist_nn.push_back(hist.second);
            simulation_params_.hist_times.push_back(t);

            // Гистограмма частот аллелей
            auto freq_hist = ComputeHistogram(simulation_params_.f_site[t]);
            simulation_params_.freq_hist_xx.push_back(freq_hist.first);
            simulation_params_.freq_hist_nn.push_back(freq_hist.second);
            simulation_params_.freq_hist_times.push_back(t);
        }
        // 9. Сохранение меток родителей для первого и последнего локусов
        for (int i = 0; i < start_params_.N; ++i) {
            simulation_params_.P1[i][t] = simulation_params_.P[i][0];
            simulation_params_.PL[i][t] = simulation_params_.P[i][start_params_.L - 1];
        }
        // 10. Сохранение приспособленностей
        for (int i = 0; i < start_params_.N; ++i) {
            simulation_params_.W[i][t] = w[i];
        }
    }

    vector<int> BrokenStick(const vector<double>& n_prog_av) {
        const int BUCKETS = min(1000, start_params_.N); // Можно подобрать оптимальное значение

        vector<double> b2(start_params_.N);
        partial_sum(n_prog_av.begin(), n_prog_av.end(), b2.begin());

        vector<double> b1(start_params_.N);
        b1[0] = 0.0;
        for (int i = 1; i < start_params_.N; ++i) b1[i] = b2[i - 1];

        // Создаем гистограмму (корзины)
        vector<int> histogram(BUCKETS, 0);
        for (int i = 0; i < start_params_.N; ++i) {
            double x = RNG.Random() * start_params_.N;
            int bucket = static_cast<int>(x * BUCKETS / start_params_.N);
            ++histogram[bucket];
        }

        // Распределяем по интервалам через гистограмму
        vector<int> n_prog(start_params_.N, 0);
        for (int i = 0; i < start_params_.N; ++i) {
            // Находим какие корзины попадают в интервал [b1[i], b2[i])
            int start_bucket = static_cast<int>(b1[i] * BUCKETS / start_params_.N);
            int end_bucket = static_cast<int>(b2[i] * BUCKETS / start_params_.N);

            // Распределяем содержимое корзин пропорционально
            for (int bucket = start_bucket; bucket <= end_bucket; ++bucket) {
                if (bucket < 0 || bucket >= BUCKETS) continue;

                // Вычисляем какую часть корзины занимает наш интервал
                double bucket_start = static_cast<double>(bucket) * start_params_.N / BUCKETS;
                double bucket_end = static_cast<double>(bucket + 1) * start_params_.N / BUCKETS;

                double overlap_start = max(b1[i], bucket_start);
                double overlap_end = min(b2[i], bucket_end);

                if (overlap_end > overlap_start) {
                    double overlap_ratio = (overlap_end - overlap_start) / (bucket_end - bucket_start);
                    n_prog[i] += static_cast<int>(histogram[bucket] * overlap_ratio + 0.5);
                }
            }
        }

        // Корректировка (добавляем недостающие/убираем лишние)
        int total_assigned = accumulate(n_prog.begin(), n_prog.end(), 0);
        int remaining = start_params_.N - total_assigned;

        if (remaining != 0) {
            // Простое распределение остатка
            for (int i = 0; i < abs(remaining); ++i) {
                int idx = RNG.RandomInt(0, start_params_.N - 1);
                if (remaining > 0) ++n_prog[idx];
                else if (n_prog[idx] > 0) --n_prog[idx];
            }
        }

        return n_prog;
    }

    void UpdatePopulation(const vector<int>& n_prog) {
        // Создаем временные матрицы
        simulation_params_.K_new.resize(start_params_.N, vector<int>(start_params_.L, 0));
        simulation_params_.A_new.resize(start_params_.N, vector<int>(start_params_.L, 0));
        simulation_params_.P_new.resize(start_params_.N, vector<int>(start_params_.L, 0));

        // Вычисляем кумулятивные суммы как в MATLAB
        vector<int> is(start_params_.N);
        is[0] = 0;
        for (int i = 1; i < start_params_.N; ++i) {
            is[i] = is[i - 1] + n_prog[i - 1];
        }

        // Заполняем новые матрицы
        for (int i = 0; i < start_params_.N; ++i) {
            if (n_prog[i] > 0) {
                // Для каждого потомка этого родителя
                for (int j = 0; j < n_prog[i]; ++j) {
                    int row_idx = is[i] + j;
                    if (row_idx < start_params_.N) {
                        simulation_params_.K_new[row_idx] = simulation_params_.K[i];
                        simulation_params_.A_new[row_idx] = simulation_params_.A[i];
                        simulation_params_.P_new[row_idx] = simulation_params_.P[i];
                    }
                }
            }
        }

        simulation_params_.K = move(simulation_params_.K_new);
        simulation_params_.A = move(simulation_params_.A_new);
        simulation_params_.P = move(simulation_params_.P_new);
    }

    // Функция рекомбинации
    void Recombination() {
        int npairs = round(start_params_.r * start_params_.N / 2.0);
        for (int pair_idx = 0; pair_idx < npairs; ++pair_idx) {
            // Выбираем двух случайных родителей
            int i1 = RNG.RandomInt(0, start_params_.N);
            int i2 = RNG.RandomInt(0, start_params_.N);

            // Генерируем вектор точек кроссинговера
            vector<int> xx(start_params_.L, 0);
            int count = 0;
            for (int locus = 0; locus < start_params_.L; ++locus) {
                if (RNG.Random() < static_cast<double>(start_params_.M) / start_params_.L) ++count;
                xx[locus] = count;
            }

            // Определяем, какие локусы берутся от первого родителя
            vector<bool> first(start_params_.L);
            for (int locus = 0; locus < start_params_.L; ++locus) {
                first[locus] = (round(xx[locus] / 2.0) == xx[locus] / 2.0);
            }

            // Создаем рекомбинантную последовательность
            vector<int> prog1_DNA(start_params_.L);
            vector<int> prog1_Ancestor(start_params_.L);
            vector<int> prog1_Parent(start_params_.L);

            for (int locus = 0; locus < start_params_.L; ++locus) {
                if (first[locus]) {
                    prog1_DNA[locus] = simulation_params_.K[i1][locus];
                    prog1_Ancestor[locus] = simulation_params_.A[i1][locus];
                    prog1_Parent[locus] = simulation_params_.P[i1][locus];
                }
                else {
                    prog1_DNA[locus] = simulation_params_.K[i2][locus];
                    prog1_Ancestor[locus] = simulation_params_.A[i2][locus];
                    prog1_Parent[locus] = simulation_params_.P[i2][locus];
                }
            }
            // Случайным образом заменяем одного из родителей
            if (RNG.Random() > 0.5) {
                simulation_params_.K[i1] = prog1_DNA;
                simulation_params_.A[i1] = prog1_Ancestor;
                simulation_params_.P[i1] = prog1_Parent;
            }
            else {
                simulation_params_.K[i2] = prog1_DNA;
                simulation_params_.A[i2] = prog1_Ancestor;
                simulation_params_.P[i2] = prog1_Parent;
            }
        }
    }

    // Запись наблюдаемых величин на каждом временном шаге
    void RecordObservables(const int& t, const vector<double>& w) {
        // 1. Частоты аллелей в каждом локусе
        for (int locus = 0; locus < start_params_.L; ++locus) {
            double allele_sum = accumulate(simulation_params_.K.begin(), simulation_params_.K.end(), 0.0,
                [locus](double sum, const vector<int>& genome) { return sum + genome[locus]; });
            simulation_params_.f_site[t][locus] = allele_sum / static_cast<double>(start_params_.N);
        }
        // 2. Среднее число благоприятных аллелей на геном
        double mean_freq = Mean(simulation_params_.f_site[t]);
        simulation_params_.k_av[t] = mean_freq * start_params_.L;
        // 3. Дисперсия приспособленности (нормированная на s0^2)
        vector<int> k_values(start_params_.N, 0);
        for (int i = 0; i < start_params_.N; ++i) {
            for (int j = 0; j < start_params_.L; ++j) {
                k_values[i] += simulation_params_.K[i][j];
            }
        }
        // 3. Вычисляем дисперсию k
        double mean_w = Mean(w);
        double var_w = 0.0;
        for (int i = 0; i < start_params_.N; ++i) {
            double diff = w[i] - mean_w;
            var_w += diff * diff;
        }
        if (start_params_.N > 1) var_w /= start_params_.N;
        simulation_params_.V_ark[t] = var_w / (start_params_.s0 * start_params_.s0);

        // 4. Доля локусов, в которых есть хотя бы один благоприятный аллель
        int genomes_with_at_least_one = 0;
        for (int j = 0; j < start_params_.L; ++j) {
            bool found_one = false;
            for (int i = 0; i < start_params_.N; ++i) {
                if (simulation_params_.K[i][j] == 1) {
                    found_one = true;
                    break;
                }
            }
            if (found_one) ++genomes_with_at_least_one;
        }
        simulation_params_.f_survive[t] = static_cast<double>(genomes_with_at_least_one) / start_params_.L;

        // 5. Доля пар с общим предком (выборка)
        int sample_size = round(start_params_.N * simulation_params_.f_sample);
        double total_matches = 0.0;
        for (int pair_idx = 0; pair_idx < sample_size; ++pair_idx) {
            // Используем RandomInt для генерации индексов от 0 до N-1
            int idx1 = RNG.RandomInt(0, start_params_.N - 1);
            int idx2 = RNG.RandomInt(0, start_params_.N - 1);
            int matches_in_pair = 0;
            for (int locus = 0; locus < start_params_.L; ++locus) {
                if (simulation_params_.A[idx1][locus] == simulation_params_.A[idx2][locus]) ++matches_in_pair;
            }
            total_matches += static_cast<double>(matches_in_pair) / start_params_.L;
        }
        simulation_params_.C[t] = total_matches / sample_size;

        // 6. Доля локусов, где все особи имеют общего предка (через std)
        double common_all = 0.0;
        for (int locus = 0; locus < start_params_.L; ++locus) {
            // Собираем метки предков для данного локуса
            vector<double> ancestors(start_params_.N);
            for (int i = 0; i < start_params_.N; ++i) {
                ancestors[i] = static_cast<double>(simulation_params_.A[i][locus]);
            }

            // Вычисляем среднее
            double mean = 0.0;
            for (double anc : ancestors) mean += anc;
            mean /= start_params_.N;

            // Вычисляем дисперсию (несмещенная оценка)
            double variance = 0.0;
            for (double anc : ancestors) {
                double diff = anc - mean;
                variance += diff * diff;
            }
            variance /= (start_params_.N - 1);

            // Если стандартное отклонение близко к 0 (с учетом погрешности)
            if (variance < 1e-12) ++common_all;
        }
        simulation_params_.C_all[t] = common_all / start_params_.L;
        // 7. Средняя приспособленность
        simulation_params_.mean_W[t] = mean_w;
        // Вычисляем dist (генетическое разнообразие)
        double dist_sum = 0.0;
        for (int locus = 0; locus < start_params_.L; ++locus) {
            double f = simulation_params_.f_site[t][locus];
            dist_sum += 2.0 * f * (1.0 - f);
        }
        simulation_params_.dist_over_L[t] = dist_sum / start_params_.L; // dist/L
        // Вычисляем теоретическую частоту f1site
        simulation_params_.f1_site[t] = start_params_.f0 / (start_params_.f0 + (1.0 - start_params_.f0) * exp(-start_params_.s0 * t));
    }

    void PrintProgress(int current, int total) {
        if (total == 0) return;

        int barWidth = 50;  // Ширина прогресс-бара
        float progress = static_cast<float>(current) / total;
        int filled = static_cast<int>(barWidth * progress);

        cout << "[";
        for (int i = 0; i < barWidth; ++i) {
            if (i < filled) cout << "=";
            else if (i == filled) cout << ">";
            else cout << " ";
        }
        cout << "] " << setw(3) << static_cast<int>(progress * 100.0) << "%\r";
        cout.flush();  // Принудительно выводим на экран
    }

    void PrintDebugHeader() {
        cout << "================================================================================\n";
        cout << "DEBUG OUTPUT FOR MONTE CARLO SIMULATION\n";
        cout << "Parameters: N=" << start_params_.N
            << ", L=" << start_params_.L
            << ", s0=" << start_params_.s0
            << ", r=" << start_params_.r
            << ", f0=" << start_params_.f0
            << ", muL=" << start_params_.muL
            << ", tf=" << start_params_.tf
            << ", M=" << start_params_.M
            << ", distribution=" << start_params_.distribution << "\n";
        cout << "================================================================================\n";
        cout << setw(4) << "t"
            << setw(10) << "k_av"
            << setw(10) << "k_av/L"
            << setw(10) << "Var_k"
            << setw(10) << "mean_W"
            << setw(12) << "f_survive"
            << setw(8) << "C"
            << setw(8) << "C_all"
            << setw(10) << "sum(K)"
            << setw(10) << "min_f"
            << setw(10) << "max_f"
            << setw(10) << "mean_f"
            << "\n";
        cout << "--------------------------------------------------------------------------------\n";
    }

    void PrintDebugStep(int t) {
        // Вычисляем дополнительные статистики для отладки
        double min_freq = 1.0;
        double max_freq = 0.0;
        double mean_freq = 0.0;

        for (int locus = 0; locus < start_params_.L; ++locus) {
            double freq = simulation_params_.f_site[t][locus];
            mean_freq += freq;
            if (freq < min_freq) min_freq = freq;
            if (freq > max_freq) max_freq = freq;
        }
        mean_freq /= start_params_.L;

        double sum_K = 0;
        for (int i = 0; i < start_params_.N; ++i) {
            for (int j = 0; j < start_params_.L; ++j) {
                sum_K += simulation_params_.K[i][j];
            }
        }

        cout << setw(4) << t
            << setw(10) << fixed << setprecision(2) << simulation_params_.k_av[t]
            << setw(10) << fixed << setprecision(4) << simulation_params_.k_av[t] / start_params_.L
            << setw(10) << fixed << setprecision(4) << simulation_params_.V_ark[t]
            << setw(10) << fixed << setprecision(4) << simulation_params_.mean_W[t]
            << setw(12) << fixed << setprecision(4) << simulation_params_.f_survive[t]
            << setw(8) << fixed << setprecision(4) << simulation_params_.C[t]
            << setw(8) << fixed << setprecision(4) << simulation_params_.C_all[t]
            << setw(10) << fixed << setprecision(0) << sum_K
            << setw(10) << fixed << setprecision(4) << min_freq
            << setw(10) << fixed << setprecision(4) << max_freq
            << setw(10) << fixed << setprecision(4) << mean_freq
            << "\n";
    }
};
