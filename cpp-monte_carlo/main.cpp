#include "domain.h"
#include "random.h"
#include "utils.h"

uint32_t SEED = RANDOM::Hash();
RandomGenerator RNG(SEED);

// Функция для вычисления гистограммы
pair<vector<double>, vector<double>> ComputeHistogram(const vector<double>& data, int bins = 50) {
    if (data.empty()) return { {}, {} };

    double min_val = *min_element(data.begin(), data.end());
    double max_val = *max_element(data.begin(), data.end());

    if (min_val == max_val) {
        min_val -= 0.5;
        max_val += 0.5;
    }

    double bin_width = (max_val - min_val) / bins;
    vector<double> xx(bins);
    vector<double> nn(bins, 0.0);

    // Центры бинов
    for (int i = 0; i < bins; ++i) {
        xx[i] = min_val + (i + 0.5) * bin_width;
    }

    // Подсчет
    for (double val : data) {
        int bin_idx = min(max(static_cast<int>((val - min_val) / bin_width), 0), bins - 1);
        nn[bin_idx] += 1.0;
    }

    return { xx, nn };
}

// Функция для вычисления численной скорости адаптации
void ComputeNumericalVelocity(const StartParams& params, SimulationParams& sim_params) {
    if (sim_params.T.size() < 3) {
        sim_params.V_num = 0.0;
        return;
    }

    vector<double> V_num_vec(sim_params.T.size(), 0.0);
    // Центральная разность для внутренних точек
    for (size_t i = 1; i < sim_params.T.size() - 1; ++i) {
        double dt = sim_params.T[i + 1] - sim_params.T[i - 1];
        V_num_vec[i] = (sim_params.k_av[i + 1] - sim_params.k_av[i - 1]) / dt;
    }
    // Односторонние разности для границ
    V_num_vec[0] = (sim_params.k_av[1] - sim_params.k_av[0]) / (sim_params.T[1] - sim_params.T[0]);
    V_num_vec.back() = (sim_params.k_av.back() - sim_params.k_av[sim_params.k_av.size() - 2]) / (sim_params.T.back() - sim_params.T[sim_params.T.size() - 2]);
    // Используем ТОЛЬКО последнюю треть времени (установившийся режим)
    int t_start = round(2.0 * sim_params.T.size() / 3.0);
    double sum = 0.0;
    int count = 0;

    for (size_t i = t_start; i < V_num_vec.size(); ++i) {
        if (!isnan(V_num_vec[i]) && !isinf(V_num_vec[i])) {
            sum += V_num_vec[i];
            count++;
        }
    }

    sim_params.V_num = (count > 0) ? sum / count : 0.0;
}

// Функция для вычисления аналитической скорости адаптации
double ComputeAnalyticalVelocity(const StartParams& params) {
    double Ub = params.muL * (1.0 - params.f0);
    // Вычисление первого аргумента: Nsqrt(sUb)
    double N_sqrt_sUb = static_cast<double>(params.N) * sqrt(params.s0 * Ub);
    double log_N_sqrt_sUb = log(N_sqrt_sUb);
    // Второй логарифмический член в знаменателе
    double s_Ub_ratio = params.s0 / Ub;
    double arg_log2 = s_Ub_ratio * log_N_sqrt_sUb;
    double log_arg_log2 = log(arg_log2);
    // Вычисление скорости по формуле
    double numerator = 2.0 * params.s0 * log_N_sqrt_sUb;
    double denominator = log_arg_log2 * log_arg_log2;
    double V = numerator / denominator;

    // В конце ComputeAnalyticalVelocity
    if (log_arg_log2 < 3.0 || log_N_sqrt_sUb < 3.0) {
        cout << "WARNING! Velocity may be comuted with error\n";
        cout << "  Required: log(arg_log2) >= 3\n";
        cout << "  Actual:   log(arg_log2) = " << log_arg_log2 << "\n";
        cout << "  Required: log(N_sqrt_sUb) >= 3\n";
        cout << "  Actual:   log(N_sqrt_sUb) = " << log_N_sqrt_sUb << "\n";
        cout << "  s/Ub = " << s_Ub_ratio << " (should be >> 1)\n";
    }

    return V;
}

void InitDirectory(const StartParams& params) {
    if (!FILEUTILS::FileExists(params.output_dir)) FILEUTILS::CreateDirectory(params.output_dir);
    if (!FILEUTILS::FileExists(params.plot_dir)) FILEUTILS::CreateDirectory(params.plot_dir);
}

vector<double> InitAdaptiveLandscape(const StartParams& params) {
    vector<double> s(params.L);
    switch (params.distribution)
    {
    case DistributionType::CONSTANT:
    {
        // Постоянный коэффициент отбора для всех локусов
        for (int i = 0; i < params.L; ++i)
            s[i] = params.s0;
    }
    break;
    case DistributionType::EXPONENTIAL:
    {
        // Экспоненциальное распределение : s = -s0 * ln(rand)
        for (int i = 0; i < params.L; ++i)
            s[i] = -params.s0 * log(RNG.Random());
    }
    break;
    case DistributionType::HALF_GAUSSIAN:
    {
        // Полугауссово распределение: PDF = (2/sqrt(2π)/s0) * exp(-s²/(2s0²))
            // Для получения среднего s0 нужно умножить на sqrt(pi/2)
        double factor = params.s0 * sqrt(CONST::PI / 2.0);
        for (int i = 0; i < params.L; ++i)
            s[i] = factor * fabs(RNG.Normal(0.0, 1.0));
    }
    break;
    default:
        break;
    }

    return s;
}

SimulationParams InitStartParams(const StartParams& params) {
    SimulationParams sim_params;
    sim_params.mu = params.muL / params.L;
    sim_params.T.resize(params.tf + 1);
    iota(sim_params.T.begin(), sim_params.T.end(), 0);

    sim_params.s = InitAdaptiveLandscape(params);
    sim_params.V_an = ComputeAnalyticalVelocity(params);


    // Инициализация матрицы предков A: (1:N)'*ones(1,L)
    sim_params.A.resize(params.N, vector<int>(params.L));
    for (int i = 0; i < params.N; ++i)
        for (int j = 0; j < params.L; ++j)
            sim_params.A[i][j] = i + 1;  // +1 для соответствия MATLAB (индексация с 1)

    sim_params.t_int = round(params.tf / 10.0);
    sim_params.colors = "rgbmkrgbmkrgbmkrgbmk";
    sim_params.f_sample = 0.1;

    // Инициализация популяции K: (rand(N,L) < f0) или zeros(N,L)
    sim_params.K.resize(params.N, vector<int>(params.L, 0));

    if (params.f0 != 0.0) {
        for (int i = 0; i < params.N; ++i) {
            for (int j = 0; j < params.L; ++j) {
                sim_params.K[i][j] = (RNG.Random() < params.f0) ? 1 : 0;
            }
        }
    }

    int time_steps = params.tf + 1;

    sim_params.W.resize(params.N, vector<double>(time_steps, 0.0));
    sim_params.P1.resize(params.N, vector<int>(time_steps, 0));
    sim_params.PL.resize(params.N, vector<int>(time_steps, 0));

    sim_params.f_site.resize(time_steps, vector<double>(params.L, 0.0));
    sim_params.k_av.resize(time_steps, 0.0);
    sim_params.V_ark.resize(time_steps, 0.0);
    sim_params.f_survive.resize(time_steps, 0.0);
    sim_params.C.resize(time_steps, 0.0);
    sim_params.C_all.resize(time_steps, 0.0);
    sim_params.mean_W.resize(time_steps, 0.0);

    // Инициализация дополнительных векторов
    sim_params.dist_over_L.resize(params.tf + 1, 0.0);
    sim_params.f1_site.resize(params.tf + 1, 0.0);
    sim_params.C.resize(params.tf + 1, 0.0);
    sim_params.C_all.resize(params.tf + 1, 0.0);
    sim_params.f_survive.resize(params.tf + 1, 0.0);

    return sim_params;
}

void Mutation(const StartParams& params, SimulationParams& simulation_params) {
    if (simulation_params.mu > 0) {
        for (int i = 0; i < params.N; ++i) {
            for (int j = 0; j < params.L; ++j) {
                double temp = RNG.Random();
                if (RNG.Random() < simulation_params.mu)
                    simulation_params.K[i][j] = simulation_params.K[i][j] ^ 1;  // XOR для инвертирования бита
            }
        }
    }
}

vector<int> BrokenStick(const StartParams& params, const vector<double>& n_prog_av) {
    vector<double> b1(params.N), b2(params.N);
    vector<int> n_prog(params.N, 0);

    // Рассчитываем кумулятивные суммы
    b2[0] = n_prog_av[0];
    for (int i = 1; i < params.N; ++i) {
        b2[i] = b2[i - 1] + n_prog_av[i];
    }

    b1[0] = 0.0;
    for (int i = 1; i < params.N; ++i) {
        b1[i] = b2[i - 1];
    }

    // Генерация случайных чисел X ~ U(0, N)
    vector<double> X(params.N);
    for (int i = 0; i < params.N; ++i) {
        X[i] = RNG.Random() * params.N;
    }

    // Расчет числа потомков для каждой особи
    // Более эффективная версия с поиском бинарным поиском
    for (int i = 0; i < params.N; ++i) {
        // Для каждой особи i ищем, в какой интервал попадают точки X
        for (int j = 0; j < params.N; ++j) {
            if (X[j] > b1[i] && X[j] < b2[i]) {
                n_prog[i]++;
            }
        }
    }

    return n_prog;
}

void UpdatePopulation(const StartParams& params, SimulationParams& sim_params, const vector<int>& n_prog) {
    // Создаем временные матрицы для нового поколения
    sim_params.K_new.resize(params.N, vector<int>(params.L, 0));
    sim_params.A_new.resize(params.N, vector<int>(params.L, 0));
    sim_params.P_new.resize(params.N, vector<int>(params.L, 0));

    // Создаем вектор начальных индексов для каждого родителя
    vector<int> is(params.N, 0);
    is[0] = 0;
    for (int i = 1; i < params.N; ++i) {
        is[i] = is[i - 1] + n_prog[i - 1];
    }
    // Заполняем новые матрицы потомками
    int current_pos = 0;
    for (int i = 0; i < params.N; ++i) {
        if (n_prog[i] > 0) {
            for (int j = 0; j < n_prog[i]; ++j) {
                int row_idx = is[i] + j;
                if (row_idx >= params.N) break; // Защита от выхода за границы
                // Копируем геном
                sim_params.K_new[row_idx] = sim_params.K[i];
                // Копируем метки предков
                sim_params.A_new[row_idx] = sim_params.A[i];
                // Копируем метки родителей
                sim_params.P_new[row_idx] = sim_params.P[i];
            }
        }
    }

    // Заменяем старые матрицы новыми
    sim_params.K = move(sim_params.K_new);
    sim_params.A = move(sim_params.A_new);
    sim_params.P = move(sim_params.P_new);
}

// Функция рекомбинации
void Recombination(const StartParams& params, SimulationParams& sim_params) {
    int npairs = round(params.r * params.N / 2.0);
    for (int pair_idx = 0; pair_idx < npairs; ++pair_idx) {
        // Выбираем двух случайных родителей
        int i1 = RNG.RandomInt(0, params.N - 1);
        int i2 = RNG.RandomInt(0, params.N - 1);

        // Генерируем вектор точек кроссинговера
        vector<int> xx(params.L, 0);
        int count = 0;
        for (int locus = 0; locus < params.L; ++locus) {
            if (RNG.Random() < static_cast<double>(params.M) / params.L) ++count;
            xx[locus] = count;
        }

        // Определяем, какие локусы берутся от первого родителя
        vector<bool> first(params.L);
        for (int locus = 0; locus < params.L; ++locus) {
            first[locus] = (round(xx[locus] / 2.0) == xx[locus] / 2.0);
        }

        // Создаем рекомбинантную последовательность
        vector<int> prog1_DNA(params.L);
        vector<int> prog1_Ancestor(params.L);
        vector<int> prog1_Parent(params.L);

        for (int locus = 0; locus < params.L; ++locus) {
            if (first[locus]) {
                prog1_DNA[locus] = sim_params.K[i1][locus];
                prog1_Ancestor[locus] = sim_params.A[i1][locus];
                prog1_Parent[locus] = sim_params.P[i1][locus];
            }
            else {
                prog1_DNA[locus] = sim_params.K[i2][locus];
                prog1_Ancestor[locus] = sim_params.A[i2][locus];
                prog1_Parent[locus] = sim_params.P[i2][locus];
            }
        }
        // Случайным образом заменяем одного из родителей
        if (RNG.Random() > 0.5) {
            sim_params.K[i1] = prog1_DNA;
            sim_params.A[i1] = prog1_Ancestor;
            sim_params.P[i1] = prog1_Parent;
        }
        else {
            sim_params.K[i2] = prog1_DNA;
            sim_params.A[i2] = prog1_Ancestor;
            sim_params.P[i2] = prog1_Parent;
        }
    }
}

// Запись наблюдаемых величин на каждом временном шаге
void RecordObservables(const StartParams& params, SimulationParams& sim_params, const int& t, const vector<double>& w) {
    // 1. Частоты аллелей в каждом локусе
    vector<double> locus_frequencies(params.L, 0.0);
    for (int locus = 0; locus < params.L; ++locus) {
        double sum = 0.0;
        for (int i = 0; i < params.N; ++i) {
            sum += sim_params.K[i][locus];
        }
        double freq = sum / params.N;
        sim_params.f_site[t][locus] = freq;
        locus_frequencies[locus] = freq;
    }
    // 2. Среднее число благоприятных аллелей на геном
    double mean_freq = 0.0;
    for (int locus = 0; locus < params.L; ++locus) {
        mean_freq += locus_frequencies[locus];
    }
    mean_freq /= params.L;
    sim_params.k_av[t] = mean_freq * params.L;
    // 3. Дисперсия приспособленности (нормированная на s0^2)
    double mean_w = 0.0;
    for (int i = 0; i < params.N; ++i) {
        mean_w += w[i];
    }
    mean_w /= params.N;

    double var_w = 0.0;
    for (int i = 0; i < params.N; ++i) {
        var_w += (w[i] - mean_w) * (w[i] - mean_w);
    }
    var_w /= params.N;
    sim_params.V_ark[t] = (sqrt(var_w) / params.s0) * (sqrt(var_w) / params.s0);
    // 4. Доля выживших генотипов (не все нули)
    int count_nonzero = 0;
    for (int i = 0; i < params.N; ++i) {
        bool all_zeros = true;
        for (int locus = 0; locus < params.L; ++locus) {
            if (sim_params.K[i][locus] != 0) {
                all_zeros = false;
                break;
            }
        }
        if (!all_zeros) ++count_nonzero;
    }
    sim_params.f_survive[t] = static_cast<double>(count_nonzero) / params.N;
    // 5. Доля пар с общим предком (выборка)
    int sample_size = round(params.N * 0.1); // fsample = 0.1
    double common_ancestor = 0.0;
    int total_pairs = 0;

    for (int pair_idx = 0; pair_idx < sample_size; ++pair_idx) {
        int i1 = RNG.RandomInt(0, params.N - 1);
        int i2 = RNG.RandomInt(0, params.N - 1);

        for (int locus = 0; locus < params.L; ++locus) {
            if (sim_params.A[i1][locus] == sim_params.A[i2][locus]) ++common_ancestor;
        }
        ++total_pairs;
    }
    sim_params.C[t] = common_ancestor / (total_pairs * params.L);
    // 6. Доля локусов, где все особи имеют общего предка
    double common_all = 0.0;
    for (int locus = 0; locus < params.L; ++locus) {
        bool all_same = true;
        int first_ancestor = sim_params.A[0][locus];
        for (int i = 1; i < params.N; ++i) {
            if (sim_params.A[i][locus] != first_ancestor) {
                all_same = false;
                break;
            }
        }
        if (all_same) ++common_all;
    }
    sim_params.C_all[t] = common_all / params.L;
    // 7. Средняя приспособленность
    sim_params.mean_W[t] = mean_w;

    // Вычисляем dist (генетическое разнообразие)
    double dist_sum = 0.0;
    for (int locus = 0; locus < params.L; ++locus) {
        double f = sim_params.f_site[t][locus];
        dist_sum += 2.0 * f * (1.0 - f);
    }
    sim_params.dist_over_L[t] = dist_sum / params.L; // dist/L

    // Вычисляем теоретическую частоту f1site
    sim_params.f1_site[t] = params.f0 / (params.f0 + (1.0 - params.f0) * exp(-params.s0 * t));
}

void SaveFigureData(const StartParams& params, const SimulationParams& sim_params) {
    string dir = params.output_dir + "/" + params.exp_name;
    if (!FILEUTILS::FileExists(dir)) FILEUTILS::CreateDirectory(dir);

    // 1. Сохраняем параметры
    ofstream param_file(dir + "/parameters.txt");
    param_file << "N=" << params.N << "\n";
    param_file << "L=" << params.L << "\n";
    param_file << "s0=" << params.s0 << "\n";
    param_file << "r=" << params.r << "\n";
    param_file << "f0=" << params.f0 << "\n";
    param_file << "muL=" << params.muL << "\n";
    param_file << "tf=" << params.tf << "\n";
    param_file << "M=" << params.M << "\n";
    param_file << "distribution=";
    switch (params.distribution) {
    case DistributionType::CONSTANT: param_file << "const"; break;
    case DistributionType::EXPONENTIAL: param_file << "exponential"; break;
    case DistributionType::HALF_GAUSSIAN: param_file << "halfgaussian"; break;
    }
    param_file << "\n";
    param_file << "V_num=" << sim_params.V_num << "\n";
    param_file << "V_an=" << sim_params.V_an << "\n";
    param_file.close();

    // 2. Сохраняем данные для подграфика 2 (средние величины)
    vector<vector<double>> mean_data(sim_params.T.size(), vector<double>(8));
    for (size_t i = 0; i < sim_params.T.size(); ++i) {
        mean_data[i][0] = sim_params.T[i];                      // время
        mean_data[i][1] = sim_params.k_av[i] / params.L;        // k_ср/L
        mean_data[i][2] = sqrt(sim_params.V_ark[i]) / params.L; // sqrt(Vark)/L
        mean_data[i][3] = sim_params.dist_over_L[i];            // dist/L
        mean_data[i][4] = sim_params.C[i];                      // C
        mean_data[i][5] = sim_params.f_survive[i];              // f_survive
        mean_data[i][6] = sim_params.C_all[i];                  // Call
        mean_data[i][7] = sim_params.f1_site[i];                // f1site
    }

    vector<string> headers = { "t", "kav_over_L", "sqrtVark_over_L", "dist_over_L", "C", "fsurvive", "Call", "f1site" };
    FILEUTILS::WriteCSV(dir + "/mean_observables.csv", mean_data, headers);

    // 3. Сохраняем данные для подграфика 3 (частоты по локусам)
    // Сохраняем f_site целиком
    ofstream fsite_file(dir + "/fsite_matrix.csv");
    fsite_file << "time,locus,frequency\n";
    for (size_t t = 0; t < sim_params.T.size(); ++t) {
        for (int locus = 0; locus < params.L; ++locus) {
            fsite_file << sim_params.T[t] << "," << locus << "," << sim_params.f_site[t][locus] << "\n";
        }
    }
    fsite_file.close();

    // 4. Сохраняем гистограммы приспособленностей
    FILEUTILS::WriteHistogramData(dir + "/fitness_histograms.csv", sim_params.hist_times, sim_params.fitness_hist_xx, sim_params.fitness_hist_nn);
    // 5. Сохраняем гистограммы частот аллелей
    FILEUTILS::WriteHistogramData(dir + "/frequency_histograms.csv", sim_params.freq_hist_times, sim_params.freq_hist_xx, sim_params.freq_hist_nn);
    // 6. Сохраняем средние приспособленности
    FILEUTILS::WriteVector(dir + "/mean_fitness.csv", sim_params.mean_W);
    // 7. Сохраняем теоретические частоты
    FILEUTILS::WriteVector(dir + "/f1site_theoretical.csv", sim_params.f1_site);

    cout << "Data saved in directory: " << dir << endl;
}

void SimulationStep(const StartParams& params, SimulationParams& sim_params, const int& t) {
    // 1. Мутации    
    Mutation(params, sim_params);
    // 2. Начальные метки родителей: P = (1:N)'*ones(1,L)
    // Каждый локус имеет родителя - саму особь
    sim_params.P.resize(params.N, vector<int>(params.L));
    for (int i = 0; i < params.N; ++i) {
        for (int j = 0; j < params.L; ++j) {
            sim_params.P[i][j] = i + 1;
        }
    }

    // 2. Начальные метки родителей: P = (1:N)'*ones(1,L)
    // Каждый локус имеет родителя - саму особь
    sim_params.P.resize(params.N, vector<int>(params.L));
    // Для каждого вектора в P заполняем его значениями от 1 до N
    for_each(execution::par, sim_params.P.begin(), sim_params.P.end(), [i = 0, L = params.L](auto& row) mutable { iota(row.begin(), row.end(), ++i); });

    // 3. Расчет приспособленностей: w = K * s'
    vector<double> w(params.N);
    for (int i = 0; i < params.N; ++i) {
        double fitness = 0.0;
        for (int j = 0; j < params.L; ++j) {
            fitness += sim_params.K[i][j] * sim_params.s[j];
        }
        w[i] = fitness;
    }
    // 4. Расчет среднего числа потомков: nprogav = exp(w) / mean(exp(w))
    vector<double> exp_w(params.N);
    double sum_exp_w = 0.0;

    for (int i = 0; i < params.N; ++i) {
        exp_w[i] = exp(w[i]);
        sum_exp_w += exp_w[i];
    }

    double mean_exp_w = sum_exp_w / params.N;
    vector<double> n_prog_av(params.N);
    for (int i = 0; i < params.N; ++i) {
        n_prog_av[i] = exp_w[i] / mean_exp_w;
    }
    // 5. Метод "сломанной палки": расчет кумулятивных сумм
    vector<int> n_prog = BrokenStick(params, n_prog_av);
    // 6. Обновление популяции (отбор)
    UpdatePopulation(params, sim_params, n_prog);
    // 7. Рекомбинация
    if (params.r > 0) Recombination(params, sim_params);
    // 8. Запись наблюдаемых величин
    RecordObservables(params, sim_params, t, w);

    if (sim_params.t_int > 0 && t % sim_params.t_int == 0) {
        // Гистограмма приспособленностей
        vector<double> w_current(params.N);
        for (int i = 0; i < params.N; ++i) {
            w_current[i] = w[i];
        }
        auto hist = ComputeHistogram(w_current);
        sim_params.fitness_hist_xx.push_back(hist.first);
        sim_params.fitness_hist_nn.push_back(hist.second);
        sim_params.hist_times.push_back(t);

        // Гистограмма частот аллелей
        vector<double> freqs_current(params.L);
        for (int locus = 0; locus < params.L; ++locus) {
            freqs_current[locus] = sim_params.f_site[t][locus];
        }
        auto freq_hist = ComputeHistogram(freqs_current);
        sim_params.freq_hist_xx.push_back(freq_hist.first);
        sim_params.freq_hist_nn.push_back(freq_hist.second);
        sim_params.freq_hist_times.push_back(t);
    }

    // 9. Сохранение меток родителей для первого и последнего локусов
    for (int i = 0; i < params.N; ++i) {
        sim_params.P1[i][t] = sim_params.P[i][0];
        sim_params.PL[i][t] = sim_params.P[i][params.L - 1];
    }
    // 10. Сохранение приспособленностей
    for (int i = 0; i < params.N; ++i) {
        sim_params.W[i][t] = w[i];
    }
}

void RunSimulation(const StartParams& params) {
    InitDirectory(params);
    SimulationParams simulation_params = InitStartParams(params);

    for (const auto& t : simulation_params.T) {
        SimulationStep(params, simulation_params, t);
        if (t % (params.tf / 3) == 0) cout << "Generation " << t << " from " << params.tf << " done\n";
    }
    ComputeNumericalVelocity(params, simulation_params);
    cout << "V_an = " << simulation_params.V_an << "\n";
    cout << "V_num = " << simulation_params.V_num << "\n";

    // Сохраняем данные для графиков
    SaveFigureData(params, simulation_params);
}

int main() {
    cout << setprecision(3);
    cout << "Monte Carlo Population Genetics Simulation\n";
    cout << "=========================================\n";

    // Основной эксперимент с постоянным распределением отбора
    StartParams params;
    params.exp_name = "test";
    params.N = 1000;
    params.L = 300;
    params.s0 = 0.1;
    params.r = 0;
    params.f0 = 0;
    params.muL = 0.01;
    params.tf = 150;
    params.M = 3;
    // Запускаем симуляцию
    RunSimulation(params);

    system(("python fig1.py --dir " + params.output_dir + "/" + params.exp_name).c_str());

    return 0;
}