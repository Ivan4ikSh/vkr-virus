#pragma once
//uint32_t SEED = RANDOM::Hash();
//RandomGenerator RNG(SEED);
RandomGenerator RNG(1);

void InitDirectory(const StartParams& params) {
    FILEUTILS::CreateDirectory(params.output_dir);
    FILEUTILS::CreateDirectory(params.output_dir + "/" + params.exp_name);
    FILEUTILS::CreateDirectory(params.plot_dir);
}

template<typename T>
double Mean(const vector<T>& v) {
    return accumulate(v.begin(), v.end(), 0.0) / v.size();
}

template<typename T>
double StdDev(const vector<T>& v) {
    if (v.size() <= 1) return 0.0;
    double m = Mean(v);
    double acc = 0.0;
    for (double x : v) acc += (x - m) * (x - m);
    return sqrt(acc / (v.size() - 1));
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

    for (double i = t_start; i < V_num_vec.size(); ++i) {
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
        cout << "WARNING! Velocity may be computed with error\n";
        cout << "  Required: log(arg_log2) >= 3\n";
        cout << "  Actual:   log(arg_log2) = " << log_arg_log2 << "\n";
        cout << "  Required: log(N_sqrt_sUb) >= 3\n";
        cout << "  Actual:   log(N_sqrt_sUb) = " << log_N_sqrt_sUb << "\n";
        cout << "  s/Ub = " << s_Ub_ratio << " (should be >> 1)\n";
    }

    return V;
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

vector<double> ComputeFitness(const StartParams& params, SimulationParams& sim_params) {
    vector<double> w(params.N);
  
    for (int i = 0; i < params.N; ++i) {
        double fitness = 0.0;
        for (int j = 0; j < params.L; ++j) {
            fitness += sim_params.K[i][j] * sim_params.s[j];
        }
        w[i] = fitness;
    }
    return w;
}

vector<double> ComputeDescNumAv(const StartParams& params, SimulationParams& sim_params, const vector<double>& w) {
    int N = params.N;
    vector<double> n_prog_av(N);

    // 1. Находим максимальное w для масштабирования
    double max_w = w[0];
    for (int i = 1; i < N; ++i) {
        if (w[i] > max_w) max_w = w[i];
    }

    // 2. Масштабируем, чтобы избежать переполнения exp(w)
    double sum_exp = 0.0;
    vector<double> exp_w(N);

    for (int i = 0; i < N; ++i) {
        exp_w[i] = exp(w[i] - max_w);  // Вычитаем максимум для стабильности
        sum_exp += exp_w[i];
    }

    double mean_exp = sum_exp / N;

    // 3. Вычисляем среднее число потомков
    for (int i = 0; i < N; ++i) {
        n_prog_av[i] = exp_w[i] / mean_exp;
    }

    return n_prog_av;
}

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


void Mutation(const StartParams& params, SimulationParams& simulation_params) {
    if (simulation_params.mu <= 0) return;

    int total_mutations = 0;
    for (int i = 0; i < params.N; ++i) {
        for (int j = 0; j < params.L; ++j) {
            if (RNG.Random() < simulation_params.mu) 
                simulation_params.K[i][j] ^= 1; // XOR для инвертирования бита
        }
    }
}