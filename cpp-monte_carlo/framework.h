#pragma once
uint32_t SEED = RANDOM::Hash();
RandomGenerator RNG(SEED);

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

// Правая часть уравнения (51)
double RightSide51(double V, double s, double Ub) {
    double log_term = log(V / (exp(1.0) * Ub));
    double term1 = V / (2.0 * s) * (log_term * log_term + 1.0);
    double term2 = 0.5 * log((s * s * s * Ub) / (V * V * log(V / Ub)));
    return term1 - term2;
}

// Правая часть уравнения (52)
double RightSide52(double V, double s, double Ub) {
    double log_term = log(V / (exp(1.0) * Ub));
    double term1 = V / (2.0 * s) * (log_term * log_term + 1.0);
    double term2 = 0.5 * log((s * s * Ub) / (V * log(V / Ub)));
    return term1 - term2;
}

// Метод бисекции для решения уравнения f(V) = target
double Bisection(double V_low, double V_high, double target, double (*func)(double, double, double), double s, double Ub, double tol = 1e-8, int max_iter = 1000) {
    double f_low = func(V_low, s, Ub) - target;
    double f_high = func(V_high, s, Ub) - target;

    // Проверка, что функция меняет знак на отрезке
    if (f_low * f_high > 0) {
        // Если оба значения положительны, возможно, V_low слишком велико
        if (f_low > 0 && f_high > 0) return V_low;
        // Если оба отрицательны, V_high слишком мало
        return V_high;
    }

    for (int iter = 0; iter < max_iter; ++iter) {
        double V_mid = (V_low + V_high) / 2.0;
        if (V_high - V_low < tol) return V_mid;
        double f_mid = func(V_mid, s, Ub) - target;
        
        if (f_mid == 0.0) return V_mid;
        else if (f_low * f_mid < 0) {
            V_high = V_mid;
            f_high = f_mid;
        }
        else {
            V_low = V_mid;
            f_low = f_mid;
        }
    }

    return (V_low + V_high) / 2.0;
}

double ComputeAnalyticalVelocity(const StartParams& params) {
    double Ub = params.muL * (1.0 - params.f0);
    double s = params.s0;
    double N = static_cast<double>(params.N);
    double target = log(N);

    // Проверка применимости теории
    if (Ub <= 0 || s <= 0) {
        cerr << "ERROR: Ub or s must be positive" << endl;
        return 0.0;
    }
    // Определяем границы поиска V
    double V_min = Ub * 1.001;  // V должно быть > Ub для положительного логарифма
    double V_max = 1000.0 * s;  // Начальная верхняя граница
    // Увеличиваем V_max до тех пор, пока правая часть (51) не станет больше target
    while (RightSide51(V_max, s, Ub) < target && V_max < 1e20) {
        V_max *= 2.0;
    }

    if (V_max >= 1e20) {
        cerr << "ERROR: Cannot find suitable V_max, N might be too large" << endl;
        return 0.0;
    }
    // Решаем уравнение (51)
    double V1 = Bisection(V_min, V_max, target, RightSide51, s, Ub);
    // Проверяем условие длинного хвоста
    if (V1 * log(V1 / Ub) < s) {
        cout << "WARNING: Condition V*ln(V/Ub) >> s might not be satisfied." << endl;
        cout << "  V*ln(V/Ub) = " << V1 * log(V1 / Ub) << ", s = " << s << endl;
    }

    // Если V1 > s, используем результат (51) для широкого распределения
    if (V1 > s) return V1;
    // Иначе решаем уравнение (52) для узкого распределения
    else {
        // Решаем уравнение (52)
        double V2 = Bisection(V_min, V_max, target, RightSide52, s, Ub);
        // Проверяем условие для формулы (52)
        if (V2 < s && V2 > s / log(V2 / Ub)) return V2;
        else {
            // Если условия не выполняются, используем V1 с предупреждением
            cout << "WARNING: Conditions for formula (52) not satisfied." << endl;
            cout << "  Using V from formula (51): V = " << V1 << endl;
            return V1;
        }
    }
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