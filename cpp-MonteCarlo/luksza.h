#pragma once
#include "domain.h"
#include "framework.h"
#include "random.h"

class LukszaModel {
public:
    LukszaModel() = default;
    LukszaModel(const StartParams& start_params, const LukszaParams& luksza_params) : L_(start_params.L), mu_(start_params.muL / start_params.L), params_(luksza_params) {
        if (params_.L_ep > L_) {
            cerr << "Error: epitope length exceeds total length\n";
            params_.L_ep = L_;
            L_ne_ = 0;
        }
        L_ne_ = L_ - params_.L_ep;

        cout << "LukszaModel initialized: L = " << L_ << ", L_ep = " << params_.L_ep << ", L_ne = " << L_ne_
            << ", sigma_ne = " << params_.sigma_ne << ", D0 = " << params_.D0
            << ", T_memory = " << params_.T_memory << ", mu = " << mu_ << "\n";
    }
    // Основная функция: вычисление фитнеса для популяции
    // Принимает K напрямую (vector<vector<int>>)
    vector<double> CalculateFitness(const vector<vector<int>>& K, const vector<vector<int>>& A, const vector<vector<int>>& K_prev) {
        int N = K.size();
        if (N == 0) return vector<double>();
        
        // 1. Преобразуем K в строки и подсчитываем частоты
        current_population_ = ConvertNLMapToStrings(K);
        vector<string> prev_population = ConvertNLMapToStrings(K_prev);
        UpdateAncestors(A, prev_population);
        
        // 2. Вычисляем сырой фитнес для каждого уникального гаплотипа
        unordered_map<string, double> raw_fitness_map; // хэш: гаплотип -> сырой фитнес
        unordered_map<string, double> freq_map = current_freq_; // копия частот

        // Для каждого уникального гаплотипа
        for (const auto& [hap, freq] : freq_map) {
            // 2a. Мутационная нагрузка
            double mutational_load = ComputeMutationalLoad(hap);
            // 2b. Иммунная нагрузка (суммируем по всей истории)
            double immune_load = ComputeImmuneLoad(hap);
            // 2c. Сырой фитнес (без нормировки f0)
            double raw_fitness = -mutational_load - immune_load;
            raw_fitness_map[hap] = raw_fitness;
        }

        // 3. ВЫЧИСЛЯЕМ НОРМИРОВОЧНУЮ КОНСТАНТУ f0
           // Согласно формуле (19) из статьи: Σ x_i * exp(f_i) = 1
           // Где f_i = raw_fitness_i + f0
           // Поэтому: Σ x_i * exp(raw_fitness_i + f0) = 1
           // => exp(f0) * Σ x_i * exp(raw_fitness_i) = 1
           // => f0 = -log(Σ x_i * exp(raw_fitness_i))

        double sum_exp_raw = 0.0;
        for (const auto& [hap, freq] : freq_map) {
            sum_exp_raw += freq * exp(raw_fitness_map[hap]);
        }

        // Защита от деления на 0 (если sum_exp_raw очень маленький)
        if (sum_exp_raw <= 0.0) {
            cerr << "LukszaModel Warning: sum_exp_raw <= 0, setting to 1.0\n";
            sum_exp_raw = 1.0;
        }

        double f0 = -log(sum_exp_raw);

        // 4. ВЫЧИСЛЯЕМ ФИНАЛЬНЫЙ ФИТНЕС ДЛЯ КАЖДОГО УНИКАЛЬНОГО ГАПЛОТИПА
        unordered_map<string, double> final_fitness_map;
        for (const auto& [hap, raw_fit] : raw_fitness_map) {
            final_fitness_map[hap] = raw_fit + f0;
        }

        // 5. СОЗДАЕМ ВЕКТОР ФИТНЕСА ДЛЯ ВСЕХ ОСОБЕЙ
        vector<double> fitness(N);
        for (int i = 0; i < N; ++i) {
            const string& hap = current_population_[i];
            fitness[i] = exp(final_fitness_map[hap]);
        }

        current_freq_ = ComputeFreq();
        UpdateHistory();

        return fitness;
    }
    // Очистка истории
    void ClearHistory() { history_.clear(); }

    void PrintInfo() const {
        cout << "=== LukszaModel Info ===\n";
        cout << "L: " << L_ << "\n";
        cout << "L_ep: " << params_.L_ep << "\n";
        cout << "L_ne: " << L_ne_ << "\n";
        cout << "sigma_ne: " << params_.sigma_ne << "\n";
        cout << "D0: " << params_.D0 << "\n";
        cout << "T_memory: " << params_.T_memory << "\n";
        cout << "mu: " << mu_ << "\n";
        cout << "History size: " << history_.size() << "\n";
        cout << "=======================\n";
    }

private:
    typedef unordered_map<string, double> CurFreqCondition;
    
    LukszaParams params_;
    int L_;      // Общая длина гаплотипа
    double mu_;  // Частота мутаций на локус
    int L_ne_;   // Количество неэпитопных позиций
    // История популяции
    deque<CurFreqCondition> history_;
    unordered_map<string, unordered_map<string, int>> ancestor_freqs_;
    unordered_map<string, string> ancestors_map_;

    vector<string> current_population_;
    CurFreqCondition current_freq_;

    // Обновление истории после завершения поколения
    void UpdateHistory() {
        if (current_population_.size() == 0) return;
        // Добавляем в историю
        history_.push_back(current_freq_);
        // Ограничиваем глубину истории
        if (history_.size() > params_.T_memory) history_.pop_front();
    }

    void UpdateAncestors(const vector<vector<int>>& A, const vector<string>& prev_population) {
        ancestors_map_.clear();
        // Проверяем, есть ли предыдущее поколение
        if (prev_population.empty()) {
            // Нет предков для первого поколения
            for (const string& hap : current_population_)
                ancestors_map_[hap] = "";  // пустая строка = нет предка
            return;
        }
        // Для каждой особи находим доминантного предка
        for (size_t i = 0; i < current_population_.size(); ++i) {
            string current_hap = current_population_[i];
            // Собираем всех предков по локусам
            unordered_map<int, int> ancestor_counts;  // индекс предка → количество локусов
            for (int locus = 0; locus < L_; ++locus) {
                int ancestor_idx = A[i][locus];
                if (ancestor_idx >= 0 && ancestor_idx < static_cast<int>(prev_population.size())) ++ancestor_counts[ancestor_idx];
            }
            // Находим наиболее частого предка
            int dominant_idx = -1;
            int max_count = 0;
            for (const auto& [idx, count] : ancestor_counts) {
                if (count > max_count) {
                    max_count = count;
                    dominant_idx = idx;
                }
            }
            // Сохраняем гаплотип предка
            if (dominant_idx >= 0) ancestors_map_[current_hap] = prev_population[dominant_idx];
            else ancestors_map_[current_hap] = "";
        }
    }

    unordered_map<string, double> ComputeFreq() {
        unordered_map<string, int> counts;
        for (const auto& hap : current_population_) {
            ++counts[hap];
        }

        unordered_map<string, double> current_freq;
        int N = current_population_.size();
        for (const auto& [hap, count] : counts) {
            current_freq[hap] = static_cast<double>(count) / N;
        }
        return current_freq;
    }

    // Преобразование vector<vector<int>> в vector<string>
    vector<string> ConvertNLMapToStrings(const vector<vector<int>>& K) {
        int N = K.size();
        if (N == 0) return vector<string>();

        vector<string> population(N);
        for (int i = 0; i < N; ++i) {
            string hap(L_, '0');
            for (int j = 0; j < L_; ++j) {
                hap[j] = (K[i][j] == 1) ? '1' : '0';
            }
            population[i] = hap;
        }
        return population;
    }
    // Вычисление мутационной нагрузки L(a_i)
    double ComputeMutationalLoad(const string& hap) {
        if (!ancestors_map_.count(hap) || ancestors_map_[hap].empty()) return 0.0;  // нет предка
        int d_ne = HammingDistanceNonEpitope(hap, ancestors_map_[hap]);
        return params_.sigma_ne * d_ne;
    }
    // Вычисление нагрузки перекрёстного иммунитета
    double ComputeImmuneLoad(const string& hap) {
        double immune_load = 0.0;
        for (const auto& hist_block : history_) {
            for (const auto& [hist_hap, hist_freq] : hist_block) {
                int d_ep = HammingDistanceEpitope(hap, hist_hap);
                double C = exp(-d_ep / params_.D0);
                immune_load += hist_freq * C;
            }
        }
        return immune_load;
    }

    // Расстояние по эпитопным позициям (первые L_ep позиций)
    int HammingDistanceEpitope(const string& a, const string& b) {
        int dist = 0;
        for (int i = 0; i < params_.L_ep; ++i) {
            if (a[i] != b[i]) ++dist;
        }
        return dist;
    }
    // Расстояние по неэпитопным позициям (оставшиеся L_ne позиций)
    int HammingDistanceNonEpitope(const string& a, const string& b) {
        int dist = 0;
        for (int i = params_.L_ep; i < L_; ++i) {
            if (a[i] != b[i]) ++dist;
        }
        return dist;
    }
};