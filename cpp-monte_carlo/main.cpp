#include "domain.h"

class RecombinationSimulator {
public:
    RecombinationSimulator(const Args& population_args) : params(population_args) {
        LOG::Log("INIT", "RecombinationSimulator", "constructor called");
    }

    // Основной метод запуска симуляции
    DATA::Output RunRecombTrain() {
        InitializePopulationData();
        DATA::Output output;

        double& s0 = params.adaptive_landscape;
        double& muL = params.mutation_freq;
        double& r = params.recomb_cnt;
        int& M = params.cross_cnt;
        int& N = params.population_cnt;
        int& L = params.locus_cnt;
        int& tf = params.full_time_in_epoch;
        double& f0 = params.start_good_allele_freq;

        output.theoretical_velocity = ComputeTheoreticalVelocity(N, s0, L, f0, muL);

        vector<string> colors = { "r", "g", "b", "m", "k", "r", "g", "b", "m", "k",
                                  "r", "g", "b", "m", "k", "r", "g", "b", "m", "k" };

        for (int t = 0; t < tf; ++t) {
            if (t % data.tint == 0) LOG::Log("PROGRESS", "RunRecombTrain", "generation " + to_string(t) + "/" + to_string(tf));

            ApplyMutations(data.K, data.mu, N, L);

            vector<vector<int>> P(N, vector<int>(L));
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < L; ++j) P[i][j] = i;
            }

            vector<double> w = CalculateFitness(data.K, data.s);
            vector<double> exp_w = CalculateExponentialFitness(w);
            vector<double> nprogav = CalculateProgenyAverages(exp_w);
            vector<int> nprog = CalculateProgenyCounts(nprogav, N);

            UpdatePopulationWithSelection(data.K, data.A, P, data.Knew, data.Anew, data.Pnew, nprog);
            ApplyRecombination(data.K, data.A, P, r, M);
            CalculateAllStatistics(w, P, t);

            output.time_points.push_back(t);
            output.kav.push_back(data.kav[t]);
            output.Vark.push_back(data.Vark[t]);
            output.fsurvive.push_back(data.fsurvive[t]);
            output.C.push_back(data.C[t]);
            output.Call.push_back(data.Call[t]);
            output.meanW.push_back(data.meanW[t]);
            output.fsite.push_back(data.fsite[t]);

            if (t % data.tint == 0) {
                vector<double> w_current = CalculateFitness(data.K, data.s);
                vector<double> w_normalized(w_current.size());
                for (size_t i = 0; i < w_current.size(); ++i) {
                    w_normalized[i] = w_current[i] / s0;
                }

                double min_w = *min_element(w_normalized.begin(), w_normalized.end());
                double max_w = *max_element(w_normalized.begin(), w_normalized.end());
                int bins = 50;
                auto histogram = ComputeHistogram(w_normalized, bins, min_w, max_w);

                output.fitness_histograms.push_back(histogram);
                output.histogram_times.push_back(t);

                int color_index = (t / data.tint) % colors.size();
                output.histogram_colors.push_back(colors[color_index]);
            }
        }

        return output;
    }

    // Функции для сохранения результатов
    void SaveResults(const DATA::Output& results, const string& base_name) {
        LOG::Log("SAVE", "SaveResults", "saving all results");

        // Сохраняем начальную матрицу локусов
        SaveInitialLocusMatrix("data/" + base_name);

        SAVE::Parameters(params, results.theoretical_velocity, "data/" + base_name + "_params.csv");
        SAVE::Statistics(results, params, "data/" + base_name + "_stats.csv");
        SAVE::FitnessHistograms(results, "data/" + base_name + "_histograms.csv");
        SAVE::AlleleFrequencies(results, "data/" + base_name + "_frequencies.csv");
    }

    // Функции для построения графиков
    void PlotFitnessHistograms(const string& base_name) {
        LOG::Log("PLOT", "PlotFitnessHistograms", "calling Python script for fitness histograms");
        string command = "python scripts/plot_fig1.py " + base_name;
        system(command.c_str());
    }

    void PlotStatistics(const string& base_name) {
        LOG::Log("PLOT", "PlotStatistics", "calling Python script for statistics");
        string command = "python scripts/plot_fig4.py " + base_name;
        system(command.c_str());
    }

private:
    Args params;
    DATA::Population data;

    void SaveInitialLocusMatrix(const string& base_name) {
        string filename = base_name + "_initial_locus_matrix.bin";
        ofstream file(filename, ios::binary);

        if (!file.is_open()) {
            LOG::Log("ERROR", "SaveInitialLocusMatrixOptimized", "Cannot open file: " + filename);
            return;
        }

        int N = data.K.size();
        int L = (N > 0) ? data.K[0].size() : 0;

        file.write(reinterpret_cast<const char*>(&N), sizeof(N));
        file.write(reinterpret_cast<const char*>(&L), sizeof(L));

        int total_bits = N * L;
        int total_bytes = (total_bits + 7) / 8;
        vector<unsigned char> buffer(total_bytes, 0);

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < L; ++j) {
                int bit_position = i * L + j;
                int byte_index = bit_position / 8;
                int bit_index = bit_position % 8;

                if (data.K[i][j] == 1) buffer[byte_index] |= (1 << bit_index);
            }
        }

        file.write(reinterpret_cast<const char*>(buffer.data()), buffer.size());

        file.close();
        LOG::Log("SAVE", "SaveInitialLocusMatrixOptimized", "Initial locus matrix saved in optimized binary format to: " + filename);
    }
    // Вспомогательные методы
    double ComputeTheoreticalVelocity(int N, double s, int L, double f0, double muL) {
        double Ub = muL * (1.0 - f0);
        if (Ub > 0 || s > 0) {
            double N_sqrt_sUb = N * sqrt(s * Ub);
            if (N_sqrt_sUb > 1.0) {
                double log_N_sqrt_sUb = log(N_sqrt_sUb);
                double arg_log2 = s / Ub * log_N_sqrt_sUb;
                if (arg_log2 > 1.0) {
                    double V = 2.0 * s * log_N_sqrt_sUb / (log(arg_log2) * log(arg_log2));
                    if (V >= 0.0 || V <= s * L) return V;
                }
            }
        }
        return 0.0;
    }

    vector<double> ComputeHistogram(const vector<double>& data, int bins, double min_val, double max_val) {
        vector<int> counts(bins, 0);
        double bin_width = (max_val - min_val) / bins;

        for (double value : data) {
            int bin_index = min(bins - 1, max(0, static_cast<int>((value - min_val) / bin_width)));
            ++counts[bin_index];
        }

        vector<double> result;
        for (int i = 0; i < bins; ++i) {
            double bin_center = min_val + (i + 0.5) * bin_width;
            double frequency = static_cast<double>(counts[i]) / data.size();
            result.push_back(bin_center);
            result.push_back(frequency);
        }
        return result;
    }

    void ApplyConstantDistribution(vector<double>& s, double s0) {
        for (auto& s_i : s) s_i = s0;
    }

    void ApplyExponentialDistribution(vector<double>& s, double s0) {
        for (auto& s_i : s) s_i = -s0 * log(RANDOM::ScaleToRange01(rand()));
    }

    void ApplyHalfGaussianDistribution(vector<double>& s, double s0) {
        for (auto& s_i : s) s_i = s0 * sqrt(CONST::PI / 2.0) * abs(RANDOM::ScaleToRange01(rand()));
    }

    void DistributionSelection(const Distribution& distribution_s, vector<double>& s, const double& s0) {
        switch (distribution_s) {
        case Distribution::CONSTANT: ApplyConstantDistribution(s, s0); break;
        case Distribution::EXPONENTIAL: ApplyExponentialDistribution(s, s0); break;
        case Distribution::HALF_GAUSSIAN: ApplyHalfGaussianDistribution(s, s0); break;
        default: break;
        }
    }

    void InitializeParameters() {
        data.mu = params.mutation_freq / params.locus_cnt;
        data.tint = round(static_cast<double>(params.full_time_in_epoch) / 10.0);
        data.fsample = 0.1;
    }

    void InitializeAncestryMatrix(int N, int L) {
        data.A.assign(N, {});
        for (int i = 0; i < N; ++i) data.A[i].assign(L, i);
    }

    void InitializeDNA(int N, int L, double f0) {
        data.K.resize(N, vector<int>(L, 0));
        if (f0 != 0) {
            for (auto& row : data.K) {
                std::generate(row.begin(), row.end(), [&]() {
                    return RANDOM::ScaleToRange01(rand()) < f0 ? 1 : 0;
                    });
            }
        }
    }

    void ReserveMemoryForTemporalData(int N, int L, int time_steps) {
        data.W.resize(N, vector<double>(time_steps, 0.0));
        data.P1.resize(N, vector<double>(time_steps, 0.0));
        data.PL.resize(N, vector<double>(time_steps, 0.0));
        data.fsite.resize(time_steps, vector<double>(L, 0.0));
    }

    void ReserveMemoryForStatistics(int time_steps) {
        data.kav.resize(time_steps, 0.0);
        data.Vark.resize(time_steps, 0.0);
        data.fsurvive.resize(time_steps, 0.0);
        data.C.resize(time_steps, 0.0);
        data.Call.resize(time_steps, 0.0);
        data.meanW.resize(time_steps, 0.0);
    }

    void ReserveMemoryForNewGenerations(int N, int L) {
        data.Knew.resize(N, vector<int>(L, 0));
        data.Anew.resize(N, vector<int>(L, 0));
        data.Pnew.resize(N, vector<int>(L, 0));
    }

    void InitializePopulationData() {
        const double& s0 = params.adaptive_landscape;
        const int& L = params.locus_cnt;
        const int& N = params.population_cnt;
        const int& tf = params.full_time_in_epoch;
        const double& f0 = params.start_good_allele_freq;

        InitializeParameters();
        data.s.resize(L);
        DistributionSelection(params.dist_s, data.s, s0);
        InitializeAncestryMatrix(N, L);
        InitializeDNA(N, L, f0);

        int time_steps = tf;
        ReserveMemoryForTemporalData(N, L, time_steps);
        ReserveMemoryForStatistics(time_steps);
        ReserveMemoryForNewGenerations(N, L);

        LOG::Log("INIT", "InitializePopulationData", "completed");
    }

    void ApplyMutations(vector<vector<int>>& K, double mu, int N, int L) {
        if (mu > 0) {
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < L; ++j) {
                    if (RANDOM::ScaleToRange01(rand()) < mu) K[i][j] = K[i][j] ^ 1;
                }
            }
        }
    }

    vector<double> CalculateFitness(const vector<vector<int>>& K, const vector<double>& s) {
        int N = K.size();
        int L = K[0].size();
        vector<double> w(N, 0.0);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < L; ++j) {
                w[i] += K[i][j] * s[j];
            }
        }
        return w;
    }

    vector<double> CalculateExponentialFitness(const vector<double>& w) {
        vector<double> exp_w;
        for (auto& w_i : w) exp_w.push_back(exp(w_i));
        return exp_w;
    }

    vector<double> CalculateProgenyAverages(const vector<double>& exp_w) {
        double sum = 0.0;
        for (double val : exp_w) sum += val;
        double mean_exp_w = sum / exp_w.size();

        vector<double> nprogav(exp_w.size());
        for (size_t i = 0; i < exp_w.size(); ++i) {
            nprogav[i] = exp_w[i] / mean_exp_w;
        }
        return nprogav;
    }

    vector<int> CalculateProgenyCounts(const vector<double>& nprogav, int N) {
        vector<double> b2(N);
        b2[0] = nprogav[0];
        for (int i = 1; i < N; ++i) {
            b2[i] = b2[i - 1] + nprogav[i];
        }

        vector<int> nprog(N, 0);
        for (int i = 0; i < N; ++i) {
            double X = RANDOM::ScaleToRange01(rand()) * b2[N - 1];
            for (int j = 0; j < N; ++j) {
                if ((j == 0 && X < b2[0]) || (j > 0 && X >= b2[j - 1] && X < b2[j])) {
                    ++nprog[j];
                    break;
                }
            }
        }
        return nprog;
    }

    void UpdatePopulationWithSelection(vector<vector<int>>& K, vector<vector<int>>& A, vector<vector<int>>& P,
        vector<vector<int>>& Knew, vector<vector<int>>& Anew, vector<vector<int>>& Pnew,
        const vector<int>& nprog) {
        int N = nprog.size();
        int current_row = 0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < nprog[i]; ++j) {
                Knew[current_row] = K[i];
                Anew[current_row] = A[i];
                Pnew[current_row] = P[i];
                ++current_row;
            }
        }
        K = Knew;
        A = Anew;
        P = Pnew;
    }

    vector<bool> GenerateCrossoverMask(int L, int M) {
        vector<bool> use_parent1(L, true);
        bool current_parent = true;
        for (int locus = 0; locus < L; locus++) {
            if (RANDOM::ScaleToRange01(rand()) < static_cast<double>(M) / L) current_parent = !current_parent;
            use_parent1[locus] = current_parent;
        }
        return use_parent1;
    }

    void PerformRecombination(vector<vector<int>>& K, vector<vector<int>>& A, vector<vector<int>>& P, int parent1, int parent2, const vector<bool>& use_parent1, int L) {
        vector<int> recombinant_DNA(L);
        vector<int> recombinant_ancestor(L);
        vector<int> recombinant_parent(L);

        for (int locus = 0; locus < L; locus++) {
            if (use_parent1[locus]) {
                recombinant_DNA[locus] = K[parent1][locus];
                recombinant_ancestor[locus] = A[parent1][locus];
                recombinant_parent[locus] = P[parent1][locus];
            }
            else {
                recombinant_DNA[locus] = K[parent2][locus];
                recombinant_ancestor[locus] = A[parent2][locus];
                recombinant_parent[locus] = P[parent2][locus];
            }
        }

        if (rand() % 2 == 0) {
            K[parent1] = recombinant_DNA;
            A[parent1] = recombinant_ancestor;
            P[parent1] = recombinant_parent;
        }
        else {
            K[parent2] = recombinant_DNA;
            A[parent2] = recombinant_ancestor;
            P[parent2] = recombinant_parent;
        }
    }

    void ApplyRecombination(vector<vector<int>>& K, vector<vector<int>>& A, vector<vector<int>>& P, double r, int M) {
        int N = P.size();
        int L = P[0].size();
        int npairs = static_cast<int>(round(r * N / 2.0));

        for (int pair = 0; pair < npairs; ++pair) {
            int parent1 = rand() % N;
            int parent2 = rand() % N;
            vector<bool> crossover_mask = GenerateCrossoverMask(L, M);
            PerformRecombination(K, A, P, parent1, parent2, crossover_mask, L);
        }
    }

    void CalculateAlleleFrequencies(vector<vector<double>>& fsite, const vector<vector<int>>& K, int t) {
        int N = K.size();
        int L = K[0].size();

        for (int j = 0; j < L; ++j) {
            double sum = 0.0;
            for (int i = 0; i < N; ++i) {
                sum += K[i][j];
            }
            fsite[t][j] = sum / N;
        }
    }

    void CalculateAverageAlleles(vector<double>& kav, const vector<vector<int>>& K, int t) {
        int N = K.size();
        int L = K[0].size();

        double total_K = 0.0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < L; ++j) {
                total_K += K[i][j];
            }
        }
        kav[t] = total_K / N;
    }

    void CalculateFitnessVariance(vector<double>& Vark, const vector<double>& w, int t, double s0) {
        double meanW_val = 0.0;
        for (double w_val : w) meanW_val += w_val;
        meanW_val /= w.size();

        double varianceW = 0.0;
        for (double w_val : w) varianceW += (w_val - meanW_val) * (w_val - meanW_val);
        varianceW /= w.size();

        Vark[t] = varianceW / (s0 * s0);
    }

    void CalculateSurvivalFraction(vector<double>& fsurvive, const vector<vector<int>>& K, int t) {
        int N = K.size();
        int L = K[0].size();

        int survive_cnt = 0;
        for (int i = 0; i < N; ++i) {
            bool has_allele = false;
            for (int j = 0; j < L; ++j) {
                if (K[i][j] != 0) {
                    has_allele = true;
                    break;
                }
            }
            if (has_allele) survive_cnt++;
        }
        fsurvive[t] = static_cast<double>(survive_cnt) / N;
    }

    void CalculateAncestorCorrelation(vector<double>& C, const vector<vector<int>>& A, int t, double fsample) {
        int N = A.size();
        int L = A[0].size();
        int xx = static_cast<int>(round(N * fsample));

        double total_matches = 0.0;
        for (int k = 0; k < xx; ++k) {
            int idx1 = rand() % N;
            int idx2 = rand() % N;
            int matches = 0;
            for (int j = 0; j < L; ++j) {
                if (A[idx1][j] == A[idx2][j]) matches++;
            }
            total_matches += static_cast<double>(matches);
        }
        C[t] = total_matches / (xx * L);
    }

    void CalculateLocusVariation(vector<double>& Call, const vector<vector<int>>& A, int t) {
        int N = A.size();
        int L = A[0].size();

        double no_variation_cnt = 0.0;
        for (int j = 0; j < L; ++j) {
            bool all_same = true;
            int first_val = A[0][j];
            for (int i = 1; i < N; ++i) {
                if (A[i][j] != first_val) {
                    all_same = false;
                    break;
                }
            }
            if (all_same) no_variation_cnt += 1.0;
        }
        Call[t] = no_variation_cnt / L;
    }

    void StoreFitnessData(vector<vector<double>>& W, vector<double>& meanW, const vector<double>& w, int t) {
        double meanW_val = 0.0;
        for (size_t i = 0; i < w.size(); ++i) {
            W[i][t] = w[i];
            meanW_val += w[i];
        }
        meanW_val /= w.size();
        meanW[t] = meanW_val;
    }

    void StoreParentLabels(vector<vector<double>>& P1, vector<vector<double>>& PL, const vector<vector<int>>& P, int t) {
        int N = P.size();
        for (int i = 0; i < N; ++i) {
            P1[i][t] = P[i][0];
            PL[i][t] = P[i][P[i].size() - 1];
        }
    }

    void CalculateAllStatistics(const vector<double>& w, vector<vector<int>>& P, int t) {
        CalculateAlleleFrequencies(data.fsite, data.K, t);
        CalculateAverageAlleles(data.kav, data.K, t);
        CalculateFitnessVariance(data.Vark, w, t, data.s[0]);
        CalculateSurvivalFraction(data.fsurvive, data.K, t);
        CalculateAncestorCorrelation(data.C, data.A, t, data.fsample);
        CalculateLocusVariation(data.Call, data.A, t);
        StoreFitnessData(data.W, data.meanW, w, t);
        StoreParentLabels(data.P1, data.PL, P, t);
    }
};

int main() {
    LOG::Log("MAIN", "main", "program started");

    Args monte_carlo_params;
    monte_carlo_params.dist_s = Distribution::CONSTANT;
    monte_carlo_params.recomb_cnt = 0;
    monte_carlo_params.cross_cnt = 3;
    monte_carlo_params.adaptive_landscape = 0.1;
    monte_carlo_params.locus_cnt = 300;
    monte_carlo_params.population_cnt = 1000;
    monte_carlo_params.full_time_in_epoch = 500;
    monte_carlo_params.start_good_allele_freq = 0;
    monte_carlo_params.mutation_freq = 0.1;
    monte_carlo_params.exp_name = "experiment_1";

    RecombinationSimulator simulator(monte_carlo_params);
    DATA::Output results = simulator.RunRecombTrain();

    string base_name = "results_" + monte_carlo_params.exp_name;

    simulator.SaveResults(results, base_name);

    simulator.PlotFitnessHistograms(base_name);
    simulator.PlotStatistics(base_name);

    LOG::Log("MAIN", "main", "program completed");
    return 0;
}