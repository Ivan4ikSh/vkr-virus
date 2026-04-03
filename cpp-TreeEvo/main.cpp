#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

using namespace std;
namespace fs = filesystem;

namespace CONST {
    const double EPS = 1e-9;
    const std::string OUTPUT = "out";
    const std::string DATA = "data";
}

std::string NAME = "test";

// ================== Model Parameters ==================
struct ModelParams {
    int seed = 1;
    int iter = 100;
    double R0 = 1.8;
    double Ub = 1e-3;
    double a = 14.0;
    int64_t N = 1e8;
    double dt = 0.5;
    int M = 5000;
    double T0 = 500.0;
    int max_nodes = 8000;
};

// ================== Tree Topology ==================
struct Node {
    int id;
    int depth;
    int creation_step;
    double i;
    double r;
    bool is_infected;

    Node* parent;
    Node* l_child;
    Node* r_child;

    Node(int node_id, int node_depth, int step_born, double initial_i, double initial_r, Node* parent_node = nullptr)
        : id(node_id), depth(node_depth), creation_step(step_born), i(initial_i), r(initial_r), is_infected(initial_i > CONST::EPS), parent(parent_node), l_child(nullptr), r_child(nullptr) {}
};

class TreeTopology {
public:
    TreeTopology() {
        root_ = new Node(0, 0, 0, 0.0, 1.0 - 1e-4);
        all_nodes_.push_back(root_);
        node_count_ = 1;

        Node* curr = root_;
        for (int i = 1; i <= 14; ++i) {
            Node* child = new Node(node_count_++, i, 0, 0.0, 0.0, curr);
            all_nodes_.push_back(child);
            curr->l_child = child;
            curr = child;
        }

        curr->i = 1e-4;
        curr->is_infected = true;
    }

    ~TreeTopology() {
        for (Node* node : all_nodes_) {
            delete node;
        }
    }

    Node* AddChild(Node* parent, int step_born, double initial_i, double initial_r) {
        Node* child = new Node(node_count_++, parent->depth + 1, step_born, initial_i, initial_r, parent);
        all_nodes_.push_back(child);

        if (!parent->l_child) parent->l_child = child;
        else if (!parent->r_child) parent->r_child = child;

        return child;
    }

    const std::vector<Node*>& GetAllNodes() const { return all_nodes_; }
    Node* GetRoot() const { return root_; }
    int GetNodeCount() const { return node_count_; }
private:
    Node* root_;
    int node_count_;
    std::vector<Node*> all_nodes_;
};

// ================== Simulator ==================
class TreeEvolutionSimulator {
public:
    TreeEvolutionSimulator(const ModelParams& params) : params_(params), rng_(params.seed) {
        int max_possible_depth = params_.M + 2000;
        K_table_.resize(max_possible_depth, 0.0);
        for (int d = 0; d < max_possible_depth; ++d) {
            if (d == 0) K_table_[d] = 0.0;
            else K_table_[d] = 1.0 - std::exp(-static_cast<double>(d) / params_.a);
        }
    }

    void Run() {
        Node* root = tree_.GetRoot();
        std::cout << "Starting simulation...\n";

        int T0_step = static_cast<int>(std::round(params_.T0 / params_.dt));

        for (int step = 1; step <= params_.M; ++step) {
            Step(step);

            if (step == T0_step) {
                double sum_pos = 0.0;
                double sI = 0.0;
                for (Node* node : tree_.GetAllNodes()) {
                    sum_pos += node->i * node->depth;
                    sI += node->i;
                }
                if (sI > 0) meanI0_ = sum_pos / sI;
            }

            if (step % params_.iter == 0) std::cout << "Step: " << step << " | Nodes: " << tree_.GetNodeCount() << "\n";
            if (step == params_.M / 4 || step == params_.M / 2 || step == 3 * params_.M / 4 || step == params_.M) SaveSnapshot(step);
        }

        double sum_pos = 0.0;
        double sI = 0.0;
        for (Node* node : tree_.GetAllNodes()) {
            sum_pos += node->i * node->depth;
            sI += node->i;
        }

        if (sI > 0) {
            double meanI = sum_pos / sI;
            double Tmax = params_.M * params_.dt;
            double delta_T = Tmax - params_.T0;

            if (delta_T > 0.0) speed_ = (meanI - meanI0_) / delta_T;
            else speed_ = 0.0;

            double varI = 0.0;
            for (Node* node : tree_.GetAllNodes()) {
                varI += node->i * (node->depth - meanI) * (node->depth - meanI);
            }
            sdI_ = std::sqrt(varI / sI);
            final_finf_ = sI;
        }

        std::cout << "\nSimulation Complete. Speed (c): " << speed_ << " | sdI: " << sdI_ << "\n";
    }

    void SaveData() {
        const std::string output_dir = CONST::OUTPUT + "/" + NAME;
        const std::string data_dir = output_dir + "/" + CONST::DATA;
        std::cout << "Saving timeseries data to: " << data_dir << "\n";

        std::ofstream out(data_dir + "/timeseries.csv");
        out << "step,total_i\n";
        for (size_t i = 0; i < total_i_history_.size(); ++i) {
            out << (i + 1) << "," << total_i_history_[i] << "\n";
        }
        out.close();

        std::ofstream w_out(data_dir + "/wave.csv");
        w_out << "step,depth,i\n";
        for (size_t st = 0; st < wave_history_.size(); ++st) {
            for (const auto& kv : wave_history_[st]) {
                w_out << (st + 1) << "," << kv.first << "," << kv.second << "\n";
            }
        }
        w_out.close();

        std::ofstream f_obs(data_dir + "/observables.txt");
        f_obs << "R0=" << params_.R0 << "\n";
        f_obs << "Ub=" << params_.Ub << "\n";
        f_obs << "a=" << params_.a << "\n";
        f_obs << "N=" << params_.N << "\n";
        f_obs << "dt=" << params_.dt << "\n";
        f_obs << "M=" << params_.M << "\n";
        f_obs << "seed=" << params_.seed << "\n";

        f_obs << "speed=" << speed_ << "\n";
        f_obs << "sdI=" << sdI_ << "\n";
        f_obs << "finf=" << final_finf_ << "\n";
        f_obs.close();
    }
private:
    ModelParams params_;
    TreeTopology tree_;
    std::vector<double> total_i_history_;
    std::vector<std::map<int, double>> wave_history_;
    std::mt19937 rng_;
    std::vector<double> K_table_;

    double meanI0_ = 0.0;
    double speed_ = 0.0;
    double sdI_ = 0.0;
    double final_finf_ = 0.0;

    double GetStochasticCorrection(const double& val) {
        const int64_t& N = params_.N;
        const double lambda = val * N;

        if (lambda >= 10.0) return val;

        const int xm = static_cast<int>(std::round(6.0 * lambda));
        if (xm == 0) return 0.0;

        std::uniform_real_distribution<double> dist(0.0, 1.0);
        int count = 0;
        for (int n = 0; n < xm; ++n) {
            if (xm * dist(rng_) < lambda) ++count;
        }

        return static_cast<double>(count) / N;
    }

    void Step(int step) {
        const double& dt = params_.dt;
        const double& R0 = params_.R0;
        const auto& nodes = tree_.GetAllNodes();

        int max_depth = 0;
        for (Node* n : nodes) {
            if (n->depth > max_depth) max_depth = n->depth;
        }

        std::vector<double> sum_I(max_depth + 1, 0.0);
        std::vector<double> sum_R(max_depth + 1, 0.0);

        for (Node* n : nodes) {
            sum_I[n->depth] += n->i;
            sum_R[n->depth] += n->r;
        }

        std::map<int, double> current_wave;
        for (int d = 0; d <= max_depth; ++d) {
            if (sum_I[d] > 1e-9) current_wave[d] = sum_I[d];
        }
        wave_history_.push_back(current_wave);

        std::vector<double> Q_depth(max_depth + 1, 0.0);
        std::vector<double> P_depth(max_depth + 1, 0.0);

        for (int d_j = 0; d_j <= max_depth; ++d_j) {
            for (int d_k = 0; d_k <= max_depth; ++d_k) {
                int dist = std::abs(d_j - d_k);
                double k_val = K_table_[dist];
                Q_depth[d_j] += k_val * sum_R[d_k];
                P_depth[d_j] += k_val * sum_I[d_k];
            }
        }

        double current_total_i = 0.0;

        for (Node* n : nodes) {
            double I_old = n->i;
            double R_old = n->r;

            double Q_val = Q_depth[n->depth];
            double P_val = P_depth[n->depth];

            double I_new = I_old * (1.0 + dt * (R0 * Q_val - 1.0));
            double R_new = R_old * (1.0 - dt * R0 * P_val) + dt * I_old;

            n->r = std::max(0.0, R_new);

            n->i = GetStochasticCorrection(std::max(0.0, I_new));
            n->is_infected = (n->i > CONST::EPS);

            current_total_i += n->i;
        }

        total_i_history_.push_back(current_total_i);

        size_t current_nodes_count = nodes.size();
        for (size_t j = 0; j < current_nodes_count; ++j) {
            Node* curr = nodes[j];
            if (!curr->is_infected) continue;
            if (curr->i * params_.N < 1.0) continue;

            double expected_mutants_frac = dt * params_.Ub * curr->i;
            if (expected_mutants_frac > 0.0) Mutation(expected_mutants_frac, curr, step);
        }
    }

    // 4. ПРИНИМАЕМ step СЮДА
    void Mutation(double expected_mutants_frac, Node* curr, int step) {
        double actual_mutants_frac = GetStochasticCorrection(expected_mutants_frac);

        if (actual_mutants_frac > 0.0) {
            int total_mutants = static_cast<int>(std::round(actual_mutants_frac * params_.N));

            if (total_mutants > 0) {
                int mutants_l = std::binomial_distribution<int>(total_mutants, 0.5)(rng_);
                int mutants_r = total_mutants - mutants_l;

                double frac_l = static_cast<double>(mutants_l) / params_.N;
                double frac_r = static_cast<double>(mutants_r) / params_.N;

                double total_frac_to_move = frac_l + frac_r;
                if (total_frac_to_move > curr->i) {
                    double scale = curr->i / total_frac_to_move;
                    frac_l *= scale;
                    frac_r *= scale;
                }

                curr->i -= (frac_l + frac_r);

                if (frac_l > 0.0) {
                    if (curr->l_child) {
                        curr->l_child->i += frac_l;
                        curr->l_child->is_infected = true;
                    }
                    else tree_.AddChild(curr, step, frac_l, 0.0);
                }

                if (frac_r > 0.0) {
                    if (curr->r_child) {
                        curr->r_child->i += frac_r;
                        curr->r_child->is_infected = true;
                    }
                    else tree_.AddChild(curr, step, frac_r, 0.0);
                }
            }
        }
    }

    void SaveSnapshot(int step) {
        const std::string output_dir = CONST::OUTPUT + "/" + NAME;
        const std::string data_dir = output_dir + "/" + CONST::DATA;
        std::ofstream out(data_dir + "/snapshot_" + std::to_string(step) + ".csv");

        out << "id,parent_id,depth,creation_step,i,r,is_infected\n";

        for (Node* node : tree_.GetAllNodes()) {
            int parent_id = (node->parent != nullptr) ? node->parent->id : -1;
            out << node->id << ","
                << parent_id << ","
                << node->depth << ","
                << node->creation_step << ","
                << node->i << ","
                << node->r << ","
                << (node->is_infected ? 1 : 0) << "\n";
        }
        out.close();
    }
};

// ================== Execution and Management ==================
void InitDirectory(const std::string& name) {
    const std::string output_dir = CONST::OUTPUT + "/" + name;
    const std::string data_dir = output_dir + "/" + CONST::DATA;
    fs::create_directories(data_dir);

    for (const auto& entry : fs::directory_iterator(data_dir)) {
        fs::remove_all(entry.path());
    }
}

void TEST(const ModelParams& params, const std::string& name) {
    NAME = name;
    InitDirectory(NAME);
    TreeEvolutionSimulator sim(params);
    sim.Run();

    const std::string output_dir = CONST::OUTPUT + "/" + NAME;
    const std::string data_dir = output_dir + "/" + CONST::DATA;
    sim.SaveData();

    std::string cmd = "python tree_visualize.py -d " + data_dir + " -o " + output_dir;
    std::cout << "Executing: " << cmd << "\n";
    system(cmd.c_str());
}

int main() {
    ModelParams params;
    params.M = 2'500;
    params.T0 = params.M / 10;
    params.iter = 100;
    params.R0 = 1.8;
    params.Ub = 1e-3;
    params.a = 14.0;
    params.seed = 1;

    TEST(params, "exp1");

    return 0;
}