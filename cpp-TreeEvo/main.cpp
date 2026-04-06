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

// ================== Бинарная топология с антигенными координатами ==================
struct Node {
    int id;
    int depth;
    int creation_step;
    double i;
    double r;
    bool is_infected;

    double edge_length; // Длина прыжка мутации 
    double absolute_x;  // Абсолютная координата в антигенном пространстве

    Node* parent;
    Node* l_child;
    Node* r_child;

    Node(int node_id, int node_depth, int step_born, double initial_i, double initial_r, double edge_len, Node* parent_node = nullptr)
        : id(node_id), depth(node_depth), creation_step(step_born), i(initial_i), r(initial_r),
        is_infected(initial_i * 1e8 >= 1.0), edge_length(edge_len), parent(parent_node), l_child(nullptr), r_child(nullptr) {

        if (parent_node) absolute_x = parent_node->absolute_x + edge_length;
        else absolute_x = 0.0;
    }
};

class TreeTopology {
public:
    TreeTopology() {
        root_ = new Node(0, 0, 0, 0.0, 1.0 - 1e-4, 0.0);
        all_nodes_.push_back(root_);
        node_count_ = 1;

        Node* curr = root_;
        for (int i = 1; i <= 12; ++i) {
            Node* child = new Node(node_count_++, i, 0, 0.0, 0.0, 1.0, curr);
            all_nodes_.push_back(child);
            curr->l_child = child;
            curr = child;
        }
        curr->i = 1e-4;
        curr->is_infected = true;
    }

    ~TreeTopology() { for (Node* node : all_nodes_) delete node; }

    Node* AddChild(Node* parent, int step_born, double initial_i, double initial_r, double edge_len) {
        Node* child = new Node(node_count_++, parent->depth + 1, step_born, initial_i, initial_r, edge_len, parent);
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
    TreeEvolutionSimulator(const ModelParams& params) : params_(params), rng_(params.seed) {}

    void Run() {
        std::cout << "Starting binary simulation...\n";
        int T0_step = static_cast<int>(std::round(params_.T0 / params_.dt));

        for (int step = 1; step <= params_.M; ++step) {
            Step(step);
            if (step == T0_step) {
                double sPos = 0.0;
                double sI = 0.0;
                for (Node* n : tree_.GetAllNodes()) { sPos += n->i * n->absolute_x; sI += n->i; }
                if (sI > 0) meanI0_ = sPos / sI;
            }
            if (step % params_.iter == 0) std::cout << "Step: " << step << " | Nodes: " << tree_.GetNodeCount() << "\n";
            if (step == params_.M / 4 || step == params_.M / 2 || step == 3 * params_.M / 4 || step == params_.M) SaveSnapshot(step);
        }

        double sPos = 0.0;
        double sI = 0.0;
        for (Node* n : tree_.GetAllNodes()) { sPos += n->i * n->absolute_x; sI += n->i; }
        if (sI > 0) {
            double meanI = sPos / sI;
            double delta_T = (params_.M * params_.dt) - params_.T0;
            speed_ = (delta_T > 0) ? (meanI - meanI0_) / delta_T : 0.0;
            final_finf_ = sI;
        }
        std::cout << "\nComplete. Speed: " << speed_ << " | Nodes: " << tree_.GetNodeCount() << "\n";
    }

    void SaveData() {
        const std::string output_path = CONST::OUTPUT + "/" + NAME;
        const std::string data_dir = output_path + "/" + CONST::DATA;

        // 1. Сохранение временного ряда общей инфекции
        std::ofstream out(data_dir + "/timeseries.csv");
        out << "step,total_i\n";
        for (size_t i = 0; i < total_i_history_.size(); ++i) {
            out << (i + 1) << "," << total_i_history_[i] << "\n";
        }
        out.close();

        // 2. Сохранение данных бегущей волны (wave.csv)
        std::ofstream w_out(data_dir + "/wave.csv");
        w_out << "step,depth,i\n";
        for (size_t st = 0; st < wave_history_.size(); ++st) {
            for (const auto& kv : wave_history_[st]) {
                w_out << (st + 1) << "," << kv.first << "," << kv.second << "\n";
            }
        }
        w_out.close();

        // 3. Сохранение ВСЕХ параметров в observables.txt
        std::ofstream f_obs(data_dir + "/observables.txt");
        // Входные параметры модели
        f_obs << "R0=" << params_.R0 << "\n";
        f_obs << "Ub=" << params_.Ub << "\n";
        f_obs << "a=" << params_.a << "\n";
        f_obs << "N=" << params_.N << "\n";
        f_obs << "seed=" << params_.seed << "\n";

        // Результирующие показатели
        f_obs << "speed=" << speed_ << "\n";
        f_obs << "finf=" << final_finf_ << "\n";

        f_obs.close();
        std::cout << "Data saved to " << data_dir << std::endl;
    }

private:
    ModelParams params_;
    TreeTopology tree_;
    std::vector<double> total_i_history_;
    std::vector<std::map<int, double>> wave_history_;
    std::mt19937 rng_;
    double meanI0_ = 0.0;
    double speed_ = 0.0;
    double final_finf_ = 0.0;

    double GetDistance(Node* u, Node* v) {
        double d = 0.0;
        Node* currU = u, * currV = v;
        while (currU->depth > currV->depth) { d += currU->edge_length; currU = currU->parent; }
        while (currV->depth > currU->depth) { d += currV->edge_length; currV = currV->parent; }
        while (currU != currV) { d += currU->edge_length + currV->edge_length; currU = currU->parent; currV = currV->parent; }
        return d;
    }

    double GetStochasticCorrection(const double& val) {
        double lambda = val * params_.N;
        if (lambda >= 10.0) return val;

        int xm = static_cast<int>(std::round(6.0 * lambda));
        if (xm == 0) return 0.0;

        std::uniform_real_distribution<double> dist(0.0, 1.0);
        int count = 0;
        for (int n = 0; n < xm; ++n) if (xm * dist(rng_) < lambda) ++count;
        return static_cast<double>(count) / params_.N;
    }

    void Step(int step) {
        const auto& nodes = tree_.GetAllNodes();
        int num_nodes = nodes.size();

        std::map<int, double> current_wave;
        for (Node* n : nodes) {
            if (n->i > CONST::EPS) current_wave[n->depth] += n->i;
        }
        wave_history_.push_back(current_wave);

        std::vector<double> Q(num_nodes, 0.0);
        std::vector<double> P(num_nodes, 0.0);

        for (int j = 0; j < num_nodes; ++j) {
            for (int k = 0; k < num_nodes; ++k) {
                if (j == k) continue;
                double dist = GetDistance(nodes[j], nodes[k]);
                double k_val = 1.0 - std::exp(-dist / params_.a);
                if (nodes[j]->absolute_x > nodes[k]->absolute_x) Q[j] += k_val * nodes[k]->r;
                if (nodes[k]->absolute_x > nodes[j]->absolute_x) P[j] += k_val * nodes[k]->i;
            }
        }

        double total_i = 0.0;
        for (int j = 0; j < num_nodes; ++j) {
            double I_new = nodes[j]->i * (1.0 + params_.dt * (params_.R0 * Q[j] - 1.0));
            nodes[j]->r = std::max(0.0, nodes[j]->r * (1.0 - params_.dt * params_.R0 * P[j]) + params_.dt * nodes[j]->i);
            nodes[j]->i = GetStochasticCorrection(std::max(0.0, I_new));
            if (nodes[j]->i * params_.N < 1.0) nodes[j]->i = 0.0;
            nodes[j]->is_infected = (nodes[j]->i > CONST::EPS);
            total_i += nodes[j]->i;
        }
        total_i_history_.push_back(total_i);

        // Мутации в бинарном дереве
        for (int j = 0; j < num_nodes; ++j) {
            if (!nodes[j]->is_infected) continue;
            double expected = params_.dt * params_.Ub * nodes[j]->i * params_.N;
            int actual = std::binomial_distribution<int>(static_cast<int>(std::round(expected * 2.0)), 0.5)(rng_);
            if (actual > 0) Mutation(actual, nodes[j], step);
        }
    }

    void Mutation(int total_mutants, Node* curr, int step) {
        // Делим мутантов между левой и правой веткой
        int m_l = std::binomial_distribution<int>(total_mutants, 0.5)(rng_);
        int m_r = total_mutants - m_l;

        double f_l = (double)m_l / params_.N, f_r = (double)m_r / params_.N;
        if (f_l + f_r > curr->i) { double s = curr->i / (f_l + f_r); f_l *= s; f_r *= s; }
        curr->i -= (f_l + f_r);

        std::uniform_real_distribution<double> lens(0.5, 1.5);
        if (f_l > 0.0) {
            if (curr->l_child) { curr->l_child->i += f_l; curr->l_child->is_infected = true; }
            else tree_.AddChild(curr, step, f_l, 0.0, lens(rng_));
        }
        if (f_r > 0.0) {
            if (curr->r_child) { curr->r_child->i += f_r; curr->r_child->is_infected = true; }
            else tree_.AddChild(curr, step, f_r, 0.0, lens(rng_));
        }
    }

    void SaveSnapshot(int step) {
        const std::string data_dir = CONST::OUTPUT + "/" + NAME + "/" + CONST::DATA;
        std::ofstream out(data_dir + "/snapshot_" + std::to_string(step) + ".csv");
        out << "id,parent_id,depth,creation_step,i,r,is_infected,edge_len,absolute_x\n";
        for (Node* n : tree_.GetAllNodes()) {
            int p = n->parent ? n->parent->id : -1;
            out << n->id << "," << p << "," << n->depth << "," << n->creation_step << "," << n->i << "," << n->r << "," << n->is_infected << "," << n->edge_length << "," << n->absolute_x << "\n";
        }
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
    params.M = 5'000;
    params.T0 = params.M / 10;
    params.iter = params.M / 10;
    params.R0 = 1.8;
    params.Ub = 1e-4;
    params.a = 14.0;
    params.seed = 1;

    TEST(params, "exp1");

    return 0;
}