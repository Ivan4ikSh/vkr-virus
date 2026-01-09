#include "domain.h"
#include "monte_carlo.h"

int main() {
    cout << setprecision(3);
    cout << "Monte Carlo Population Genetics Simulation\n";

    {
        // Основной эксперимент с постоянным распределением отбора
        cout << "=========================================\n";
        StartParams params;
        params.exp_name = "test1";
        params.N = 1000;
        params.L = 100;
        params.s0 = 0.05;
        params.r = 0;
        params.f0 = 0;
        params.muL = 0.1;
        params.tf = 150;
        params.M = 3;
        // Запускаем симуляцию
        MonteCarlo MC(params);
        MC.RunSimulation();
        cout << "=========================================\n";
        system(("python fig1.py --dir " + params.output_dir + "/" + params.exp_name + " --output " + params.plot_dir + "/" + params.exp_name).c_str());
        cout << "=========================================\n";
    }
    {
        // Основной эксперимент с постоянным распределением отбора
        cout << "=========================================\n";
        StartParams params;
        params.exp_name = "test2";
        params.N = 1000;
        params.L = 100;
        params.s0 = 0.1;
        params.r = 0;
        params.f0 = 0;
        params.muL = 0.1;
        params.tf = 150;
        params.M = 3;
        // Запускаем симуляцию
        MonteCarlo MC(params);
        MC.RunSimulation();
        cout << "=========================================\n";
        system(("python fig1.py --dir " + params.output_dir + "/" + params.exp_name + " --output " + params.plot_dir + "/" + params.exp_name).c_str());
        cout << "=========================================\n";
    }
    {
        // Основной эксперимент с постоянным распределением отбора
        cout << "=========================================\n";
        StartParams params;
        params.exp_name = "test3";
        params.N = 1000;
        params.L = 100;
        params.s0 = 0.2;
        params.r = 0;
        params.f0 = 0;
        params.muL = 0.1;
        params.tf = 150;
        params.M = 3;
        // Запускаем симуляцию
        MonteCarlo MC(params);
        MC.RunSimulation();
        cout << "=========================================\n";
        system(("python fig1.py --dir " + params.output_dir + "/" + params.exp_name + " --output " + params.plot_dir + "/" + params.exp_name).c_str());
        cout << "=========================================\n";
    }
    return 0;
}