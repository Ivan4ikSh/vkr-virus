#include "domain.h"
#include "monte_carlo.h"

int main() {
    cout << setprecision(3);
    cout << "Monte Carlo Population Genetics Simulation\n";

    {
        // Основной эксперимент с постоянным распределением отбора
        cout << "=========================================\n";
        // Флаг на использование Модели Лукши (по умолчанию false)
        CONST::USE_LUKSZA = true;
        // Параметры модели Монте-Карло
        StartParams params;
        params.exp_name = "test";
        params.N = 1000;      // Размер популяции (можно увеличить до 10^4 для лучшей статистики)
        params.L = 100;       // Более реалистичная длина для HA белка (бинарное представление)
        params.s0 = 0.05;      // Не используется в модели Луки
        params.r = 0.0;       // Нет рекомбинации у гриппа
        params.f0 = 0.01;
        params.muL = 0.5;     // В среднем 1 мутация на гаплотип за поколение (реалистичнее)
        params.tf = 300;      // 500 поколений для наблюдения долговременной динамики
        params.M = 3;         // Начальное разнообразие

        // Параметры модели Лукши-Лассига
        LukszaParams l_params;
        l_params.D0 = 14.0;           // Соответствует статье
        l_params.L_ep = 50;           // длина эпитопной части (33%)
        l_params.T_memory = 7;        // иммунной памяти
        l_params.sigma_ne = 0.5;      // Оптимально по статье

        // Запускаем симуляцию
        MonteCarlo MC(params, l_params);
        auto start = chrono::high_resolution_clock::now();
        MC.RunSimulation();
        auto end = chrono::high_resolution_clock::now();
        cout << "Simulation time: " << chrono::duration_cast<chrono::seconds>(end - start).count() << "s\n\n";
        cout << "=========================================\n";
        system(("python fig1.py --dir " + params.output_dir + "/" + params.exp_name + " --output " + params.plot_dir + "/" + params.exp_name).c_str());
        cout << "=========================================\n";
    }
    
    return 0;
}