#pragma once

#include <algorithm>
#include <chrono>
#include <execution>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using namespace std;

namespace CONST {
    const double PI = 3.1415926535897932;
}

// Типы распределений для коэффициентов отбора
enum DistributionType {
    EXPONENTIAL,    // Экспоненциальное: s = -s0 * log(rand)
    CONSTANT,       // Постоянное: s = s0
    HALF_GAUSSIAN   // Полугауссово: s = s0 * sqrt(pi/2) * abs(randn)
};

// Структура для хранения всех параметров симуляции
struct StartParams {
    // Основные параметры
    int N;                    // Размер популяции
    int L;                    // Количество локусов
    double s0;                // Средний коэффициент отбора
    double r;                 // Частота рекомбинаций
    double f0;                // Начальная частота благоприятного аллеля
    double muL;               // Частота мутаций на геном
    int tf;                   // Полное время в поколениях
    int M;                    // Число кроссинговеров

    // Распределение коэффициентов отбора
    DistributionType distribution;

    // Флаги и дополнительные параметры
    bool track_genealogy;// Отслеживать генеалогию
    string exp_name;     // Имя эксперимента
    string output_dir;   // Директория для сохранения файлов
    string plot_dir;     // Директория для сохранения графиков

    // Конструктор с параметрами по умолчанию
    StartParams() :
        N(500), L(200), s0(0.1), r(0), f0(0), muL(0.01), tf(100), M(3),
        distribution(DistributionType::CONSTANT),
        track_genealogy(false),
        output_dir("results"), plot_dir("graphics"),
        exp_name("default_experiment") {
    }
};

struct SimulationParams {
    // Вычисленные параметры
    double mu;      // частота мутаций на локус
    vector<int> T;  // времена (0, 1, 2, ..., tf)
    // Массив коэффициентов отбора для каждого локуса
    vector<double> s;
    // Аналитическая и численная скорости адаптации
    double V_an;
    double V_num;
    // Настройки для графиков/анализа
    int t_int;           // интервал времени для графиков
    double f_sample;     // процент выборки для пар
    // Состояние популяции
    vector<vector<int>> K;     // бинарные последовательности ДНК (N x L)
    vector<vector<int>> A;     // матрица предков (N x L)
    vector<vector<int>> P;     // матрица родителей (N x L)
    // Для временного хранения новых состояний
    vector<vector<int>> K_new;
    vector<vector<int>> A_new;
    vector<vector<int>> P_new;
    // Наблюдаемые величины во времени
    vector<vector<double>> W;       // приспособленности (N x (tf+1))
    vector<vector<int>> P1;         // родители для 1-го локуса (N x (tf+1))
    vector<vector<int>> PL;         // родители для последнего локуса (N x (tf+1))
    vector<vector<double>> f_site;  // частоты аллелей по локусам ((tf+1) x L)
    vector<double> k_av;            // среднее число аллелей на геном
    vector<double> V_ark;           // дисперсия числа аллелей
    vector<double> f_survive;       // доля выживших
    vector<double> C;               // доля пар с общим предком
    vector<double> C_all;           // ?
    vector<double> mean_W;          // средняя приспособленность

    string colors;
    vector<double> dist_over_L;             // добавьте это для хранения dist/L
    vector<double> f1_site;                 // теоретические частоты
    vector<vector<double>> fitness_hist_xx; // центры бинов для гистограмм
    vector<vector<double>> fitness_hist_nn; // высоты бинов
    vector<int> hist_times;                 // моменты времени для гистограмм
    vector<vector<double>> freq_hist_xx;    // центры бинов для частот
    vector<vector<double>> freq_hist_nn;    // высоты бинов для частот
    vector<int> freq_hist_times;            // моменты времени для гистограмм частот
};