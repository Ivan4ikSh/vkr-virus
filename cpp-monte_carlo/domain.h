#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

namespace RANDOM {
    uint32_t Hash(uint32_t state) {
        state ^= 2747636219u;
        state *= 2654435769u;
        state ^= state >> 16;
        state *= 2654435769u;
        state ^= state >> 16;
        return state * 2654435769u;
    }

    float ScaleToRange01(uint32_t value) {
        return static_cast<float>(Hash(value)) / 4294967295.0f;
    }
}

namespace CONST {
    const double PI = 3.1415926535897932;
    const bool YES_GENEALOGY = true;
}

enum class Distribution {
    CONSTANT,      // Все мутации имеют одинаковый эффект
    EXPONENTIAL,   // Экспоненциальное распределение
    HALF_GAUSSIAN, // Полу-гауссово распределение
    UNIFORM        // Равномерное распределение
};

struct Args {
    Distribution dist_s;            // distribution_s - одно из распределений
    double recomb_cnt;              // r - частота рекомбинаций на геном
    int cross_cnt;                  // M - число кроссинговеров
    double adaptive_landscape;      // s0 - средний коэффициент отбора
    int locus_cnt;                  // L - число локусов
    int population_cnt;             // N - численность популяции
    int full_time_in_epoch;         // tf - полное время в поколениях
    double start_good_allele_freq;  // f0 - начальная частота благоприятного аллеля
    double mutation_freq;           // muL - частота мутаций на геном

    string exp_name;
};

namespace DATA {
    struct Population {
        vector<double> s;  // коэффициенты отбора
        vector<vector<int>> A;  // матрица предков
        vector<vector<int>> K;  // бинарные последовательности ДНК

        // Временные данные
        vector<vector<double>> W;  // приспособленность
        vector<vector<double>> P1, PL;  // метки родителей
        vector<vector<double>> fsite;  // частота аллелей по локусам

        // Статистики
        vector<double> kav, Vark, fsurvive, C, Call, meanW;

        // Новые поколения
        vector<vector<int>> Knew, Anew, Pnew;

        // Параметры
        int tint = 0;  // интервал времени для графиков
        double fsample = 0.0;  // процент выборки для пар
        double mu = 0;  // частота мутаций на локус
    };

    struct Output {
        vector<double> time_points;
        vector<double> kav;
        vector<double> Vark;
        vector<double> fsurvive;
        vector<double> C;
        vector<double> Call;
        vector<double> meanW;
        vector<vector<double>> fitness_histograms;
        vector<vector<double>> fsite;
        double theoretical_velocity;

        vector<int> histogram_times; // времена для гистограмм
        vector<string> histogram_colors; // цвета для гистограмм
    };
}

namespace SAVE {
    void Statistics(const DATA::Output& output, const Args& params, const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Ошибка открытия файла: " << filename << std::endl;
            return;
        }

        // Простой заголовок без специальных символов
        file << "Time,kav,Vark,fsurvive,C,Call,meanW\n";
        for (size_t i = 0; i < output.time_points.size(); ++i) {
            file << output.time_points[i] << ","
                << output.kav[i] << ","
                << output.Vark[i] << ","
                << output.fsurvive[i] << ","
                << output.C[i] << ","
                << output.Call[i] << ","
                << output.meanW[i] << "\n";
        }
        file.close();
    }

    void FitnessHistograms(const DATA::Output& output, const string& filename) {
        ofstream file(filename);
        if (file.is_open()) {
            // Сохраняем данные для построения бегущей волны
            for (size_t i = 0; i < output.fitness_histograms.size(); ++i) {
                file << output.histogram_times[i] << "," << output.histogram_colors[i];
                const auto& hist = output.fitness_histograms[i];
                for (size_t j = 0; j < hist.size(); j += 2) {
                    file << "," << hist[j] << "," << hist[j + 1]; // bin_center, frequency
                }
                file << "\n";
            }
        }
    }

    void AlleleFrequencies(const DATA::Output& output, const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Ошибка открытия файла: " << filename << std::endl;
            return;
        }

        // Заголовок
        file << "Time";
        for (size_t locus = 0; locus < output.fsite[0].size(); ++locus) {
            file << ",Locus" << locus + 1;
        }
        file << "\n";

        // Данные
        for (size_t t = 0; t < output.fsite.size(); ++t) {
            file << output.time_points[t];
            for (size_t locus = 0; locus < output.fsite[t].size(); ++locus) {
                file << "," << output.fsite[t][locus];
            }
            file << "\n";
        }
        file.close();
    }

    void Parameters(const Args& params, double theoretical_velocity, const string& filename) {
        ofstream file(filename);
        if (file.is_open()) {
            file << "theoretical_velocity=" << theoretical_velocity << endl;
            file << "N=" << params.population_cnt << endl;
            file << "L=" << params.locus_cnt << endl;
            file << "s0=" << params.adaptive_landscape << endl;
            file << "r=" << params.recomb_cnt << endl;
            file << "f0=" << params.start_good_allele_freq << endl;
            file << "muL=" << params.mutation_freq << endl;
            file << "tf=" << params.full_time_in_epoch << endl;
            file << "M=" << params.cross_cnt << endl;
            file << "run=1" << endl;
            file << "distribution_s=";
            switch (params.dist_s) {
            case Distribution::CONSTANT: file << "const" << endl; break;
            case Distribution::EXPONENTIAL: file << "exponential" << endl; break;
            case Distribution::HALF_GAUSSIAN: file << "halfgaussian" << endl; break;
            default: file << "constant" << endl; break;
            }
        }
    }
}

// Функция для получения текущего времени в формате [HH:MM:SS]
string GetCurrentTime() {
    auto now = chrono::system_clock::now();
    auto time_t = chrono::system_clock::to_time_t(now);
    auto tm = *localtime(&time_t);

    stringstream ss;
    ss << "[" << put_time(&tm, "%H:%M:%S") << "]";
    return ss.str();
}

// Функция логгирования
void Log(const string& stage, const string& function_name, const string& process_name) {
    cout << GetCurrentTime() << " : " << stage << " " << function_name << " " << process_name << endl;
}
