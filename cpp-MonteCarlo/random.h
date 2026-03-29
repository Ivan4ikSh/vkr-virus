#pragma once

namespace RANDOM {
    uint32_t Hash() {
        srand(time(NULL));
        uint32_t state = rand();

        state ^= 2747636219u;
        state *= 2654435769u;
        state ^= state >> 16;
        state *= 2654435769u;
        state ^= state >> 16;
        return state * 2654435769u;
    }

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

// Класс для генерации случайных чисел
class RandomGenerator {
public:
    RandomGenerator(uint32_t seed = 0) : generator(seed), uniform_dist(0.0, 1.0) {}

    double Random() { return uniform_dist(generator); }
    void SetSeed(int seed) { generator.seed(seed); }
    // Генерация экспоненциального распределения
    double Exponential(double lambda) {
        exponential_distribution<double> dist(lambda);
        return dist(generator);
    }
    // Генерация нормального распределения
    double Normal(double mean, double stddev) {
        normal_distribution<double> dist(mean, stddev);
        return dist(generator);
    }
    // Генерация случайного целого числа в диапазоне [min, max]
    int RandomInt(int min, int max) {
        uniform_int_distribution<int> dist(min, max);
        return dist(generator);
    }
private:
    mt19937 generator;
    uniform_real_distribution<double> uniform_dist;
};
