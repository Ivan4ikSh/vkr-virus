#pragma once
#include "domain.h"

namespace FILEUTILS {
    void CreateDirectory(const string& path) { filesystem::create_directories(path); }
    bool FileExists(const string& filename) { return filesystem::exists(filename); }

    bool WriteCSV(const string& filename, const vector<vector<double>>& data, const vector<string>& headers) {
        ofstream file(filename);
        if (!file.is_open()) return false;

        // Записываем заголовки
        if (!headers.empty()) {
            for (size_t i = 0; i < headers.size(); ++i) {
                file << headers[i];
                if (i < headers.size() - 1) file << ",";
            }
            file << "\n";
        }

        // Записываем данные
        for (const auto& row : data) {
            for (size_t i = 0; i < row.size(); ++i) {
                file << scientific << setprecision(6) << row[i];
                if (i < row.size() - 1) file << ",";
            }
            file << "\n";
        }

        file.close();
        return true;
    }

    bool WriteVector(const string& filename, const vector<double>& data) {
        ofstream file(filename);
        if (!file.is_open()) return false;

        file << scientific << setprecision(6);
        for (const auto& val : data) {
            file << val << "\n";
        }

        file.close();
        return true;
    }

    bool WriteMatrix(const string& filename, const vector<vector<double>>& matrix) {
        ofstream file(filename);
        if (!file.is_open()) return false;

        file << scientific << setprecision(6);
        for (const auto& row : matrix) {
            for (size_t i = 0; i < row.size(); ++i) {
                file << row[i];
                if (i < row.size() - 1) file << ",";
            }
            file << "\n";
        }

        file.close();
        return true;
    }

    bool WriteHistogramData(const string& filename, const vector<int>& times, const vector<vector<double>>& xx, const vector<vector<double>>& nn) {
        ofstream file(filename);
        if (!file.is_open()) return false;

        file << "time,bin_center,count\n";
        file << scientific << setprecision(6);

        for (size_t i = 0; i < times.size(); ++i) {
            for (size_t j = 0; j < xx[i].size(); ++j) {
                file << times[i] << "," << xx[i][j] << "," << nn[i][j] << "\n";
            }
        }

        file.close();
        return true;
    }
}
