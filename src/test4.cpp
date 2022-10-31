#include "tester.hpp"

size_t Tester::moc_found = 0;
std::mt19937 gen;

inline bool find_nlogm(std::vector<std::vector<int64>>& matrix, int64 k) {
    size_t n = matrix.size();
    size_t m = matrix[0].size();
    if (n <= m) {
        for (size_t i = 0; i < n; i++) {
            int64 l = 0;
            int64 r = m - 1;
            while (l <= r) {
                int64 mid = (l + r) / 2;
                if (matrix[i][mid] == k) {
                    return true;
                }
                else if (matrix[i][mid] < k) {
                    l = mid + 1;
                }
                else {
                    r = mid - 1;
                }
            }
        }
        return false;
    }
    else {
        for (size_t i = 0; i < m; i++) {
            int64 l = 0;
            int64 r = n - 1;
            while (l <= r) {
                int64 mid = (l + r) / 2;
                if (matrix[mid][i] == k) {
                    return true;
                }
                else if (matrix[mid][i] < k) {
                    l = mid + 1;
                }
                else {
                    r = mid - 1;
                }
            }
        }
        return false;
    }
}

inline bool find_improved_nlogm(std::vector<std::vector<int64>>& matrix, int64 k) {
    size_t n = matrix.size();
    size_t m = matrix[0].size();
    if (n <= m) {
        size_t lower_bound = m - 1;
        for (size_t i = 0; i < n;) {
            int64 l = 0;
            int64 r = lower_bound;
            while (l <= r) {
                int64 mid = (l + r) / 2;
                if (matrix[i][mid] >= k) {
                    r = mid - 1;
                    lower_bound = mid;
                }
                else {
                    l = mid + 1;
                }
            }
            if (matrix[i][lower_bound] == k) {
                return true;
            }
            else {
                i++;
            }
        }
        return false;
    }
    else {
        size_t lower_bound = n - 1;
        for (size_t j = 0; j < m;) {
            int64 l = 0;
            int64 r = lower_bound;
            while (l <= r) {
                int64 mid = (l + r) / 2;
                if (matrix[mid][j] >= k) {
                    r = mid - 1;
                    lower_bound = mid;
                }
                else {
                    l = mid + 1;
                }
            }
            if (matrix[lower_bound][j] == k) {
                return true;
            }
            else {
                j++;
            }
        }
        return false;
    }
}

inline bool find_n_plus_m(std::vector<std::vector<int64>>& matrix, int64 k) {
    int64 i = 0;
    int64 j = matrix[0].size() - 1;
    while (i < matrix.size() && j >= 0) {
        if (matrix[i][j] == k) {
            return true;
        }
        else if (matrix[i][j] > k) {
            j--;
        }
        else {
            i++;
        }
    }
    return false;
}

inline bool find_improved_n_plus_m(std::vector<std::vector<int64>>& matrix, int64 k) {
    int64 i = 0;
    int64 j = matrix[0].size() - 1;
    int64 mult = 1;
    while (i < matrix.size() && j >= 0) {
        if (matrix[i][j] > k) {
            j -= mult;
            if (j >= 0 && matrix[i][j] < k) {
                j = std::lower_bound(matrix[i].begin() + j, matrix[i].begin() + j + mult, k) - matrix[i].begin();
                i++;
                mult = 1;
            }
            else {
                mult *= 2;
            }
        }
        else if (matrix[i][j] < k) {
            i++;
        }
        else {
            return true;
        }

    }
    return false;
}

int main() {

    gen.seed(static_cast<unsigned>(std::chrono::high_resolution_clock::now().time_since_epoch().count()));

    std::ofstream fout;
    fout.open("task2.xls");
    fout << std::fixed << std::setprecision(9);

    std::vector<std::vector<int64>> matrix;
    size_t m = 1ll << 13;

    fout << "M size\ttime (seconds)\n";

    fout << "nlogm\n";
    for (size_t n = 1; n <= m; n <<= 1) {
        fout << n << "\t";
        Tester::generate2(matrix, n, m);
        int64 target = matrix[n - 1][0] + 1;
        fout << Tester::test(100, find_nlogm, matrix, target) << "\n";
    }
    fout << "\n\n";

    fout << "n + m\n";
    for (size_t n = 1; n <= m; n <<= 1) {
        fout << n << "\t";
        Tester::generate2(matrix, n, m);
        int64 target = matrix[n - 1][0] + 1;
        fout << Tester::test(100, find_n_plus_m, matrix, target) << "\n";
    }
    fout << "\n\n";

    fout << "n + m (imp)\n";
    for (size_t n = 1; n <= m; n <<= 1) {
        fout << n << "\t";
        Tester::generate2(matrix, n, m);
        int64 target = matrix[n - 1][0] + 1;
        fout << Tester::test(100, find_improved_n_plus_m, matrix, target) << "\n";
    }
        
    fout.close();

    return 0;
}