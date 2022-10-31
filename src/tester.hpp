#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <random>
#include <iomanip>
#include <fstream>
#include <Windows.h>
//
//namespace LEGACY {
//
//    std::mt19937 gen;
//    typedef long long int64;
//
//    constexpr size_t TESTS_FOR_PART = 10;
//    constexpr size_t NUMBERS = 10;
//    constexpr size_t TESTS = 10;
//    constexpr size_t UNIQUE_MATRICES = 1;
//
//    void ShowConsoleCursor(bool showFlag)
//    {
//        HANDLE out = GetStdHandle(STD_OUTPUT_HANDLE);
//
//        CONSOLE_CURSOR_INFO     cursorInfo;
//
//        GetConsoleCursorInfo(out, &cursorInfo);
//        cursorInfo.bVisible = showFlag; // set the cursor visibility
//        SetConsoleCursorInfo(out, &cursorInfo);
//    }
//
//    enum class Density {
//        sparse,
//        medium,
//        saturated
//    };
//
//    enum class Size {
//        tiny,
//        medium,
//        big,
//        large
//    };
//
//    enum class Output {
//        console,
//        txt
//    };
//
//    inline void generate_legacy(std::vector<std::vector<int64>>& matrix, size_t n, size_t m, Density density = Density::medium) {
//
//        size_t range = 1;
//        switch (density) {
//        case Density::sparse:
//            range = 10'000;
//            break;
//        case Density::medium:
//            range = 100;
//            break;
//        case Density::saturated:
//            range = 10;
//            break;
//        default:
//            break;
//        }
//
//        matrix.resize(n);
//        for (size_t i = 0; i < n; i++) {
//            matrix[i].resize(m);
//        }
//
//        matrix[0][0] = gen() % range;
//        for (size_t i = 1; i < n; i++) {
//            matrix[i][0] = matrix[i - 1][0] + gen() % range + 1;
//        }
//
//        for (size_t j = 1; j < m; j++) {
//            matrix[0][j] = matrix[0][j - 1] + gen() % range + 1;
//        }
//
//        for (size_t i = 1; i < n; i++) {
//            for (size_t j = 1; j < m; j++) {
//                matrix[i][j] = max(matrix[i - 1][j], matrix[i][j - 1]) + gen() % range + 1;
//            }
//        }
//    }
//
//    inline bool find_nlogm(std::vector<std::vector<int64>>& matrix, int64 k) {
//        size_t n = matrix.size();
//        size_t m = matrix[0].size();
//        if (n <= m) {
//            for (size_t i = 0; i < n; i++) {
//                int64 l = 0;
//                int64 r = m - 1;
//                while (l <= r) {
//                    int64 mid = (l + r) / 2;
//                    if (matrix[i][mid] == k) {
//                        return true;
//                    }
//                    else if (matrix[i][mid] < k) {
//                        l = mid + 1;
//                    }
//                    else {
//                        r = mid - 1;
//                    }
//                }
//            }
//            return false;
//        }
//        else {
//            for (size_t i = 0; i < m; i++) {
//                int64 l = 0;
//                int64 r = n - 1;
//                while (l <= r) {
//                    int64 mid = (l + r) / 2;
//                    if (matrix[mid][i] == k) {
//                        return true;
//                    }
//                    else if (matrix[mid][i] < k) {
//                        l = mid + 1;
//                    }
//                    else {
//                        r = mid - 1;
//                    }
//                }
//            }
//            return false;
//        }
//    }
//
//    inline bool find_n_plus_m(std::vector<std::vector<int64>>& matrix, int64 k) {
//        int64 i = 0;
//        int64 j = matrix[0].size() - 1;
//        while (i < matrix.size() && j >= 0) {
//            if (matrix[i][j] == k) {
//                return true;
//            }
//            else if (matrix[i][j] > k) {
//                j--;
//            }
//            else {
//                i++;
//            }
//        }
//        return false;
//    }
//
//    inline bool find_improved_nlogm(std::vector<std::vector<int64>>& matrix, int64 k) {
//        size_t n = matrix.size();
//        size_t m = matrix[0].size();
//        if (n <= m) {
//            size_t lower_bound = m - 1;
//            for (size_t i = 0; i < n;) {
//                int64 l = 0;
//                int64 r = lower_bound;
//                while (l <= r) {
//                    int64 mid = (l + r) / 2;
//                    if (matrix[i][mid] >= k) {
//                        r = mid - 1;
//                        lower_bound = mid;
//                    }
//                    else {
//                        l = mid + 1;
//                    }
//                }
//                if (matrix[i][lower_bound] == k) {
//                    return true;
//                }
//                else {
//                    i++;
//                }
//            }
//            return false;
//        }
//        else {
//            size_t lower_bound = n - 1;
//            for (size_t j = 0; j < m;) {
//                int64 l = 0;
//                int64 r = lower_bound;
//                while (l <= r) {
//                    int64 mid = (l + r) / 2;
//                    if (matrix[mid][j] >= k) {
//                        r = mid - 1;
//                        lower_bound = mid;
//                    }
//                    else {
//                        l = mid + 1;
//                    }
//                }
//                if (matrix[lower_bound][j] == k) {
//                    return true;
//                }
//                else {
//                    j++;
//                }
//            }
//            return false;
//        }
//    }
//
//    template <typename Func, typename... Args>
//    double test(size_t& found, size_t count, const Func& func, Args&... args) {
//        auto start = std::chrono::high_resolution_clock::now();
//        for (size_t i = 0; i < count; i++) {
//            found += static_cast<size_t>(func(args...));
//        }
//        auto end = std::chrono::high_resolution_clock::now();
//        std::chrono::duration<double> res = end - start;
//
//        return res.count() / static_cast<double>(count);
//    }
//
//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    // DESCRIPTION
//    /*
//        small       medium          big                 large
//        1x1'000     1x100'000       1x10'000'000        1x100'000'000
//        2x500       2x50'000        2x5'000'000         2x50'000'000
//        4x250       4x25'000        4x2'500'000         4x25'000'000
//        8x125       8x12'500        8x1'250'000         8x12'500'000
//        16x64       16x6'250        16x625'000          16x6'250'000
//        32x32       32x3'125        32x312'500          32x3'125'000
//        b.o.        50x2'000        50x200'000          50x2'000'000
//                    100x1'000       100x100'000         100x1'000'000
//                    200x500         200x50'000          200x500'000
//                    320x320         400x25'000          400x250'000
//                    b.o.            800x12'500          800x125'000
//                                    1'600x6'250         1'600x62'500
//                                    3'200x3'200         3'200x31'250
//                                    b.o.                5'000x20'000
//                                                        10'000x10'000
//                                                        b.o.
//
//        measurement system: 10 random numbers from each different part of matrix
//        small (2*2 parts)
//        medium (4*4 parts)
//        big, large (8*8 parts)
//
//        total tests complexity:
//        C(small) = 11 rectangles * 1 unique_matrices * 4 parts * 10 times * 10 numbers * 3 algorithms * 10 tests      = 13'200 tests
//        C(medium) = 19 rectangles * 1 unique_matrices * 16 parts * 10 times * 10 numbers * 3 algorithms * 10 tests    = 91'200 tests
//        C(big) = 25 rectangles * 1 unique_matrices * 64 parts * 10 times * 10 numbers * 3 algorithms * 10 tests       = 120'000 tests
//        C(large) = 29 rectangles * 1 unique_matrices * 64 parts * 10 times * 10 numbers * 3 algorithms * 10 tests     = 556'800 tests
//                                                                                                                      _______________
//                                                                                                                      = 781'200 tests
//    */
//    // DESCRIPTION
//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    std::vector<std::pair<size_t, size_t>> tiny_rectangles = {
//        {1, 1'000},
//        {2, 500},
//        {4, 250},
//        {8, 125},
//        {16, 64},
//        {32, 32},
//        {64, 16},
//        {125, 8},
//        {250, 4},
//        {500, 2},
//        {1'000, 1}
//    };
//    const size_t tiny_parts = tiny_rectangles.size() * 4 - 4 * 3;
//
//    std::vector<std::pair<size_t, size_t>> medium_rectangles = {
//        {1, 100'000},
//        //{2, 50'000},
//        {4, 25'000},
//        //{8, 12'500},
//        {16, 6'250},
//        //{32, 3'125},
//        {50, 2'000},
//        {100, 1'000},
//        {200, 500},
//        {320, 320},
//        {500, 200},
//        {1'000, 100},
//        {2'000, 50},
//        //{3'125, 32},
//        {6'250, 16},
//        //{12'500, 8},
//        {25'000, 4},
//        //{50'000, 2},
//        {100'000, 1}
//    };
//    const size_t medium_parts = medium_rectangles.size() * 16 - 2 * 15;
//
//    std::vector<std::pair<size_t, size_t>> big_rectangles = {
//        {1, 10'000'000},
//        //{2, 5'000'000},
//        //{4, 2'500'000},
//        //{8, 1'250'000},
//        {16, 625'000},
//        //{32, 312'500},
//        //{50, 200'000},
//        {100, 100'000},
//        //{200, 50'000},
//        {400, 25'000},
//        //{800, 12'500},
//        {1'600, 6'250},
//        {2'500, 4'000},
//        {3'200, 3'200},
//        {4'000, 2'500},
//        {6'250, 1'600},
//        //{12'500, 800},
//        {25'000, 400},
//        //{50'000, 200},
//        //{100'000, 100},
//        //{200'000, 50},
//        //{312'500, 32},
//        //{625'000, 16},
//        //{1'250'000, 8},
//        //{2'500'000, 4},
//        //{5'000'000, 2},
//        {10'000'000, 1}
//    };
//    const size_t big_parts = big_rectangles.size() * 64 - 2 * 63;
//
//    std::vector<std::pair<size_t, size_t>> large_rectangles = {
//        {1, 100'000'000},
//        //{2, 50'000'000},
//        //{4, 25'000'000},
//        {8, 12'500'000},
//        //{16, 6'250'000},
//        //{32, 3'125'000},
//        {50, 2'000'000},
//        //{100, 1'000'000},
//        //{200, 500'000},
//        {400, 250'000},
//        //{800, 125'000},
//        //{1'600, 62'500},
//        {3'200, 32'000},
//        {5'000, 20'000},
//        {8'000, 12'500},
//        {10'000, 10'000},
//        {12'500, 8'000},
//        {20'000, 5'000},
//        {32'000, 3'200},
//        //{62'500, 1'600},
//        //{125'000, 800},
//        {250'000, 400},
//        //{500'000, 200},
//        //{1'000'000, 100},
//        {2'000'000, 50},
//        //{3'125'000, 32},
//        //{6'250'000, 16},
//        {12'500'000, 8},
//        //{25'000'000, 4},
//        //{50'000'000, 2},
//        {100'000'000, 1}
//    };
//    const size_t large_parts = large_rectangles.size() * 64 - 2 * 63;
//
//    std::string get_file_name(Size test_set_size, Density density) {
//        std::string res;
//        switch (test_set_size) {
//        case Size::tiny:
//            res += "Tiny_";
//            break;
//        case Size::medium:
//            res += "Medium_";
//            break;
//        case Size::big:
//            res += "Big_";
//            break;
//        case Size::large:
//            res += "Large_";
//            break;
//        default:
//            break;
//        }
//
//        switch (density) {
//        case Density::sparse:
//            res += "sparce_";
//            break;
//        case Density::medium:
//            res += "medium_";
//            break;
//        case Density::saturated:
//            res += "saturated_";
//            break;
//        default:
//            break;
//        }
//
//        res += "benchmark.txt";
//        return res;
//    }
//
//    void run_test_set_legacy(std::vector<std::vector<int64>>& matrix, Size test_set_size,
//        Output output, Density density = Density::medium) {
//
//        ShowConsoleCursor(false);
//
//        std::basic_ostream<char>* out(&std::cout);
//
//        std::ofstream fout;
//        if (output == Output::txt) {
//            fout.open(get_file_name(test_set_size, density));
//            if (!fout.is_open()) {
//                std::cout << "Can't open txt" << std::endl;
//                exit(EXIT_FAILURE);
//            }
//            out = &fout;
//        }
//
//        *out << std::fixed << std::setprecision(9);
//
//        std::shared_ptr<std::vector<std::pair<size_t, size_t>>> rect = nullptr;
//        size_t partition = 0;
//        size_t total_partitions = 0;
//
//        switch (test_set_size) {
//        case Size::tiny:
//            rect = std::make_shared<std::vector<std::pair<size_t, size_t>>>(tiny_rectangles);
//            partition = 2;
//            total_partitions = tiny_parts;
//            break;
//        case Size::medium:
//            rect = std::make_shared<std::vector<std::pair<size_t, size_t>>>(medium_rectangles);
//            partition = 4;
//            total_partitions = medium_parts;
//            break;
//        case Size::big:
//            rect = std::make_shared<std::vector<std::pair<size_t, size_t>>>(big_rectangles);
//            partition = 8;
//            total_partitions = big_parts;
//            break;
//        case Size::large:
//            rect = std::make_shared<std::vector<std::pair<size_t, size_t>>>(large_rectangles);
//            partition = 8;
//            total_partitions = large_parts;
//            break;
//        default:
//            break;
//        }
//
//        size_t all_suites = UNIQUE_MATRICES * total_partitions;
//        size_t current_suites = 0;
//
//        size_t vin[3] = {};
//        std::vector<std::vector<size_t>> vin_per_part(partition * partition, { 0, 0, 0 });
//
//        for (size_t shape = 0; shape < rect->size(); shape++) {
//            size_t n = (*rect)[shape].first;
//            size_t m = (*rect)[shape].second;
//            *out << "Shape: n = " << n << "; m = " << m << " ==================================================\n";
//            for (size_t unique_matrix = 0; unique_matrix < UNIQUE_MATRICES; unique_matrix++) {
//                generate_legacy(matrix, n, m, density);
//                std::vector<std::pair<int64, int64>> parts;
//
//                if (n >= partition && m >= partition) {
//                    for (size_t row = 0; row < partition; row++) {
//                        for (size_t col = 0; col < partition; col++) {
//                            parts.push_back({
//                                matrix[row * (n / partition)][col * (m / partition)],
//                                matrix[(row + 1) * (n / partition) - 1][(col + 1) * (m / partition) - 1]
//                                });
//                        }
//                    }
//                }
//                else {
//                    parts.push_back({ matrix[0][0], matrix[n - 1][m - 1] });
//                }
//
//                for (size_t p = 0; p < parts.size(); p++) {
//                    if (output != Output::console) {
//                        current_suites++;
//                        std::cout << "\r                \r" <<
//                            min(100.0, (static_cast<double>(current_suites) / all_suites) * 100.0) << '%';
//                    }
//
//                    *out << "Part " << p + 1 << '\n';
//                    auto& part = parts[p];
//                    double average_time_nlogm = 0.0;
//                    double average_time_n_plus_m = 0.0;
//                    double average_time_improved_nlogm = 0.0;
//
//                    size_t total_found = 0;
//
//                    for (size_t part_test = 0; part_test < TESTS_FOR_PART; part_test++) {
//                        for (size_t number = 0; number < NUMBERS; number++) {
//                            int64 random_number = gen() % (part.second - part.first + 1) + part.first;
//                            average_time_nlogm += test(total_found, TESTS, find_nlogm, matrix, random_number);
//                            average_time_n_plus_m += test(total_found, TESTS, find_n_plus_m, matrix, random_number);
//                            average_time_improved_nlogm += test(total_found, TESTS, find_improved_nlogm, matrix, random_number);
//                        }
//
//                        if (average_time_nlogm < average_time_n_plus_m && average_time_nlogm < average_time_improved_nlogm) {
//                            vin[0]++;
//                            vin_per_part[p][0]++;
//                        }
//                        else if (average_time_n_plus_m < average_time_nlogm && average_time_n_plus_m < average_time_improved_nlogm) {
//                            vin[1]++;
//                            vin_per_part[p][1]++;
//                        }
//                        else {
//                            vin[2]++;
//                            vin_per_part[p][2]++;
//                        }
//                    }
//
//                    *out << "found: " <<
//                        static_cast<double>(total_found / 3) / (NUMBERS * TESTS_FOR_PART * TESTS) * 100.0 << "%\n";
//                    *out << "nlogm:\t" << average_time_nlogm / (NUMBERS * TESTS_FOR_PART) << "s.\n";
//                    *out << "n + m:\t" << average_time_n_plus_m / (NUMBERS * TESTS_FOR_PART) << "s.\n";
//                    *out << "nlogm+:\t" << average_time_improved_nlogm / (NUMBERS * TESTS_FOR_PART) << "s.\n";
//                    *out << std::endl;
//                }
//            }
//
//            *out << "nlogm/n+m/nlogm+\n";
//            for (size_t i = 0; i < partition; i++) {
//                for (size_t j = 0; j < partition; j++) {
//                    for (size_t k = 0; k < 3; k++) {
//                        *out << vin_per_part[i * partition + j][k];
//                        if (k != 2) {
//                            *out << '/';
//                        }
//                    }
//                    *out << '\t';
//                }
//                *out << '\n';
//            }
//            *out << "\n\n";
//
//            /*std::cout <<
//                "nlogm\n"
//                "n+m\n"
//                "nlogm+\n";
//            for (size_t k = 0; k < partition; k++) {
//                for (size_t j = 0; j < 3; j++) {
//                    for (size_t i = 0; i < partition; i++) {
//                        *out << vin_per_part[i][j] << "\t\t";
//                    }
//                    *out << '\n';
//                }
//                *out << '\n';
//            }
//            *out << "\n\n";*/
//
//            vin_per_part.assign(partition * partition, { 0, 0, 0 });
//        }
//
//        *out << "\n\nTOTAL STATISTICS\n";
//        *out << "nlogm is better in:\t" << vin[0] << " cases;\n";
//        *out << "n + m is better in:\t" << vin[1] << " cases;\n";
//        *out << "nlogm+ is better in:\t" << vin[2] << " cases;\n";
//        *out << "#####################################################################################################\n\n";
//
//        for (size_t i = 0; i < 3; i++) {
//            vin[i] = 0;
//        }
//
//        ShowConsoleCursor(true);
//
//        if (output != Output::console) {
//            std::cout << '\n';
//        }
//    }
//
//}

typedef long long int64;

class Tester {
private:
    static size_t moc_found; // nedeed for compiler to do not optimize the code without side effects.

public:

    static void generate(std::vector<std::vector<int64>>& matrix, size_t n, size_t m) {
        matrix.resize(n);
        for (size_t i = 0; i < n; i++) {
            matrix[i].resize(m);
        }

        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < m; j++) {
                matrix[i][j] = ((m / n) * i + j) * 2;
            }
        }
    }

    static void generate2(std::vector<std::vector<int64>>& matrix, size_t n, size_t m) {
        matrix.resize(n);
        for (size_t i = 0; i < n; i++) {
            matrix[i].resize(m);
        }

        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < m; j++) {
                matrix[i][j] = ((m / n) * i * j) * 2;
            }
        }
    }

    template <typename Func, typename... Args>
    static size_t test(size_t count, const Func& func, Args&... args) {
        auto start = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < count; i++) {
            moc_found += static_cast<size_t>(func(args...));
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

        return nanoseconds.count() / count;
    }



};
