#include <catch2/catch_all.hpp>

#include <ccd.hpp>

#include <gmp.h>
#include <filesystem>
#include <fstream>

enum Case { POINT_TRIANGLE, EDGE_EDGE };

std::pair<Eigen::MatrixXd, std::vector<bool>>
read_rational_csv(const std::string& filename);

TEST_CASE("Dataset", "[ccd][point-triangle][edge-edge][dataset][!mayfail]")
{
    using namespace ccd;
    namespace fs = std::filesystem;

    const std::string folder_name = GENERATE(
        "erleben-sliding-spike", "erleben-spike-wedge", "erleben-sliding-wedge",
        "erleben-wedge-crack", "erleben-spike-crack", "erleben-wedges",
        "erleben-cube-cliff-edges", "erleben-spike-hole",
        "erleben-cube-internal-edges", "erleben-spikes", "unit-tests", "chain",
        "cow-heads", "golf-ball", "mat-twist");
    const Case test_case = GENERATE(POINT_TRIANGLE, EDGE_EDGE);

    const fs::path path = fs::path(SAMPLE_QUERIES_DIR) / folder_name
        / (test_case == POINT_TRIANGLE ? "vertex-face" : "edge-edge");
    REQUIRE(fs::exists(path));
    REQUIRE(fs::is_directory(path));

    for (const auto& f : fs::directory_iterator(path)) {
        if (f.path().extension() != ".csv") {
            continue;
        }

        CAPTURE(f.path().string());
        const auto& [queries, results] = read_rational_csv(f.path().string());
        assert(queries.rows() % 8 == 0 && queries.cols() == 3);

        const int n_queries = queries.rows() / 8;
        for (int i = 0; i < n_queries; i++) {
            const Eigen::Matrix<double, 8, 3> query =
                queries.middleRows<8>(8 * i);
            const bool expected_result = results[i * 8];

            bool result;
            double toi;
            switch (test_case) {
            case POINT_TRIANGLE:
                result = point_triangle_ccd(
                    query.row(0), query.row(1), query.row(2), query.row(3),
                    query.row(4), query.row(5), query.row(6), query.row(7),
                    toi);
                break;
            case EDGE_EDGE:
                result = edge_edge_ccd(
                    query.row(0), query.row(1), query.row(2), query.row(3),
                    query.row(4), query.row(5), query.row(6), query.row(7),
                    toi);
            }

            CHECK(result == expected_result);
        }
    }
}

double
double_from_rational_strings(const std::string& num, const std::string& denom)
{
    static mpq_t value;
    std::string tmp = num + "/" + denom;
    mpq_set_str(value, tmp.c_str(), 10);
    return mpq_get_d(value);
}

std::pair<Eigen::MatrixXd, std::vector<bool>>
read_rational_csv(const std::string& filename)
{
    Eigen::MatrixXd queries;
    std::vector<bool> results;

    // be careful, there are n lines which means there are n/8 queries, but has
    // n results, which means results are duplicated
    std::vector<Eigen::Vector3d> vs;
    std::ifstream f(filename);
    assert(f.is_open());

    while (f) {
        std::string s;
        if (!getline(f, s))
            break;

        if (s[0] == '#') {
            continue;
        }

        std::istringstream ss(s);
        // the first six are one vetex, the seventh is the result
        std::array<std::string, 7> record;
        int c = 0;
        while (ss) {
            std::string line;
            if (!getline(ss, line, ','))
                break;
            record[c++] = line;
        }
        vs.emplace_back(
            double_from_rational_strings(record[0], record[1]),
            double_from_rational_strings(record[2], record[3]),
            double_from_rational_strings(record[4], record[5]));
        results.push_back(std::stoi(record[6]));
    }
    queries.resize(vs.size(), 3);
    for (int i = 0; i < vs.size(); i++) {
        for (int j = 0; j < 3; j++) {
            queries(i, j) = vs[i][j];
        }
    }

    return std::make_pair(queries, results);
}