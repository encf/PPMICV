#include "P0.hpp"
#include "P1.hpp"
#include "utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <random>
#include <utility>
#include <numeric>


bool P0::readAndProcessCsv(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "error : cannot open file " << filename << std::endl;
        return false;
    }

    dataColumns_.clear();
    dataColumns_.resize(4); 

    std::string line;
    if (std::getline(file, line)) {}

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> rowData;

        while (std::getline(ss, cell, ',')) {
            rowData.push_back(cell);
        }

        if (rowData.size() == 5) {
            try {
                dataColumns_[0].push_back(static_cast<uint32_t>(std::stoul(rowData[1])));
                dataColumns_[1].push_back(static_cast<uint32_t>(std::stoul(rowData[2])));
                dataColumns_[2].push_back(static_cast<uint32_t>(std::stoul(rowData[3])));

                double lastColValue = std::stod(rowData[4]);
                lastColValue *= 1024.0; 
                dataColumns_[3].push_back(static_cast<uint32_t>(std::round(lastColValue)));

            } catch (const std::exception& e) {
                std::cerr << "warning : " << line << std::endl;
            }
        }
    }

    file.close();
    return true;
}

const std::vector<std::vector<uint32_t>>& P0::getDataColumns() const {
    return dataColumns_;
}

void P0::printData() const {
    if (dataColumns_.empty() || dataColumns_[0].empty()) {
        std::cout << "data is empty." << std::endl;
        return;
    }

    size_t num_rows = dataColumns_[0].size();
    size_t num_cols = dataColumns_.size();

    for (size_t i = 0; i < num_rows; ++i) {
        for (size_t j = 0; j < num_cols; ++j) {
            std::cout << dataColumns_[j][i] << "\t";
        }
        std::cout << std::endl;
    }
}


const std::vector<std::vector<uint32_t>>& P0::getRetainedShareColumns() const {
    return retainedColumnShares_;
}

std::vector<std::vector<uint32_t>> P0::createSecretShares() {
    if (dataColumns_.empty() || dataColumns_[0].empty()) {
        std::cerr << "error : data is empty" << std::endl;
        return {};
    }

    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_int_distribution<uint32_t> dist;

    size_t num_cols = dataColumns_.size();
    size_t num_rows = dataColumns_[0].size();

    std::vector<std::vector<uint32_t>> randomColumnShares(num_cols);
    retainedColumnShares_.assign(num_cols, {});

    for (size_t j = 0; j < num_cols; ++j) {
        randomColumnShares[j].reserve(num_rows);
        retainedColumnShares_[j].reserve(num_rows);

        for (size_t i = 0; i < num_rows; ++i) {
            uint32_t random_val = dist(eng);
            uint32_t original_val = dataColumns_[j][i];
            uint32_t retained_val = original_val - random_val;

            randomColumnShares[j].push_back(random_val);
            retainedColumnShares_[j].push_back(retained_val);
        }
    }

    return randomColumnShares;
}

void P0::setArrays(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b, const std::vector<uint32_t>& c) { array_a_ = a; array_b_ = b; array_c_ = c; }
const std::vector<uint32_t>& P0::getArrayA() const { return array_a_; }
const std::vector<uint32_t>& P0::getArrayB() const { return array_b_; }
const std::vector<uint32_t>& P0::getArrayC() const { return array_c_; }


void P0::receiveP1Share(const P1RandomShare& p1_share) {
    received_p1_share_ = p1_share;
}

const P1RandomShare& P0::getReceivedP1Share() const {
    return received_p1_share_;
}

std::vector<std::vector<uint32_t>> P0::prepareComparisonShares(
    const std::vector<uint32_t>& x0_share,
    const std::vector<uint32_t>& t_values,
    const std::vector<std::vector<uint32_t>>& permutations,
    const std::vector<uint32_t>& mul_rands,
    int debug_index
) {
    size_t n = x0_share.size();
    std::vector<std::vector<uint32_t>> all_final_vectors(n);

    for (size_t i = 0; i < n; ++i) {
        uint32_t x = x0_share[i];
        if (t_values[i] == 1) {
            x = -x;
        }

        bool is_debug = (i == debug_index);
        is_debug = 0;
        if (is_debug) {
            std::cout << "\n--- DEBUG P0 (i=" << i << ") ---" << std::endl;
            std::cout << "  Step 1: Received x0 share: " << x0_share[i] << std::endl;
            std::cout << "  Step 2: t=" << t_values[i] << ", x after t-mask: " << x << std::endl;
        }

        std::vector<uint32_t> masked_vector(34);
        uint32_t u = (t_values[i] == 1) ? -1 : 1;
        u = u - t_values[i];
        masked_vector[0] = u + 3 * x - 1;
        if (is_debug) {
            print_vector_full("Step 3/4.1 : masked vector", masked_vector);
        }


        std::vector<uint32_t> u_vec(33);
        uint32_t temp = x;
        for(int j=0; j<33; ++j) {
            u_vec[j] = temp;
            temp = temp >> 1;
        }
        if (is_debug) {
            print_vector_full("Step 3/4.2 : u_vec", u_vec);
        }



        uint32_t suffix_sum = 0;
        for (int j = 32; j >= 0; --j) {
            suffix_sum += u_vec[j];
            masked_vector[j + 1] = suffix_sum - 1;
        }

        if (is_debug) {
            print_vector_full("Step 3/4 (Suffix Sum): Vector before mult-mask", masked_vector);
        }

        for (int j = 0; j < 34; ++j) {
            masked_vector[j] *= mul_rands[j];
        }

        if (is_debug) {
            print_vector_full("Step 5a: Vector after mult-mask", masked_vector);
        }

        const auto& current_permutation = permutations[i];
        std::vector<uint32_t> shuffled_vector(34);
        for (int j = 0; j < 34; ++j) {
            shuffled_vector[j] = masked_vector[current_permutation[j]];
        }
        
        if (is_debug) {
            print_vector_full("Step 5b: Final shuffled vector sent to P2", shuffled_vector);
        }
        
        all_final_vectors[i] = shuffled_vector;
    }
    return all_final_vectors;
}

void P0::receiveComparisonResultShare(const std::vector<uint32_t>& r_prime_share) {
    r_prime_share_p0_ = r_prime_share;
}

std::vector<uint32_t> P0::computeFinalComparisonResult(const std::vector<uint32_t>& t_values) {
    size_t n = t_values.size();
    std::vector<uint32_t> final_result_share(n);
    for (size_t i = 0; i < n; ++i) {
        uint32_t t = t_values[i];
        uint32_t r_prime = r_prime_share_p0_[i];
        // r = t + r' - 2*t*r'
        final_result_share[i] = t + r_prime - 2 * t * r_prime;
    }
    return final_result_share;
}

const std::vector<uint32_t>& P0::getComparisonResultShare() const {
    return r_prime_share_p0_;
}



void P0::receiveEqualityTValues(const std::vector<uint32_t>& t0, const std::vector<uint32_t>& t1) {
    t0_share_ = t0;
    t1_share_ = t1;
}

void P0::receiveEqualityResultShare(const std::vector<uint32_t>& d0_share) {
    d0_final_share_ = d0_share;
}

// P0's version of the final calculation
std::vector<uint32_t> P0::computeFinalEqualityResult() {
    size_t n = d0_final_share_.size();
    std::vector<uint32_t> result_share(n);
    for (size_t i = 0; i < n; ++i) {
        // Step 6: P0 computes t = t0 ^ t1
        uint32_t t = t0_share_[i] ^ t1_share_[i];
        // Step 7: P0 computes r0 = 1 - (t + d0 - 2*t*d0)
        uint32_t d0 = d0_final_share_[i];
        result_share[i] = 1 - (t + d0 - 2 * t * d0);
    }
    return result_share;
}


std::pair<std::vector<uint32_t>, std::vector<uint32_t>> P0::compute_d_e_shares(
    const std::vector<uint32_t>& x0,
    const std::vector<uint32_t>& y0
) {
    size_t n = x0.size();
    // P0 uses the triple shares it received, which are stored in array_a_ etc.
    const auto& a0 = this->getArrayA();
    const auto& b0 = this->getArrayB();

    std::vector<uint32_t> d0(n), e0(n);
    for (size_t i = 0; i < n; ++i) {
        d0[i] = x0[i] - a0[i];
        e0[i] = y0[i] - b0[i];
    }
    return std::make_pair(d0, e0);
}


std::vector<uint32_t> P0::compute_z0_share(
    const std::vector<uint32_t>& d_plain,
    const std::vector<uint32_t>& e_plain
) {
    size_t n = d_plain.size();
    const auto& a0 = this->getArrayA();
    const auto& b0 = this->getArrayB();
    const auto& c0 = this->getArrayC();

    std::vector<uint32_t> z0(n);
    for (size_t i = 0; i < n; ++i) {
        uint32_t d = d_plain[i];
        uint32_t e = e_plain[i];
        // P0's formula includes the public d*e term
        z0[i] = (d * e) + (d * b0[i]) + (e * a0[i]) + c0[i];
    }
    return z0;
}

uint32_t P0::computeDotProductShare(
    const std::vector<uint32_t>& d_plain,
    const std::vector<uint32_t>& e_plain
) {
    // Step 3: Compute the shares of the element-wise product z0 = d*e + d*b0 + e*a0 + c0
    std::vector<uint32_t> z0 = this->compute_z0_share(d_plain, e_plain);

    // Step 5: Locally sum the shares of the element-wise product
    return std::accumulate(z0.begin(), z0.end(), (uint32_t)0);
}

void P0::readAndShareMedicalData(const std::string& filename, P1& p1) {
    std::vector<uint32_t> type_m, code_m, date_m, cost_m;
    std::ifstream file(filename);
    std::string line;
    std::getline(file, line); // Skip header

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> rowData;
        while (std::getline(ss, cell, ',')) {
            rowData.push_back(cell);
        }
        if (rowData.size() == 5) {
            type_m.push_back(stoul(rowData[1]));
            code_m.push_back(stoul(rowData[2]));
            date_m.push_back(stoul(rowData[3]));
            
            double lastColValue = std::stod(rowData[4]);
            lastColValue *= 1024.0;
            cost_m.push_back(static_cast<uint32_t>(std::round(lastColValue)));
            //cost_m.push_back(stoul(rowData[4]));
        }
    }
    file.close();

    size_t n = type_m.size();
    type_m_share_0_.resize(n); code_m_share_0_.resize(n);
    date_m_share_0_.resize(n); cost_m_share_0_.resize(n);

    std::vector<uint32_t> type_m_share_1(n), code_m_share_1(n), date_m_share_1(n), cost_m_share_1(n);
    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_int_distribution<uint32_t> dist;

    for (size_t i = 0; i < n; ++i) {
        type_m_share_0_[i] = dist(eng); type_m_share_1[i] = type_m[i] - type_m_share_0_[i];
        code_m_share_0_[i] = dist(eng); code_m_share_1[i] = code_m[i] - code_m_share_0_[i];
        date_m_share_0_[i] = dist(eng); date_m_share_1[i] = date_m[i] - date_m_share_0_[i];
        cost_m_share_0_[i] = dist(eng); cost_m_share_1[i] = cost_m[i] - cost_m_share_0_[i];
    }

    p1.receiveMedicalDataShares({type_m_share_1, code_m_share_1, date_m_share_1, cost_m_share_1});
}


void P0::receivePolicyDataShares(const std::vector<uint32_t>& policy_shares) {
    type_share_0_.push_back(policy_shares[0]);
    code_share_0_.push_back(policy_shares[1]);
    start_date_share_0_.push_back(policy_shares[2]);
    end_date_share_0_.push_back(policy_shares[3]);
    max_cost_share_0_ = policy_shares[4];
    ratio_share_0_ = policy_shares[5];
}

void P0::expandPolicyData(size_t n) {
    type_share_0_.resize(n, type_share_0_[0]);
    code_share_0_.resize(n, code_share_0_[0]);
    start_date_share_0_.resize(n, start_date_share_0_[0]);
    end_date_share_0_.resize(n, end_date_share_0_[0]);
}
