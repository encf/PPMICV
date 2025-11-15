#include "P1.hpp"
#include "P0.hpp"
#include "utils.hpp"
#include <random> // Required for random number generation
#include <iostream>
#include <utility>
#include <numeric>

// Constructor and original getters remain the same
P1::P1() : type_(1), code_(410620009), start_date_(20240101), end_date_(20251231), max_cost_(10000*1024), ratio_(1) {}
uint32_t P1::getType() const { return type_; }
uint32_t P1::getCode() const { return code_; }
uint32_t P1::getStartDate() const { return start_date_; }
uint32_t P1::getEndDate() const { return end_date_; }
uint32_t P1::getMaxCost() const { return max_cost_; }


// --- Implementation of New Secret Sharing Functionality ---

P1RandomShare P1::createSecretShares() {
    // Set up the random number generator
    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_int_distribution<uint32_t> dist;

    // Create a struct to hold the random share that will be returned
    P1RandomShare random_share;

    // 1. Share for 'type'
    random_share.type = dist(eng);
    retained_type_ = type_ - random_share.type;

    // 2. Share for 'code'
    random_share.code = dist(eng);
    retained_code_ = code_ - random_share.code;

    // 3. Share for 'start_date'
    random_share.start_date = dist(eng);
    retained_start_date_ = start_date_ - random_share.start_date;

    // 4. Share for 'end_date'
    random_share.end_date = dist(eng);
    retained_end_date_ = end_date_ - random_share.end_date;

    // 5. Share for 'max_cost'
    random_share.max_cost = dist(eng);
    retained_max_cost_ = max_cost_ - random_share.max_cost;
    
    return random_share;
}

// --- Implementation of Getters for the retained share ---

uint32_t P1::getRetainedType() const { return retained_type_; }
uint32_t P1::getRetainedCode() const { return retained_code_; }
uint32_t P1::getRetainedStartDate() const { return retained_start_date_; }
uint32_t P1::getRetainedEndDate() const { return retained_end_date_; }
uint32_t P1::getRetainedMaxCost() const { return retained_max_cost_; }



void P1::receiveP0Share(const std::vector<std::vector<uint32_t>>& p0_share) {
    received_p0_share_ = p0_share;
}

const std::vector<std::vector<uint32_t>>& P1::getReceivedP0Share() const {
    return received_p0_share_;
}


void P1::receiveTriplesShare(const std::vector<uint32_t>& a_share,
                             const std::vector<uint32_t>& b_share,
                             const std::vector<uint32_t>& c_share) {
    triples_a_share_ = a_share;
    triples_b_share_ = b_share;
    triples_c_share_ = c_share;
}

const std::vector<uint32_t>& P1::getTriplesShareA() const { return triples_a_share_; }
const std::vector<uint32_t>& P1::getTriplesShareB() const { return triples_b_share_; }
const std::vector<uint32_t>& P1::getTriplesShareC() const { return triples_c_share_; }

std::vector<std::vector<uint32_t>> P1::prepareComparisonShares(
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
            std::cout << "\n--- DEBUG P1 (i=" << i << ") ---" << std::endl;
            std::cout << "  Step 1: Received x0 share: " << x0_share[i] << std::endl;
            std::cout << "  Step 2: t=" << t_values[i] << ", x after t-mask: " << x << std::endl;
        }

        std::vector<uint32_t> masked_vector(34);
        uint32_t u = t_values[i];
        masked_vector[0] = u + 3 * x;
        if (is_debug) {
            print_vector_full("Step 3/4.1 : masked vector", masked_vector);
        }

        std::vector<uint32_t> u_vec(33);
        uint32_t temp = x;
        uint32_t temp2 = -temp;
        u_vec[0] = x;
        for(int j=1; j<33; ++j) {
            temp2 = temp2 >> 1;
            u_vec[j] = -temp2;
        }
        if (is_debug) {
            print_vector_full("Step 3/4.2 : u_vec", u_vec);
        }

        uint32_t suffix_sum = 0;
        for (int j = 32; j >= 0; --j) {
            suffix_sum += u_vec[j];
            masked_vector[j + 1] = suffix_sum;
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

void P1::receiveComparisonResultShare(const std::vector<uint32_t>& r_prime_share) {
    r_prime_share_p1_ = r_prime_share;
}

std::vector<uint32_t> P1::computeFinalComparisonResult(const std::vector<uint32_t>& t_values) {
    size_t n = t_values.size();
    std::vector<uint32_t> final_result_share(n);
    for (size_t i = 0; i < n; ++i) {
        uint32_t t = t_values[i];
        uint32_t r_prime = r_prime_share_p1_[i];
        // r = r' - 2*t*r'
        final_result_share[i] = r_prime - 2 * t * r_prime;
    }
    return final_result_share;
}

const std::vector<uint32_t>& P1::getComparisonResultShare() const {
    return r_prime_share_p1_;
}


void P1::receiveEqualityTValues(const std::vector<uint32_t>& t0, const std::vector<uint32_t>& t1) {
    t0_share_ = t0;
    t1_share_ = t1;
}

void P1::receiveEqualityResultShare(const std::vector<uint32_t>& d1_share) {
    d1_final_share_ = d1_share;
}


std::vector<uint32_t> P1::computeFinalEqualityResult() {
    size_t n = d1_final_share_.size();
    std::vector<uint32_t> result_share(n);
    for (size_t i = 0; i < n; ++i) {
        uint32_t t = t0_share_[i] ^ t1_share_[i];
        uint32_t d1 = d1_final_share_[i];
        result_share[i] = -(d1 - 2 * t * d1);
    }
    return result_share;
}


std::pair<std::vector<uint32_t>, std::vector<uint32_t>> P1::compute_d_e_shares(
    const std::vector<uint32_t>& x1,
    const std::vector<uint32_t>& y1
) {
    size_t n = x1.size();
    // P1 uses the triple shares it received, which are stored in triples_a_share_ etc.
    const auto& a1 = this->getTriplesShareA();
    const auto& b1 = this->getTriplesShareB();

    std::vector<uint32_t> d1(n), e1(n);
    for (size_t i = 0; i < n; ++i) {
        d1[i] = x1[i] - a1[i];
        e1[i] = y1[i] - b1[i];
    }
    return std::make_pair(d1, e1);
}

std::vector<uint32_t> P1::compute_z1_share(
    const std::vector<uint32_t>& d_plain,
    const std::vector<uint32_t>& e_plain
) {
    size_t n = d_plain.size();
    const auto& a1 = this->getTriplesShareA();
    const auto& b1 = this->getTriplesShareB();
    const auto& c1 = this->getTriplesShareC();
    
    std::vector<uint32_t> z1(n);
    for (size_t i = 0; i < n; ++i) {
        uint32_t d = d_plain[i];
        uint32_t e = e_plain[i];
        // P1's formula does NOT include the public d*e term
        z1[i] = (d * b1[i]) + (e * a1[i]) + c1[i];
    }
    return z1;
}

uint32_t P1::computeDotProductShare(
    const std::vector<uint32_t>& d_plain,
    const std::vector<uint32_t>& e_plain
) {
    // Compute the shares of the element-wise product z1 = d*b1 + e*a1 + c1
    std::vector<uint32_t> z1 = this->compute_z1_share(d_plain, e_plain);

    // Locally sum the shares of the element-wise product
    return std::accumulate(z1.begin(), z1.end(), (uint32_t)0);
}


void P1::sharePolicyData(P0& p0) {
    std::vector<uint32_t> all_data = {type_, code_, start_date_, end_date_, max_cost_, ratio_};
    std::vector<uint32_t> shares_0(6), shares_1(6);

    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_int_distribution<uint32_t> dist;

    for (int i = 0; i < 6; ++i) {
        shares_1[i] = dist(eng);
        shares_0[i] = all_data[i] - shares_1[i];
    }

    type_share_1_.push_back(shares_1[0]);
    code_share_1_.push_back(shares_1[1]);
    start_date_share_1_.push_back(shares_1[2]);
    end_date_share_1_.push_back(shares_1[3]);
    max_cost_share_1_ = shares_1[4];
    ratio_share_1_ = shares_1[5];

    p0.receivePolicyDataShares(shares_0);
}

void P1::receiveMedicalDataShares(const std::vector<std::vector<uint32_t>>& medical_shares) {
    type_m_share_1_ = medical_shares[0];
    code_m_share_1_ = medical_shares[1];
    date_m_share_1_ = medical_shares[2];
    cost_m_share_1_ = medical_shares[3];
}

void P1::expandPolicyData(size_t n) {
    type_share_1_.resize(n, type_share_1_[0]);
    code_share_1_.resize(n, code_share_1_[0]);
    start_date_share_1_.resize(n, start_date_share_1_[0]);
    end_date_share_1_.resize(n, end_date_share_1_[0]);
}
