#include "P2.hpp"
#include "P0.hpp" 
#include "P1.hpp"
#include "utils.hpp"
#include <random>
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <utility>


void P2::generateAndDistributeTriples(int n, P0& p0, P1& p1) {
    std::cout << "P2: generate " << n << " triples..." << std::endl;

    a_.resize(n);
    b_.resize(n);
    c_.resize(n);

    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_int_distribution<uint32_t> dist;

    for (int i = 0; i < n; ++i) {
        a_[i] = dist(eng);
        b_[i] = dist(eng);
    }

    for (int i = 0; i < n; ++i) {
        c_[i] = a_[i] * b_[i];
    }
    std::cout << "P2: triples (a, b, c)" << std::endl;


    std::vector<uint32_t> a0(n), b0(n), c0(n); 
    std::vector<uint32_t> a1(n), b1(n), c1(n);

    for (int i = 0; i < n; ++i) {
        a0[i] = dist(eng);
        a1[i] = a_[i] - a0[i];

        b0[i] = dist(eng);
        b1[i] = b_[i] - b0[i];

        c0[i] = dist(eng);
        c1[i] = c_[i] - c0[i];
    }

    p0.setArrays(a0, b0, c0);
    p1.receiveTriplesShare(a1, b1, c1);
}

// Getter
const std::vector<uint32_t>& P2::getOriginalA() const { return a_; }
const std::vector<uint32_t>& P2::getOriginalB() const { return b_; }
const std::vector<uint32_t>& P2::getOriginalC() const { return c_; }

std::vector<uint32_t> P2::performSecureComparison(
    const std::vector<uint32_t>& x0, P0& p0,
    const std::vector<uint32_t>& x1, P1& p1,
    int debug_index
) {
    size_t n = x0.size();
    std::random_device rd;
    std::mt19937 eng(rd());

    // --- Step 1 & 5: Generate all shared randomness ---
    std::uniform_int_distribution<uint32_t> t_dist(0, 1);
    std::vector<uint32_t> t_values(n);
    for (size_t i = 0; i < n; ++i) t_values[i] = t_dist(eng);

    std::uniform_int_distribution<uint32_t> r_dist(1, -1);
    std::vector<uint32_t> mul_rands(34);
    for (int i = 0; i < 34; ++i) mul_rands[i] = r_dist(eng);

    std::vector<std::vector<uint32_t>> permutations(n, std::vector<uint32_t>(34));
    std::vector<uint32_t> base_perm(34);
    std::iota(base_perm.begin(), base_perm.end(), 0);
    for (size_t i = 0; i < n; ++i) {
        permutations[i] = base_perm;
        std::shuffle(permutations[i].begin(), permutations[i].end(), eng);
    }
    
    std::cout << "\n--- DEBUG P2 (Shared Randomness for i=" << debug_index << ") ---" << std::endl;
    std::cout << "    Generated t: " << t_values[debug_index] << std::endl;
    print_vector_full("Multiplicative randoms R", mul_rands);
    print_vector_full("Permutation P", permutations[debug_index]);
    
    // --- P0 and P1 prepare and send data ---
    auto masked_vectors_p0 = p0.prepareComparisonShares(x0, t_values, permutations, mul_rands, debug_index);
    auto masked_vectors_p1 = p1.prepareComparisonShares(x1, t_values, permutations, mul_rands, debug_index);

    // --- P2 reconstructs and checks for zero ---
    std::vector<uint32_t> r_prime(n);
    for (size_t i = 0; i < n; ++i) {
        bool is_debug = (i == debug_index);
        is_debug = 0;
        
        std::vector<uint32_t> reconstructed_shuffled_vector(34);
        for(int j=0; j<34; ++j) {
            reconstructed_shuffled_vector[j] = masked_vectors_p0[i][j] + masked_vectors_p1[i][j];
        }

        if(is_debug) {
            std::cout << "\n--- DEBUG P2 (Reconstruction for i=" << i << ") ---" << std::endl;
            print_vector_full("Reconstructed shuffled vector", reconstructed_shuffled_vector);
        }

        bool zero_found = false;
        for (const auto& val : reconstructed_shuffled_vector) {
            if (val == 0) {
                zero_found = true;
                break;
            }
        }
        r_prime[i] = zero_found ? 1 : 0;
    }
    
    std::cout << "\n--- DEBUG P2 ---" << std::endl;
    std::cout << "    Step 6: Intermediate result r' computed for all inputs: ";
    for(auto val : r_prime) std::cout << val << " ";
    std::cout << std::endl;

    // --- P2 shares and distributes r' ---
    std::uniform_int_distribution<uint32_t> share_dist;
    std::vector<uint32_t> r_prime_p0(n);
    std::vector<uint32_t> r_prime_p1(n);
    for (size_t i = 0; i < n; ++i) {
        r_prime_p0[i] = share_dist(eng);
        r_prime_p1[i] = r_prime[i] - r_prime_p0[i];
    }
    p0.receiveComparisonResultShare(r_prime_p0);
    p1.receiveComparisonResultShare(r_prime_p1);
    
    std::cout << "    Step 6: r' share for P0 (r'_0): ";
    for(auto val : r_prime_p0) std::cout << val << " ";
    std::cout << std::endl;
    std::cout << "    Step 6: r' share for P1 (r'_1): ";
    for(auto val : r_prime_p1) std::cout << val << " ";
    std::cout << std::endl;

    return t_values;
}

std::pair<std::vector<uint32_t>, std::vector<uint32_t>> P2::internalRunGTEZeroProtocol(
    const std::vector<uint32_t>& x0, P0& p0,
    const std::vector<uint32_t>& x1, P1& p1,
    int debug_index
) {
    // ... The entire implementation of the GTE-Zero protocol is here ...
    // ... from generating randomness (t, R, P) ...
    // ... to calling prepareComparisonShares ...
    // ... to reconstructing and checking for zero to get r_prime ...
    // ... to sharing r_prime and sending the shares to P0/P1 ...

    size_t n = x0.size();
    std::random_device rd;
    std::mt19937 eng(rd());

    // --- Step 1 & 5: Generate all shared randomness ---
    std::uniform_int_distribution<uint32_t> t_dist(0, 1);
    std::vector<uint32_t> t_values(n);
    for (size_t i = 0; i < n; ++i) t_values[i] = t_dist(eng);

    std::uniform_int_distribution<uint32_t> r_dist(1, -1);
    std::vector<uint32_t> mul_rands(34);
    for (int i = 0; i < 34; ++i) mul_rands[i] = r_dist(eng);

    std::vector<std::vector<uint32_t>> permutations(n, std::vector<uint32_t>(34));
    std::vector<uint32_t> base_perm(34);
    std::iota(base_perm.begin(), base_perm.end(), 0);
    for (size_t i = 0; i < n; ++i) {
        permutations[i] = base_perm;
        std::shuffle(permutations[i].begin(), permutations[i].end(), eng);
    }
    
    /*
    std::cout << "\n--- DEBUG P2 (Shared Randomness for i=" << debug_index << ") ---" << std::endl;
    std::cout << "    Generated t: " << t_values[debug_index] << std::endl;
    print_vector_full("Multiplicative randoms R", mul_rands);
    print_vector_full("Permutation P", permutations[debug_index]);
    */

    // --- P0 and P1 prepare and send data ---
    auto masked_vectors_p0 = p0.prepareComparisonShares(x0, t_values, permutations, mul_rands, debug_index);
    auto masked_vectors_p1 = p1.prepareComparisonShares(x1, t_values, permutations, mul_rands, debug_index);

    // --- P2 reconstructs and checks for zero ---
    std::vector<uint32_t> r_prime(n);
    for (size_t i = 0; i < n; ++i) {
        bool is_debug = (i == debug_index);
        
        std::vector<uint32_t> reconstructed_shuffled_vector(34);
        for(int j=0; j<34; ++j) {
            reconstructed_shuffled_vector[j] = masked_vectors_p0[i][j] + masked_vectors_p1[i][j];
        }

        /*
        if(is_debug) {
            std::cout << "\n--- DEBUG P2 (Reconstruction for i=" << i << ") ---" << std::endl;
            print_vector_full("Reconstructed shuffled vector", reconstructed_shuffled_vector);
        }
        */

        bool zero_found = false;
        for (const auto& val : reconstructed_shuffled_vector) {
            if (val == 0) {
                zero_found = true;
                break;
            }
        }
        r_prime[i] = zero_found ? 1 : 0;
    }
    
    //std::cout << "\n--- DEBUG P2 ---" << std::endl;
    //std::cout << "    Step 6: Intermediate result r' computed for all inputs: ";
    //for(auto val : r_prime) std::cout << val << " ";
    //std::cout << std::endl;

    // --- P2 shares and distributes r' ---
    std::uniform_int_distribution<uint32_t> share_dist;
    std::vector<uint32_t> r_prime_p0(n);
    std::vector<uint32_t> r_prime_p1(n);
    for (size_t i = 0; i < n; ++i) {
        r_prime_p0[i] = share_dist(eng);
        r_prime_p1[i] = r_prime[i] - r_prime_p0[i];
    }
    p0.receiveComparisonResultShare(r_prime_p0);
    p1.receiveComparisonResultShare(r_prime_p1);
    
    //std::cout << "    Step 6: r' share for P0 (r'_0): ";
    //for(auto val : r_prime_p0) std::cout << val << " ";
    //std::cout << std::endl;
    //std::cout << "    Step 6: r' share for P1 (r'_1): ";
    //for(auto val : r_prime_p1) std::cout << val << " ";
    //std::cout << std::endl;

    // The ONLY CHANGE is the return value:
    return {r_prime, t_values}; // Return both r' and t
}

void P2::performSecureEquality(
    const std::vector<uint32_t>& x0, const std::vector<uint32_t>& y0, P0& p0,
    const std::vector<uint32_t>& x1, const std::vector<uint32_t>& y1, P1& p1
) {
    size_t n = x0.size();

    // --- Step 1 & 2: P0/P1 compute shares locally ---
    // In our simulation, we compute them here.
    std::vector<uint32_t> xy0(n), yx0(n);
    std::vector<uint32_t> xy1(n), yx1(n);
    for (size_t i = 0; i < n; ++i) {
        xy0[i] = x0[i] - y0[i];
        yx0[i] = y0[i] - x0[i];
        xy1[i] = x1[i] - y1[i];
        yx1[i] = y1[i] - x1[i];
    }

    // --- NEW: Combine [x-y] and [y-x] into a single batch of size 2n ---
    std::vector<uint32_t> batch0;
    batch0.reserve(2 * n);
    batch0.insert(batch0.end(), xy0.begin(), xy0.end());
    batch0.insert(batch0.end(), yx0.begin(), yx0.end());

    std::vector<uint32_t> batch1;
    batch1.reserve(2 * n);
    batch1.insert(batch1.end(), xy1.begin(), xy1.end());
    batch1.insert(batch1.end(), yx1.begin(), yx1.end());


    // --- Step 3: Run GTE-Zero protocol ONCE on the combined batch ---
    std::cout << "P2: Running GTE-Zero on the combined batch..." << std::endl;
    auto combined_results = internalRunGTEZeroProtocol(batch0, p0, batch1, p1);
    auto all_r_prime = combined_results.first;  // This has 2n elements
    auto all_t_values = combined_results.second; // This has 2n elements
    std::cout << "P2: GTE-Zero protocol complete." << std::endl;


    // --- NEW: Split the 2n-sized results back into n-sized parts ---
    std::vector<uint32_t> r_prime_xy(all_r_prime.begin(), all_r_prime.begin() + n);
    std::vector<uint32_t> r_prime_yx(all_r_prime.begin() + n, all_r_prime.end());

    std::vector<uint32_t> t0(all_t_values.begin(), all_t_values.begin() + n);
    std::vector<uint32_t> t1(all_t_values.begin() + n, all_t_values.end());


    // --- Step 4: P2 computes XOR ---
    std::vector<uint32_t> xor_result(n);
    for (size_t i = 0; i < n; ++i) {
        xor_result[i] = r_prime_xy[i] ^ r_prime_yx[i];
    }
    std::cout << "P2: XOR of GTE-Zero results computed." << std::endl;

    // --- Step 5: P2 secret shares the XOR result ---
    std::vector<uint32_t> d0(n), d1(n);
    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_int_distribution<uint32_t> dist;
    for (size_t i = 0; i < n; ++i) {
        d0[i] = dist(eng);
        d1[i] = xor_result[i] - d0[i];
    }
    p0.receiveEqualityResultShare(d0);
    p1.receiveEqualityResultShare(d1);
    std::cout << "P2: XOR result shared and sent to P0/P1." << std::endl;

    // --- Step 6 (Implicit): P2 sends t0 and t1 to P0 and P1 ---
    p0.receiveEqualityTValues(t0, t1);
    p1.receiveEqualityTValues(t0, t1);
    std::cout << "P2: t0 and t1 sent to P0/P1 for final computation." << std::endl;
}


