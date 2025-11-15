#include "P0.hpp"
#include "P1.hpp"
#include "P2.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <cstdint>
#include <random>
#include <iomanip>
#include <numeric>
#include <fstream>
#include <sstream>

void print_debug_scalar_shares(
    const std::string& name,
    uint32_t p0_share,
    uint32_t p1_share
) {
    std::cout << "    DEBUG - " << std::left << std::setw(20) << name
              << " | P0 Share: " << std::setw(12) << p0_share
              << " | P1 Share: " << std::setw(12) << p1_share
              << " | Reconstructed: " << p0_share + p1_share << std::endl;
}

// Helper to print and reconstruct shares of a VECTOR value (first/last 3 elements)
void print_debug_vector_shares(
    const std::string& name,
    const std::vector<uint32_t>& p0_shares,
    const std::vector<uint32_t>& p1_shares,
    int n_print = 10
) {
    std::cout << "    DEBUG - " << name << " (showing first/last " << n_print << " elements):" << std::endl;
    size_t n = p0_shares.size();
    for (int i = 0; i < n; ++i) {
        if (i < n_print || i >= n - n_print) {
            std::cout << "      [" << std::setw(3) << i << "] "
                      << " | P0 Share: " << std::setw(12) << p0_shares[i]
                      << " | P1 Share: " << std::setw(12) << p1_shares[i]
                      << " | Reconstructed: " << p0_shares[i] + p1_shares[i] << std::endl;
        }
        if (i == n_print && n > n_print * 2) {
            std::cout << "      ..." << std::endl;
        }
    }
}

void printColumnData(const std::string& title, const std::vector<std::vector<uint32_t>>& columns, int n_rows) {
    std::cout << title << std::endl;
    if (columns.empty() || columns[0].empty()) {
        std::cout << "data is empty"<< std::endl;
        std::cout << "----------------------------------------------------" << std::endl;
        return;
    }
    size_t num_rows = columns[0].size();
    size_t num_cols = columns.size();
    for (size_t i = 0; i < n_rows && i < num_rows; ++i) {
        for (size_t j = 0; j < num_cols; ++j) {
            std::cout << columns[j][i] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "----------------------------------------------------" << std::endl;
}

void printP1ShareData(const std::string& title, const P1RandomShare& share) {
    std::cout << title << std::endl;
    std::cout << "  Type: " << share.type << std::endl;
    std::cout << "  Code: " << share.code << std::endl;
    std::cout << "  Start Date: " << share.start_date << std::endl;
    std::cout << "  End Date: " << share.end_date << std::endl;
    std::cout << "  Max Cost: " << share.max_cost << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;
}

bool verifyP0Sharing(const std::vector<std::vector<uint32_t>>& original,
                     const std::vector<std::vector<uint32_t>>& share1,
                     const std::vector<std::vector<uint32_t>>& share2) {
    if (original.empty()) return true;
    if (original.size() != share1.size() || original.size() != share2.size()) return false;
    size_t num_cols = original.size();
    size_t num_rows = original[0].size();
    for (size_t j = 0; j < num_cols; ++j) {
        if (original[j].size() != num_rows || share1[j].size() != num_rows || share2[j].size() != num_rows) return false;
        for (size_t i = 0; i < num_rows; ++i) {
            if (share1[j][i] + share2[j][i] != original[j][i]) return false;
        }
    }
    return true;
}

bool verifyP1Sharing(const P1& original_p1, const P1RandomShare& random_share) {
    bool ok = true;
    #define VERIFY_FIELD(field, getter, retained_getter) \
    do { \
        uint32_t reconstructed = original_p1.retained_getter() + random_share.field; \
        std::cout << "  " #field ": " << original_p1.retained_getter() << " + " << random_share.field << " = " << reconstructed; \
        if (reconstructed == original_p1.getter()) { \
            std::cout << " (== " << original_p1.getter() << ", OK)" << std::endl; \
        } else { \
            std::cout << " (!= " << original_p1.getter() << ", FAILED)" << std::endl; \
            ok = false; \
        } \
    } while(0)
    
    VERIFY_FIELD(type, getType, getRetainedType);
    VERIFY_FIELD(code, getCode, getRetainedCode);
    VERIFY_FIELD(start_date, getStartDate, getRetainedStartDate);
    VERIFY_FIELD(end_date, getEndDate, getRetainedEndDate);
    VERIFY_FIELD(max_cost, getMaxCost, getRetainedMaxCost);
    
    return ok;
}

void printTripleRow(const std::string& prefix, int index,
                    const std::vector<uint32_t>& a,
                    const std::vector<uint32_t>& b,
                    const std::vector<uint32_t>& c) {
    if (a.size() > index && b.size() > index && c.size() > index) {
        std::cout << prefix << "[" << index << "]:\t"
                  << a[index] << "\t" << b[index] << "\t" << c[index] << std::endl;
    }
}

int secret_and_generate_triples()
{

    P0 p0;
    P1 p1;
    P2 p2;


    if (!p0.readAndProcessCsv("data.csv")) {
        std::cerr << "error: read csv data failed" << std::endl;
        return 1;
    }
    auto p0_random_share = p0.createSecretShares();
    p1.receiveP0Share(p0_random_share);

    auto p1_random_share = p1.createSecretShares();
    p0.receiveP1Share(p1_random_share);



    P1RandomShare p1_original_data = { p1.getType(), p1.getCode(), p1.getStartDate(), p1.getEndDate(), p1.getMaxCost() };
    
    P1RandomShare p1_retained_share = { p1.getRetainedType(), p1.getRetainedCode(), p1.getRetainedStartDate(), p1.getRetainedEndDate(), p1.getRetainedMaxCost() };



    bool p0_verified = verifyP0Sharing(p0.getDataColumns(), p0.getRetainedShareColumns(), p1.getReceivedP0Share());
    
    std::cout << std::endl;

    bool p1_verified = verifyP1Sharing(p1, p0.getReceivedP1Share());

    
     const int NUM_TRIPLES = 5;

    p2.generateAndDistributeTriples(NUM_TRIPLES, p0, p1);
    
    int index_to_print = 0;


    const auto& a0 = p0.getArrayA();
    const auto& b0 = p0.getArrayB();
    const auto& c0 = p0.getArrayC();
    const auto& a1 = p1.getTriplesShareA();
    const auto& b1 = p1.getTriplesShareB();
    const auto& c1 = p1.getTriplesShareC();

    bool all_verified = true;
    for (int i = 0; i < NUM_TRIPLES; ++i) {
        uint32_t recon_a = a0[i] + a1[i];
        uint32_t recon_b = b0[i] + b1[i];
        uint32_t recon_c = c0[i] + c1[i];


        if (recon_a * recon_b == recon_c) {
        } else {
            all_verified = false;
        }
    }
    
}

void getz()
{
    P0 p0;
    P1 p1;
    P2 p2;

    // --- SETUP ---
    const int N_VALUES = 2; // Keep it small for debugging
    const int DEBUG_INDEX = 0; // The index we want to trace
    std::vector<int32_t> original_x_signed = {-3, 5, 0, 9, 123, -21344}; // One negative, one positive

    std::cout << "====================== DEBUG SETUP ======================" << std::endl;
    std::cout << "Tracing index: " << DEBUG_INDEX << " (value: " << original_x_signed[DEBUG_INDEX] << ")" << std::endl;

    size_t n = original_x_signed.size();
    std::vector<uint32_t> original_x(n);
    for(size_t i=0; i<n; ++i) original_x[i] = static_cast<uint32_t>(original_x_signed[i]);

    std::vector<uint32_t> x0(n), x1(n);
    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_int_distribution<uint32_t> dist;
    for (size_t i = 0; i < n; ++i) {
        x0[i] = dist(eng);
        x1[i] = original_x[i] - x0[i];
    }
    uint32_t recon_x = x0[DEBUG_INDEX] + x1[DEBUG_INDEX];
    std::cout << "Original x[" << DEBUG_INDEX << "] (" << original_x[DEBUG_INDEX] << ") -> x0="
              << x0[DEBUG_INDEX] << ", x1=" << x1[DEBUG_INDEX] << " | Recon: " << recon_x << std::endl;


    // --- EXECUTE PROTOCOL WITH DEBUGGING ---
    std::cout << "\n================== PROTOCOL EXECUTION TRACE ==================" << std::endl;
    auto t_values = p2.performSecureComparison(x0, p0, x1, p1, DEBUG_INDEX);


    // --- FINAL COMPUTATION & VERIFICATION ---
    std::cout << "\n================== FINAL COMPUTATION & VERIFICATION ==================" << std::endl;
    auto final_r0 = p0.computeFinalComparisonResult(t_values);
    auto final_r1 = p1.computeFinalComparisonResult(t_values);

    for(size_t i = 0; i < n; i++) {
        bool expected = (original_x_signed[i] >= 0);
        uint32_t final_r = final_r0[i] + final_r1[i];

        std::cout << "\n--- Verifying Index " << i << " (Original: " << original_x_signed[i] << ") ---" << std::endl;
        std::cout << "  P2 sent r'_0=" << p0.getComparisonResultShare()[i] << ", r'_1=" << p1.getComparisonResultShare()[i] << std::endl;
        std::cout << "  t was: " << t_values[i] << std::endl;
        std::cout << "  P0 computed final_r0 = t + r'_0 - 2*t*r'_0 = " << final_r0[i] << std::endl;
        std::cout << "  P1 computed final_r1 =     r'_1 - 2*t*r'_1 = " << final_r1[i] << std::endl;
        std::cout << "  Reconstructed r = final_r0 + final_r1 = " << final_r << std::endl;
        std::cout << "  Expected result (>=0): " << (expected ? "1" : "0") << std::endl;
        std::cout << "  Final Status: " << ((final_r == expected) ? "PASS" : "FAIL") << std::endl;
    }
}

void eq()
{
    P0 p0;
    P1 p1;
    P2 p2;

    std::vector<uint32_t> original_x = {100, 200, 50, 400, 999};
    std::vector<uint32_t> original_y = {100, 300, 50, 401, 998};
    size_t n = original_x.size();

    std::vector<uint32_t> x0(n), y0(n), x1(n), y1(n);
    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_int_distribution<uint32_t> dist;
    for (size_t i = 0; i < n; ++i) {
        x0[i] = dist(eng);
        y0[i] = dist(eng);
        x1[i] = original_x[i] - x0[i];
        y1[i] = original_y[i] - y0[i];
    }

    for(size_t i=0; i<n; ++i) std::cout << "  x=" << original_x[i] << ", y=" << original_y[i] << std::endl;

    p2.performSecureEquality(x0, y0, p0, x1, y1, p1);

    auto r0 = p0.computeFinalEqualityResult();
    auto r1 = p1.computeFinalEqualityResult();


    for(size_t i = 0; i < n; i++) {
        uint32_t final_r = r0[i] + r1[i];
        bool expected_result = (original_x[i] == original_y[i]);

        std::cout << std::left << std::setw(10) << original_x[i]
                  << std::setw(10) << original_y[i]
                  << std::setw(18) << (expected_result ? "true (1)" : "false (0)")
                  << std::setw(15) << final_r;

        if (final_r == expected_result) {
            std::cout << "PASS" << std::endl;
        } else {
            std::cout << "FAIL" << std::endl;
        }
    }
}

void mul()
{
    P0 p0;
    P1 p1;
    P2 p2;

    std::vector<uint32_t> original_x = {10, 20, 30, 40};
    std::vector<uint32_t> original_y = {5, 6, 7, 8};
    size_t n = original_x.size();

    // Secret share x and y
    std::vector<uint32_t> x0(n), y0(n), x1(n), y1(n);
    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_int_distribution<uint32_t> dist;
    for (size_t i = 0; i < n; ++i) {
        x0[i] = dist(eng);
        y0[i] = dist(eng);
        x1[i] = original_x[i] - x0[i];
        y1[i] = original_y[i] - y0[i];
    }

    // Generate and distribute multiplication triples from P2
    std::cout << "P2 is generating and distributing multiplication triples..." << std::endl;
    p2.generateAndDistributeTriples(n, p0, p1);
    std::cout << "Triples have been distributed to P0 and P1.\n" << std::endl;

    std::cout << "--- Starting Secure Multiplication Protocol ---" << std::endl;

    // Step 1: P0 and P1 compute their d, e shares
    std::cout << "Step 1: P0 and P1 are computing their d, e shares locally..." << std::endl;
    auto p0_d_e_shares = p0.compute_d_e_shares(x0, y0);
    auto p1_d_e_shares = p1.compute_d_e_shares(x1, y1);

    // Step 2: Main (as the network) reconstructs d and e
    std::cout << "Step 2: Reconstructing plaintext d and e..." << std::endl;
    std::vector<uint32_t> d_plain(n), e_plain(n);
    for (size_t i = 0; i < n; ++i) {
        d_plain[i] = p0_d_e_shares.first[i] + p1_d_e_shares.first[i];
        e_plain[i] = p0_d_e_shares.second[i] + p1_d_e_shares.second[i];
    }
    std::cout << "  Plaintext d[0]: " << d_plain[0] << std::endl;
    std::cout << "  Plaintext e[0]: " << e_plain[0] << std::endl;

    // Step 3 & 4: P0 and P1 compute their final z shares
    std::cout << "Step 3 & 4: P0 and P1 are computing their final z shares..." << std::endl;
    auto z0 = p0.compute_z0_share(d_plain, e_plain);
    auto z1 = p1.compute_z1_share(d_plain, e_plain);
    std::cout << "Protocol finished.\n" << std::endl;


    for(size_t i = 0; i < n; i++) {
        uint32_t reconstructed_z = z0[i] + z1[i];
        uint32_t expected_result = original_x[i] * original_y[i];

        std::cout << std::left << std::setw(10) << original_x[i]
                  << std::setw(10) << original_y[i]
                  << std::setw(18) << expected_result
                  << std::setw(15) << reconstructed_z;

        if (reconstructed_z == expected_result) {
            std::cout << "PASS" << std::endl;
        } else {
            std::cout << "FAIL" << std::endl;
        }
    }
}

void inner_product()
{
    P0 p0;
    P1 p1;
    P2 p2;

    std::vector<uint32_t> original_x = {10, 20, 30, 40};
    std::vector<uint32_t> original_y = {5, 6, 7, 8};
    size_t n = original_x.size();

    std::vector<uint32_t> x0(n), y0(n), x1(n), y1(n);
    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_int_distribution<uint32_t> dist;
    for (size_t i = 0; i < n; ++i) {
        x0[i] = dist(eng);
        y0[i] = dist(eng);
        x1[i] = original_x[i] - x0[i];
        y1[i] = original_y[i] - y0[i];
    }

    // Generate and distribute multiplication triples from P2
    std::cout << "P2 is generating and distributing multiplication triples..." << std::endl;
    p2.generateAndDistributeTriples(n, p0, p1);
    std::cout << "Triples have been distributed to P0 and P1.\n" << std::endl;

    // --- 执行协议 ---
    std::cout << "--- Starting Secure Dot Product Protocol ---" << std::endl;

    // Step 1: P0 and P1 compute their d, e shares
    std::cout << "Step 1: P0 and P1 are computing their d, e shares locally..." << std::endl;
    auto p0_d_e_shares = p0.compute_d_e_shares(x0, y0);
    auto p1_d_e_shares = p1.compute_d_e_shares(x1, y1);

    // Step 2: Main (as the network) reconstructs d and e
    std::cout << "Step 2: Reconstructing plaintext d and e..." << std::endl;
    std::vector<uint32_t> d_plain(n), e_plain(n);
    for (size_t i = 0; i < n; ++i) {
        d_plain[i] = p0_d_e_shares.first[i] + p1_d_e_shares.first[i];
        e_plain[i] = p0_d_e_shares.second[i] + p1_d_e_shares.second[i];
    }

    // Step 3, 4, 5: P0 and P1 compute their final dot product shares
    std::cout << "Step 3-5: P0 and P1 are computing their final dot product shares..." << std::endl;
    uint32_t dot_product_share0 = p0.computeDotProductShare(d_plain, e_plain);
    uint32_t dot_product_share1 = p1.computeDotProductShare(d_plain, e_plain);
    std::cout << "Protocol finished.\n" << std::endl;



    // Calculate expected result using std::inner_product for elegance
    uint32_t expected_result = std::inner_product(original_x.begin(), original_x.end(), original_y.begin(), (uint32_t)0);

    // Reconstruct the result from shares
    uint32_t reconstructed_result = dot_product_share0 + dot_product_share1;

    std::cout << "Original x: [ ";
    for(const auto& v : original_x) std::cout << v << " ";
    std::cout << "]" << std::endl;

    std::cout << "Original y: [ ";
    for(const auto& v : original_y) std::cout << v << " ";
    std::cout << "]" << std::endl;

    std::cout << "\nExpected Dot Product: " << expected_result << std::endl;
    std::cout << "Reconstructed Dot Product: " << reconstructed_result << std::endl;

    if (reconstructed_result == expected_result) {
        std::cout << "\nStatus: PASS" << std::endl;
    } else {
        std::cout << "\nStatus: FAIL" << std::endl;
    }
}


double run_plaintext_verification(const std::string& filename) {
    std::cout << "--- Running Plaintext Calculation for Verification ---" << std::endl;
    uint32_t type_ = 1, code_ = 410620009, start_date_ = 20240101, end_date_ = 20251231;
    double max_cost_ = 10000, ratio_ = 1.0;

    std::vector<uint32_t> type_m, code_m, date_m;
    std::vector<double> cost_m;
    std::ifstream file(filename);
    std::string line;
    std::getline(file, line); // Skip header

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> rowData;
        while (std::getline(ss, cell, ',')) rowData.push_back(cell);
        if (rowData.size() == 5) {
            type_m.push_back(stoul(rowData[1]));
            code_m.push_back(stoul(rowData[2]));
            date_m.push_back(stoul(rowData[3]));
            cost_m.push_back(stod(rowData[4]));
        }
    }
    file.close();
    size_t n = type_m.size();

    std::vector<uint32_t> m(n);
    for(size_t i=0; i<n; ++i) {
        bool in_coverage = (type_m[i] == type_);
        bool in_codes = (code_m[i] == code_);
        bool date_ok1 = (date_m[i] >= start_date_);
        bool date_ok2 = (end_date_ >= date_m[i]);
        m[i] = in_coverage && in_codes && date_ok1 && date_ok2;
    }

    double total_cost = 0;
    for(size_t i=0; i<n; ++i) {
        total_cost += m[i] * cost_m[i];
    }

    double base_amount = total_cost * ratio_;
    bool meets_threshold = (total_cost >= max_cost_);

    uint32_t r1 = meets_threshold * max_cost_;
    double r2 = base_amount * (1 - meets_threshold);
    double r = r1 + r2;
    std::cout << "Plaintext total valid cost: " << total_cost << std::endl;
    std::cout << "Plaintext final result (r): " << r << "\n" << std::endl;
    return r;
}


void run_equality_and_store(
    const std::vector<uint32_t>& x0, const std::vector<uint32_t>& y0,
    const std::vector<uint32_t>& x1, const std::vector<uint32_t>& y1,
    P0& p0, P1& p1, P2& p2,
    std::vector<uint32_t>& result0, std::vector<uint32_t>& result1
) {
    p2.performSecureEquality(x0, y0, p0, x1, y1, p1);
    result0 = p0.computeFinalEqualityResult();
    result1 = p1.computeFinalEqualityResult();
}

void run_gte_zero_and_store(
    const std::vector<uint32_t>& z0, const std::vector<uint32_t>& z1,
    P0& p0, P1& p1, P2& p2,
    std::vector<uint32_t>& result0, std::vector<uint32_t>& result1
) {
    auto t_values = p2.internalRunGTEZeroProtocol(z0, p0, z1, p1);
    result0 = p0.computeFinalComparisonResult(t_values.second);
    result1 = p1.computeFinalComparisonResult(t_values.second);
}

void run_vector_mult_and_store(
    const std::vector<uint32_t>& x0, const std::vector<uint32_t>& y0,
    const std::vector<uint32_t>& x1, const std::vector<uint32_t>& y1,
    P0& p0, P1& p1, P2& p2,
    std::vector<uint32_t>& result0, std::vector<uint32_t>& result1
) {
    size_t n = x0.size();
    p2.generateAndDistributeTriples(n, p0, p1);
    auto p0_de = p0.compute_d_e_shares(x0, y0);
    auto p1_de = p1.compute_d_e_shares(x1, y1);
    std::vector<uint32_t> d(n), e(n);
    for(size_t i=0; i<n; ++i) {
        d[i] = p0_de.first[i] + p1_de.first[i];
        e[i] = p0_de.second[i] + p1_de.second[i];
    }
    result0 = p0.compute_z0_share(d, e);
    result1 = p1.compute_z1_share(d, e);
}

void run_scalar_mult_and_store(
    uint32_t x0, uint32_t y0, uint32_t x1, uint32_t y1,
    P0& p0, P1& p1, P2& p2,
    uint32_t& result0, uint32_t& result1
) {
    std::vector<uint32_t> res0, res1;
    run_vector_mult_and_store({x0}, {y0}, {x1}, {y1}, p0, p1, p2, res0, res1);
    result0 = res0[0];
    result1 = res1[0];
}



int main() {
    double  r_expected = run_plaintext_verification("data.csv");
    std::cout << "--- Starting Secure Multi-Party Computation ---" << std::endl;

    P0 p0;
    P1 p1;
    P2 p2;
    size_t n = 100;

    // Step 1 & 2
    std::cout << "[Step 1,2] Data loading and initial sharing..." << std::endl;
    p0.readAndShareMedicalData("data.csv", p1);
    p1.sharePolicyData(p0);
    print_debug_vector_shares("cost_m (initial share)", p0.cost_m_share_0_, p1.cost_m_share_1_);
    print_debug_scalar_shares("max_cost (initial share)", p0.max_cost_share_0_, p1.max_cost_share_1_);


    // Step 3
    std::cout << "[Step 3] Expanding policy data..." << std::endl;
    p0.expandPolicyData(n);
    p1.expandPolicyData(n);
    
    // Step 4
    std::cout << "[Step 4] Equality check: type_m == type..." << std::endl;
    run_equality_and_store(p0.type_m_share_0_, p0.type_share_0_, p1.type_m_share_1_, p1.type_share_1_, p0, p1, p2, p0.in_coverage_share_0_, p1.in_coverage_share_1_);
    run_equality_and_store(p0.type_m_share_0_, p0.type_share_0_, p1.type_m_share_1_, p1.type_share_1_, p0, p1, p2, p0.in_coverage_share_0_, p1.in_coverage_share_1_);
    print_debug_vector_shares("in_coverage", p0.in_coverage_share_0_, p1.in_coverage_share_1_);
    
    // Step 5
    std::cout << "[Step 5] Equality check: code_m == code..." << std::endl;
    run_equality_and_store(p0.code_m_share_0_, p0.code_share_0_, p1.code_m_share_1_, p1.code_share_1_, p0, p1, p2, p0.in_codes_share_0_, p1.in_codes_share_1_);
    print_debug_vector_shares("in_codes", p0.in_codes_share_0_, p1.in_codes_share_1_);


    // Step 6
    std::cout << "[Step 6] GTE check: date_m >= start_date..." << std::endl;
    std::vector<uint32_t> d1_0(n), d1_1(n);
    for(size_t i=0; i<n; ++i) {
        d1_0[i] = p0.date_m_share_0_[i] - p0.start_date_share_0_[i];
        d1_1[i] = p1.date_m_share_1_[i] - p1.start_date_share_1_[i];
    }
    run_gte_zero_and_store(d1_0, d1_1, p0, p1, p2, p0.date_ok1_share_0_, p1.date_ok1_share_1_);
    print_debug_vector_shares("date_ok1 (d1>=0)", p0.date_ok1_share_0_, p1.date_ok1_share_1_);

    // Step 7
    std::cout << "[Step 7] GTE check: end_date >= date_m..." << std::endl;
    std::vector<uint32_t> d2_0(n), d2_1(n);
    for(size_t i=0; i<n; ++i) {
        d2_0[i] = p0.end_date_share_0_[i] - p0.date_m_share_0_[i];
        d2_1[i] = p1.end_date_share_1_[i] - p1.date_m_share_1_[i];
    }
    run_gte_zero_and_store(d2_0, d2_1, p0, p1, p2, p0.date_ok2_share_0_, p1.date_ok2_share_1_);
    print_debug_vector_shares("date_ok2 (d2>=0)", p0.date_ok2_share_0_, p1.date_ok2_share_1_);
    print_debug_vector_shares("date_m (d2>=0)", p0.date_m_share_0_, p1.date_m_share_1_);
    print_debug_vector_shares("end_date (d2>=0)", p0.end_date_share_0_, p1.end_date_share_1_);
    print_debug_vector_shares("diff (d2>=0)", d2_0, d2_1);

    // Step 8
    std::cout << "[Step 8] Multiplication: m1 = in_coverage * in_codes..." << std::endl;
    run_vector_mult_and_store(p0.in_coverage_share_0_, p0.in_codes_share_0_, p1.in_coverage_share_1_, p1.in_codes_share_1_, p0, p1, p2, p0.m1_share_0_, p1.m1_share_1_);
    print_debug_vector_shares("m1 (d2>=0)", p0.m1_share_0_, p1.m1_share_1_);
    print_debug_vector_shares("in_coverage (d2>=0)", p0.in_coverage_share_0_, p1.in_coverage_share_1_);
    print_debug_vector_shares("in_codes (d2>=0)", p0.in_codes_share_0_, p1.in_codes_share_1_);

    // Step 9
    std::cout << "[Step 9] Multiplication: m2 = date_ok1 * date_ok2..." << std::endl;
    run_vector_mult_and_store(p0.date_ok1_share_0_, p0.date_ok2_share_0_, p1.date_ok1_share_1_, p1.date_ok2_share_1_, p0, p1, p2, p0.m2_share_0_, p1.m2_share_1_);
    print_debug_vector_shares("m2", p0.m2_share_0_, p1.m2_share_1_);
    print_debug_vector_shares("date_ok1 (d1>=0)", p0.date_ok1_share_0_, p1.date_ok1_share_1_);

    // Step 10
    std::cout << "[Step 10] Multiplication: m = m1 * m2..." << std::endl;
    run_vector_mult_and_store(p0.m1_share_0_, p0.m2_share_0_, p1.m1_share_1_, p1.m2_share_1_, p0, p1, p2, p0.m_share_0_, p1.m_share_1_);
    print_debug_vector_shares("m (final mask)", p0.m_share_0_, p1.m_share_1_);

    // Step 11
    std::cout << "[Step 11] Dot Product: total_cost = m . cost_m..." << std::endl;
    p2.generateAndDistributeTriples(n, p0, p1);
    auto p0_de_dp = p0.compute_d_e_shares(p0.m_share_0_, p0.cost_m_share_0_);
    auto p1_de_dp = p1.compute_d_e_shares(p1.m_share_1_, p1.cost_m_share_1_);
    std::vector<uint32_t> d_dp(n), e_dp(n);
    for(size_t i=0; i<n; ++i) {
        d_dp[i] = p0_de_dp.first[i] + p1_de_dp.first[i];
        e_dp[i] = p0_de_dp.second[i] + p1_de_dp.second[i];
    }
    p0.total_cost_share_0_ = p0.computeDotProductShare(d_dp, e_dp);
    p1.total_cost_share_1_ = p1.computeDotProductShare(d_dp, e_dp);
    print_debug_scalar_shares("total_cost", p0.total_cost_share_0_, p1.total_cost_share_1_);

    // Step 12
    std::cout << "[Step 12] Multiplication: base_amount = total_cost * ratio..." << std::endl;
    std::cout << "_________________________________________\n";
    std::cout << p0.total_cost_share_0_ << " " << p1.total_cost_share_1_<< "\n";
    std::cout << p0.total_cost_share_0_ + p1.total_cost_share_1_ << "\n";
    std::cout << "_________________________________________\n";
    std::cout << p0.ratio_share_0_ << " " << p1.ratio_share_1_ << "\n";
    std::cout << p0.ratio_share_0_ +  p1.ratio_share_1_ << "\n";
    std::cout << "_________________________________________\n";


    run_scalar_mult_and_store(p0.total_cost_share_0_, p0.ratio_share_0_, p1.total_cost_share_1_, p1.ratio_share_1_, p0, p1, p2, p0.base_amount_share_0_, p1.base_amount_share_1_);
    print_debug_scalar_shares("base_amount", p0.base_amount_share_0_, p1.base_amount_share_1_);

    // Step 13
    std::cout << "[Step 13] GTE check: total_cost >= max_cost..." << std::endl;
    uint32_t diff0 = p0.total_cost_share_0_ - p0.max_cost_share_0_;
    uint32_t diff1 = p1.total_cost_share_1_ - p1.max_cost_share_1_;
    std::vector<uint32_t> meets0, meets1;
    run_gte_zero_and_store({diff0}, {diff1}, p0, p1, p2, meets0, meets1);
    p0.meets_threshold_share_0_ = meets0[0];
    p1.meets_threshold_share_1_ = meets1[0];
    print_debug_scalar_shares("meets_threshold", p0.meets_threshold_share_0_, p1.meets_threshold_share_1_);

    // Step 14
    std::cout << "[Step 14] Multiplication: r1 = meets_threshold * max_cost..." << std::endl;
    run_scalar_mult_and_store(p0.meets_threshold_share_0_, p0.max_cost_share_0_, p1.meets_threshold_share_1_, p1.max_cost_share_1_, p0, p1, p2, p0.r1_share_0_, p1.r1_share_1_);
    print_debug_scalar_shares("r1", p0.r1_share_0_, p1.r1_share_1_);

    std::cout << "[Step 14] Multiplication: r2 = base_amount * (1 - meets_threshold)..." << std::endl;
    uint32_t one_minus_mt0 = 1 - p0.meets_threshold_share_0_;
    uint32_t one_minus_mt1 =   - p1.meets_threshold_share_1_;
    run_scalar_mult_and_store(p0.base_amount_share_0_, one_minus_mt0, p1.base_amount_share_1_, one_minus_mt1, p0, p1, p2, p0.r2_share_0_, p1.r2_share_1_);
     print_debug_scalar_shares("r2", p0.r2_share_0_, p1.r2_share_1_);

    // Step 15
    std::cout << "[Step 15] Final reconstruction..." << std::endl;
    uint32_t r0 = p0.r1_share_0_ + p0.r2_share_0_;
    uint32_t r1 = p1.r1_share_1_ + p1.r2_share_1_;
    uint32_t r_reconstructed = r0 + r1;
    double mpc_result = static_cast<double>(r_reconstructed);
    mpc_result = mpc_result / 1024;
    print_debug_scalar_shares("r (final)", r0, r1);

    // --- Final Verification ---
    std::cout << "\n====================== FINAL VERIFICATION ======================" << std::endl;
    std::cout << "Expected Plaintext Result (r): " << r_expected << std::endl;
    std::cout << "Reconstructed MPC Result (r):  " << mpc_result << std::endl;
    return 0;
}
