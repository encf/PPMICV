#ifndef P1_HPP
#define P1_HPP

#include <cstdint>
#include <vector>
#include <string>
#include <utility>

// A structure to hold the random share, making the return type clear.
struct P1RandomShare {
    uint32_t type;
    uint32_t code;
    uint32_t start_date;
    uint32_t end_date;
    uint32_t max_cost;
};

class P0;

class P1 {
public:
    // Constructor
    P1();

    // --- Original Getters for initial values ---
    uint32_t getType() const;
    uint32_t getCode() const;
    uint32_t getStartDate() const;
    uint32_t getEndDate() const;
    uint32_t getMaxCost() const;

    // --- New Secret Sharing Functionality ---
    P1RandomShare createSecretShares();

    // --- New Getters for the retained share ---
    uint32_t getRetainedType() const;
    uint32_t getRetainedCode() const;
    uint32_t getRetainedStartDate() const;
    uint32_t getRetainedEndDate() const;
    uint32_t getRetainedMaxCost() const;

    void receiveP0Share(const std::vector<std::vector<uint32_t>>& p0_share);
    const std::vector<std::vector<uint32_t>>& getReceivedP0Share() const;


    void receiveTriplesShare(const std::vector<uint32_t>& a_share,
                             const std::vector<uint32_t>& b_share,
                             const std::vector<uint32_t>& c_share);

    const std::vector<uint32_t>& getTriplesShareA() const;
    const std::vector<uint32_t>& getTriplesShareB() const;
    const std::vector<uint32_t>& getTriplesShareC() const;

    std::vector<std::vector<uint32_t>> prepareComparisonShares(
        const std::vector<uint32_t>& x_share,
        const std::vector<uint32_t>& t_values,
        const std::vector<std::vector<uint32_t>>& permutations,
        const std::vector<uint32_t>& mul_rands,
        int debug_index
    );

    void receiveComparisonResultShare(const std::vector<uint32_t>& r_prime_share);
    const std::vector<uint32_t>& getComparisonResultShare() const;
    std::vector<uint32_t> computeFinalComparisonResult(const std::vector<uint32_t>& t_values);

    void receiveEqualityTValues(const std::vector<uint32_t>& t0, const std::vector<uint32_t>& t1);
    void receiveEqualityResultShare(const std::vector<uint32_t>& d0_share);
    std::vector<uint32_t> computeFinalEqualityResult();

    std::pair<std::vector<uint32_t>, std::vector<uint32_t>> compute_d_e_shares(
        const std::vector<uint32_t>& x1,
        const std::vector<uint32_t>& y1
    );

    std::vector<uint32_t> compute_z1_share(
        const std::vector<uint32_t>& d_plain,
        const std::vector<uint32_t>& e_plain
    );

    uint32_t computeDotProductShare(
        const std::vector<uint32_t>& d_plain,
        const std::vector<uint32_t>& e_plain
    );


    void sharePolicyData(P0& p0);

    void receiveMedicalDataShares(const std::vector<std::vector<uint32_t>>& medical_shares);

    void expandPolicyData(size_t n);

    uint32_t type_, code_, start_date_, end_date_, max_cost_, ratio_;
    std::vector<uint32_t> type_m_share_1_, code_m_share_1_, date_m_share_1_, cost_m_share_1_;
    std::vector<uint32_t> type_share_1_, code_share_1_, start_date_share_1_, end_date_share_1_;
    uint32_t max_cost_share_1_;
    uint32_t ratio_share_1_;
    std::vector<uint32_t> in_coverage_share_1_;
    std::vector<uint32_t> in_codes_share_1_;
    std::vector<uint32_t> date_ok1_share_1_;
    std::vector<uint32_t> date_ok2_share_1_;
    std::vector<uint32_t> m1_share_1_;
    std::vector<uint32_t> m2_share_1_;
    std::vector<uint32_t> m_share_1_;
    uint32_t total_cost_share_1_;
    uint32_t base_amount_share_1_;
    uint32_t meets_threshold_share_1_;
    uint32_t r1_share_1_;
    uint32_t r2_share_1_;

private:
    // Original values

    // P1's retained share (Original - Random)
    uint32_t retained_type_;
    uint32_t retained_code_;
    uint32_t retained_start_date_;
    uint32_t retained_end_date_;
    uint32_t retained_max_cost_;


    std::vector<std::vector<uint32_t>> received_p0_share_;

    std::vector<uint32_t> triples_a_share_;
    std::vector<uint32_t> triples_b_share_;
    std::vector<uint32_t> triples_c_share_;

    std::vector<uint32_t> r_prime_share_p1_;

    std::vector<uint32_t> t0_share_;
    std::vector<uint32_t> t1_share_;
    std::vector<uint32_t> d1_final_share_;



};

#endif // P1_HPP
