#ifndef P0_H
#define P0_H

#include <string>
#include <vector>
#include <cstdint>
#include <utility>
#include "P1.hpp"

class P1;
class P0 {
public:
    bool readAndProcessCsv(const std::string& filename);
    const std::vector<std::vector<uint32_t>>& getDataColumns() const;
    void printData() const;

    std::vector<std::vector<uint32_t>> createSecretShares();
    const std::vector<std::vector<uint32_t>>& getRetainedShareColumns() const;

    void setArrays(const std::vector<uint32_t>& a, 
                   const std::vector<uint32_t>& b, 
                   const std::vector<uint32_t>& c);
    const std::vector<uint32_t>& getArrayA() const;
    const std::vector<uint32_t>& getArrayB() const;
    const std::vector<uint32_t>& getArrayC() const;


    void receiveP1Share(const P1RandomShare& p1_share);
    const P1RandomShare& getReceivedP1Share() const;


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
        const std::vector<uint32_t>& x0,
        const std::vector<uint32_t>& y0
    );

    std::vector<uint32_t> compute_z0_share(
        const std::vector<uint32_t>& d_plain,
        const std::vector<uint32_t>& e_plain
    );

    uint32_t computeDotProductShare(
        const std::vector<uint32_t>& d_plain,
        const std::vector<uint32_t>& e_plain
    );
    
    void readAndShareMedicalData(const std::string& filename, P1& p1);
    void receivePolicyDataShares(const std::vector<uint32_t>& policy_shares);

    void expandPolicyData(size_t n);

    std::vector<uint32_t> type_m_share_0_, code_m_share_0_, date_m_share_0_, cost_m_share_0_;
    std::vector<uint32_t> type_share_0_, code_share_0_, start_date_share_0_, end_date_share_0_;

    uint32_t max_cost_share_0_, ratio_share_0_;
    std::vector<uint32_t> in_coverage_share_0_, in_codes_share_0_, date_ok1_share_0_, date_ok2_share_0_;
    std::vector<uint32_t> m1_share_0_, m2_share_0_, m_share_0_;
    uint32_t total_cost_share_0_, base_amount_share_0_, meets_threshold_share_0_, r1_share_0_, r2_share_0_;


private:
    std::vector<std::vector<uint32_t>> dataColumns_;
    std::vector<std::vector<uint32_t>> retainedColumnShares_;

    std::vector<uint32_t> array_a_;
    std::vector<uint32_t> array_b_;
    std::vector<uint32_t> array_c_;


    P1RandomShare received_p1_share_;

    std::vector<uint32_t> r_prime_share_p0_;

    std::vector<uint32_t> t0_share_;
    std::vector<uint32_t> t1_share_;
    std::vector<uint32_t> d0_final_share_;



};
#endif
