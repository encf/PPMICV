#ifndef P2_HPP
#define P2_HPP

#include <vector>
#include <cstdint>

class P0;
class P1;

class P2 {
public:
    void generateAndDistributeTriples(int n, P0& p0, P1& p1);

    const std::vector<uint32_t>& getOriginalA() const;
    const std::vector<uint32_t>& getOriginalB() const;
    const std::vector<uint32_t>& getOriginalC() const;


    std::vector<uint32_t> performSecureComparison(
        const std::vector<uint32_t>& x0, P0& p0,
        const std::vector<uint32_t>& x1, P1& p1,
        int debug_index
    );

    void performSecureEquality(
        const std::vector<uint32_t>& x0, const std::vector<uint32_t>& y0, P0& p0,
        const std::vector<uint32_t>& x1, const std::vector<uint32_t>& y1, P1& p1
    );

    std::pair<std::vector<uint32_t>, std::vector<uint32_t>> internalRunGTEZeroProtocol(
        const std::vector<uint32_t>& z0, P0& p0,
        const std::vector<uint32_t>& z1, P1& p1,
        int debug_index = 0 // Keep debug functionality
    );

private:
    std::vector<uint32_t> a_;
    std::vector<uint32_t> b_;
    std::vector<uint32_t> c_;

};

#endif // P2_HPP
