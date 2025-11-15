#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <vector>
#include <string>
#include <cstdint> // For uint32_t

// Helper function to print a vector in a readable format
inline void print_vector_full(const std::string& title, const std::vector<uint32_t>& vec) {
    std::cout << "    " << title << " (" << vec.size() << " elements):" << std::endl;
    std::cout << "      [ ";
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i] << (i == vec.size() - 1 ? "" : ", ");
        if ((i + 1) % 8 == 0 && i < vec.size() - 1) { // Newline every 8 elements
            std::cout << "\n        ";
        }
    }
    std::cout << " ]" << std::endl;
}

#endif // UTILS_HPP
