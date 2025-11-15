# PPMICV
# Privacy-Preserving Medical Insurance Claim Verification

This repository contains the official implementation of the scheme described in the paper **"Privacy-Preserving Medical Insurance Claim Verification Scheme"**. The project utilizes secure multi-party computation (MPC) to verify medical insurance claims without compromising patient privacy.

The code implements a three-party protocol involving an insurance hospital (P0), a company (P1), and a regulatory authority (P2) to securely compute the claim payout amount.

## Overview

The core challenge in medical insurance claim verification is the conflict between the insurer's need to validate payout amounts and the paramount importance of patient privacy. Current practices often require hospitals to disclose complete, sensitive medical records, leading to significant privacy risks.

This project addresses this problem by leveraging a novel MPC-based scheme. The hospital provides encrypted medical records, the insurer provides the policy rules, and a third party facilitates the computation. The final output is the precise payout amount—revealed only to the insurer—without ever decrypting the underlying medical records.

## Features

- **Secure MPC Protocols:** Full implementation of the secure multi-party computation protocols described in the paper:
    - Secret Sharing
    - Secure Multiplication 
    - Secure Comparison 
- **Three-Party Claim Verification:** A complete simulation of the claim verification workflow between:
    - **P0:** The hospital, which inputs the patient's medical records.
    - **P1:** The insurance company, which inputs the policy rules.
    - **P2:** A trusted third party that coordinates the computation.
- **Core Logic:** The `main.cpp` file contains the primary logic that orchestrates the interaction between the three parties, following `Algorithm 5` from the paper.

## Data Generation

The medical data used in this project is synthetic data generated using [Synthea™](https://github.com/synthetichealth/synthea), an open-source patient population simulator. This allows us to test the scheme with realistic data without using real protected health information (PHI).

## Prerequisites

To build and run this project, you will need the following software installed on your system:

- **Operating System:** Ubuntu 22.04 LTS
- **Compiler:** g++ (with support for C++11 or newer)
- **Build System:** CMake 3.16 or newer

You can install these dependencies on Ubuntu with the following command:
```bash
sudo apt-get update && sudo apt-get install build-essential cmake
```
## Installation and Usage

Follow these steps to compile and run the project from your terminal:

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/encf/PPMICV/your-repo-name.git
    cd your-repo-name
    ```

2.  **Configure the project with CMake:**
    ```bash
    cmake CMakeLists.txt
    ```

3.  **Compile the source code:**
    ```bash
    make
    ```

4.  **Run the main executable:**
    ```bash
    ./MAIN
    ```
    The program will simulate the MPC protocol and print the results to the console.



## License

This project is licensed under the [MIT License](LICENSE.md). Please see the LICENSE file for details.
