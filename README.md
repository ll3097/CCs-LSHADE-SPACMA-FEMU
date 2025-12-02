# CEPSTRUM-BASED CUSTOMIZED LSHADE-SPACMA HYBRID FEMU FRAMEWORK (MATLAB Implementation)

**This repository provides the MATLAB implementation of the Cepstral Coefficients (CCs)-based Finite Element Model updating (FEMU) method using the customized LSHADE-SPACMA algorithm, a hybrid novel variant of LSHADE with CMA-ES integration. The framework is designed for efficient structural parameter identification in vibration-based structural health monitoring (SHM).**

# Key Features

1) Utilize CCs as the key cepstral characteristics for response comparison, effectively overcoming the 2-step modal updating inefficiency and decoupling stiffness–damping effects.

2) Open, reproducible, and extendable MATLAB code for SHM/FE model updating research.

# File Description

CCsFEMU_main_V1.m – Main driver program (run CC-based customized LSHADE-SPACMA FEMU algorithm).

response_simulation8d.m – Simulates structural acceleration responses.

extract_cc.m – Extracts cepstral coefficients from time-domain responses.

cc_update_cost.m – Cost function comparing experimental vs. simulated CCs.

exp_resp.mat – Example experimental dataset.

# Usage

Clone or download this repository.

Open MATLAB and run:

CCsFEMU_main_V1.m

The program will:

Load experimental response data.

Initialize structural parameters.

Run the developed algorithm for FEMU.

Print identified stiffness and damping ratios per DOF.

# Requirements

MATLAB R2021a or newer (tested).

No additional toolboxes required.

# Reference

If you use this code in your research, please cite:

1) Li L, Morgantini M, Betti R. Structural damage assessment through a new generalized autoencoder with features in the quefrency domain. Mechanical Systems and Signal Processing. 2023;184:109713.
2) Li L, Shen Z, Gan L, Xu L, Zhang H, Ye W. A novel cepstrum-based and evolutionary optimization hybrid framework for structural parameter identification. Advances in Structural Engineering. 2025;0(0).




