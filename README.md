CEPSTRUM-BASED CUSTOMIZED LSHADE-SPACMA HYBRID FEMU FRAMEWORK (MATLAB Implementation)

This repository provides the MATLAB implementation of the Cepstral Coefficients (CCs)-based Finite Element Model updating (FEMU) method using the customized LSHADE-SPACMA algorithm, a hybrid novel variant of LSHADE with CMA-ES integration. The framework is designed for efficient structural parameter identification in vibration-based structural health monitoring (SHM).

Key Features

Utilize CCs as the key cepstral characteristics for response comparison, effectively overcoming the 2-step modal updating inefficiency and decoupling stiffness–damping effects.

Open, reproducible, and extendable MATLAB code for SHM/FE model updating research.

File Description

Main.m – Main driver program (run CC-based customized LSHADE-SPACMA FEMU algorithm).

response_simulation8d.m – Simulates structural acceleration responses.

extract_cc.m – Extracts cepstral coefficients from time-domain responses.

cc_update_cost.m – Cost function comparing experimental vs. simulated CCs.

exp_resp_baseline.mat – Example experimental dataset.

Usage

Clone or download this repository.

Open MATLAB and run:

Main


The program will:

Load experimental response data.

Initialize structural parameters.

Run the developed algorithm for FEMU.

Print identified stiffness and damping ratios per DOF.


Requirements

MATLAB R2021a or newer (tested).

No additional toolboxes required.

Reference

If you use this code in your research, please cite:

[Your Paper Title Here]
[Your Name(s)], [Journal/Conference], [Year].