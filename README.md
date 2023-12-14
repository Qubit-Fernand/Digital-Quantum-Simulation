# Quantum-Simulation
Numerical Experiments for Quantum Simulation Algorithms: Trotter-Suzuki, qDRIFT, Randomized Permutaion and Symmetry Transformation.

The "Hamiltonian" folder contains different types of models: single chain homogeneous & 2x6 Heisenberg and power-law models in which each local term is shifted to semidefinite for subspace projection analysis. 

In each folder the "test.m" will process matrix data to operator norm according to given low-energy threshold. For instance, "qDRIFT/raw_1x8.mat" denotes matrix error for 1x8 homogeneous Heisenberg model induced by the qDRIFT algorithm. Correspondingly, "qDRIFT/norm_1x8_14.mat" denotes the operator norm in full space as well as low-energy subspace with \Delta = 14.

**Data Availability:** Original ".mat" data used for plotting figures in our paper is too large (>1GB) to be pushed to github. If you are interested, please contact AntiEntropy@pku.edu.cn upon reasonable request.
