global L
global H
global H1
global H2
global H3
global H_shift
global H1_shift
global H2_shift
global H3_shift
 
% Define the parameters of the Heisenberg model
% Here we suppose the external magnetic field is zero
J = 1; % Exchange interaction energy

% Define the size of the system (a 2 × L lattice with two L-length chains)
L = 4;

% Construct the Hamiltonian matrix, H1: intra-chain interaction & H2: inter-chain interaction
i = 1;
% Define the exchange interaction term between adjacent spins (1D chain)
SxSx1 = kron(kron(speye(2^(i-1)),kron(sparse([0, 1; 1, 0]), sparse([0, 1; 1, 0]))), speye(2^(L-i-1)));
SySy1 = kron(kron(speye(2^(i-1)),kron(sparse([0, -1i; 1i, 0]), sparse([0, -1i; 1i, 0]))), speye(2^(L-i-1)));
SzSz1 = kron(kron(speye(2^(i-1)),kron(sparse([1, 0; 0, -1]), sparse([1, 0; 0, -1]))), speye(2^(L-i-1)));
% Add the exchange interaction term to the Hamiltonian matrix
H1 = - J * (SxSx1 + SySy1 + SzSz1);

i = 2;
% Define the exchange interaction term between adjacent spins (1D chain)
SxSx2 = kron(kron(speye(2^(i-1)),kron(sparse([0, 1; 1, 0]), sparse([0, 1; 1, 0]))), speye(2^(L-i-1)));
SySy2 = kron(kron(speye(2^(i-1)),kron(sparse([0, -1i; 1i, 0]), sparse([0, -1i; 1i, 0]))), speye(2^(L-i-1)));
SzSz2 = kron(kron(speye(2^(i-1)),kron(sparse([1, 0; 0, -1]), sparse([1, 0; 0, -1]))), speye(2^(L-i-1)));
% Add the exchange interaction term to the Hamiltonian matrix
H2 = - J * (SxSx2 + SySy2 + SzSz2);

i = 3;
% Define the exchange interaction term between adjacent spins (1D chain)
SxSx3 = kron(kron(speye(2^(i-1)),kron(sparse([0, 1; 1, 0]), sparse([0, 1; 1, 0]))), speye(2^(L-i-1)));
SySy3 = kron(kron(speye(2^(i-1)),kron(sparse([0, -1i; 1i, 0]), sparse([0, -1i; 1i, 0]))), speye(2^(L-i-1)));
SzSz3 = kron(kron(speye(2^(i-1)),kron(sparse([1, 0; 0, -1]), sparse([1, 0; 0, -1]))), speye(2^(L-i-1)));
% Add the exchange interaction term to the Hamiltonian matrix
H3 = - J * (SxSx3 + SySy3 + SzSz3);

% Shift the Hamiltonian by the ground state energy
% Lanczos Method for sparse matrix. By default, the output is sorted in descending order of absolute value
energy1 = sort(eigs(H1, 16)); 
energy2 = sort(eigs(H2, 16)); 
energy3 = sort(eigs(H3, 16));
ground_energy1 = min(energy1);
ground_energy2 = min(energy2);
ground_energy3 = min(energy3);

% Floor function removed
H1_shift = H1 - ground_energy1 * speye(2^L);
H2_shift = H2 - ground_energy2 * speye(2^L);
H3_shift = H3 - ground_energy3 * speye(2^L);

% Derive the global Hamiltonian
H = H1 + H2 + H3;
H_shift = H1_shift+H2_shift+H3_shift;
[states, energy] = eigs(H_shift, 16);
energy = diag(energy);