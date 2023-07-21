global N;
global H;
global H1;
global H2;
global H3;
global H_shift;
global H1_shift;
global H2_shift;
global H3_shift;

% Define the parameters of the Heisenberg model
% Here we suppose the external magnetic field is zero

% Define the size of the system (a N = 2 * L lattice with two L-length chains)
L = 5; % Chain length: L
N = 2 * L; % System size

J = 1; % Exchange interaction energy

% Construct the Hamiltonian matrix, H1: intra-chain interaction & H2: inter-chain interaction
H1 = sparse(2^N,2^N);
H2 = sparse(2^N,2^N);
H3 = sparse(2^N,2^N);

% Pauli Matrix
X = sparse([0, 1; 1, 0]);
Y = sparse([0, -1i; 1i, 0]);
Z = sparse([1, 0; 0, -1]);
I = sparse([1, 0; 0, 1]);

for i = 1:2:L-1
    % Define the exchange interaction term between adjacent spins (1D chain)
    H1 = H1 - J * kron(kron(kron(speye(2^(i-1)),kron(X, X)), speye(2^(L-i-1))), speye(2^L));
    H1 = H1 - J * kron(kron(kron(speye(2^(i-1)),kron(Y, Y)), speye(2^(L-i-1))), speye(2^L));
    H1 = H1 - J * kron(kron(kron(speye(2^(i-1)),kron(Z, Z)), speye(2^(L-i-1))), speye(2^L));
    % Define the exchange interaction term between adjacent spins (1D chain)
    H1 = H1 - J * kron(speye(2^L), kron(kron(speye(2^(i-1)),kron(X, X)), speye(2^(L-i-1))));
    H1 = H1 - J * kron(speye(2^L), kron(kron(speye(2^(i-1)),kron(Y, Y)), speye(2^(L-i-1))));
    H1 = H1 - J * kron(speye(2^L), kron(kron(speye(2^(i-1)),kron(Z, Z)), speye(2^(L-i-1))));
end

for i = 2:2:L-1
    % Define the exchange interaction term between adjacent spins (1D chain)
    H2 = H2 - J * kron(kron(kron(speye(2^(i-1)),kron(X, X)), speye(2^(L-i-1))), speye(2^L));
    H2 = H2 - J * kron(kron(kron(speye(2^(i-1)),kron(Y, Y)), speye(2^(L-i-1))), speye(2^L));
    H2 = H2 - J * kron(kron(kron(speye(2^(i-1)),kron(Z, Z)), speye(2^(L-i-1))), speye(2^L));
    % Define the exchange interaction term between adjacent spins (1D chain)
    H2 = H2 - J * kron(speye(2^L), kron(kron(speye(2^(i-1)),kron(X, X)), speye(2^(L-i-1))));
    H2 = H2 - J * kron(speye(2^L), kron(kron(speye(2^(i-1)),kron(Y, Y)), speye(2^(L-i-1))));
    H2 = H2 - J * kron(speye(2^L), kron(kron(speye(2^(i-1)),kron(Z, Z)), speye(2^(L-i-1))));
end

for i = 1:L
    % Define the exchange interaction term between adjacent spins (inter-chain)
    H3 = H3 - J * kron(kron(kron(speye(2^(i-1)), X), speye(2^(L-i))), kron(kron(speye(2^(i-1)), X), speye(2^(L-i))));
    H3 = H3 - J * kron(kron(kron(speye(2^(i-1)), Y), speye(2^(L-i))), kron(kron(speye(2^(i-1)), Y), speye(2^(L-i))));
    H3 = H3 - J * kron(kron(kron(speye(2^(i-1)), Z), speye(2^(L-i))), kron(kron(speye(2^(i-1)), Z), speye(2^(L-i))));
end

% Shift the Hamiltonian by the ground state energy
% Lanczos Method for sparse matrix. By default, the output is sorted in
% descending order of norm, and the return value is a complex number.
energy1 = sort(eigs(H1, 2^N));
energy2 = sort(eigs(H2, 2^N));
energy3 = sort(eigs(H3, 2^N));
ground_energy1 = min(energy1);
ground_energy2 = min(energy2);
ground_energy3 = min(energy3);

% Floor function removed
H1_shift = H1 - ground_energy1 * speye(2^N);
H2_shift = H2 - ground_energy2 * speye(2^N);
H3_shift = H3 - ground_energy3 * speye(2^N);

% Derive the global Hamiltonian
H = H1 + H2 + H3;
H_shift = H1_shift+H2_shift+H3_shift;
[states, energy] = eigs(H_shift, 2^N);
energy = diag(energy);