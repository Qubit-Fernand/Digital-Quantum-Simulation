global N;
global H;
global H1;
global H2;
global H3;
global H_shift;
global H1_shift;
global H2_shift;
global H3_shift;

% Define the size of the system (a N = L * L square lattice with two L-length chains)
L = 3; % The length of the square lattice
N = L^2; % System size
alpha = 0; % Power of the distance

J = 1; % Exchange interaction energy

% Initialize the Hamiltonian matrices and coordinate matrix
H1 = sparse(2^N, 2^N); % Horizontal interaction
H2 = sparse(2^N, 2^N); % Vertical interaction
H3 = sparse(2^N, 2^N); % The remaining

% Pauli Matrix
X = sparse([0, 1; 1, 0]);
Y = sparse([0, -1i; 1i, 0]);
Z = sparse([1, 0; 0, -1]);
I = sparse([1, 0; 0, 1]);

%  Assign the coordinates of the sites
coords = cell(1, N); % The coordinates of the sites
for i = 1:N
    coords{i} = [mod(i-1, L)+1, floor((i-1)/L)+1];
end

% Construct the local Hamiltonians
% If L is not 3, the following needs to be modified
H1 = H1 - (J/1^alpha) * (kron(kron(kron(X, X), I), speye(2^(2*L))) + kron(kron(kron(Y, Y), I), speye(2^(2*L))) + kron(kron(kron(Z, Z), I), speye(2^(2*L))));
H1 = H1 - (J/1^alpha) * (kron(kron(I, kron(X, X)), speye(2^(2*L))) + kron(kron(I, kron(Y, Y)), speye(2^(2*L))) + kron(kron(I, kron(Z, Z)), speye(2^(2*L))));
H1 = H1 - (J/(2^alpha)) * (kron(kron(kron(X, I), X), speye(2^(2*L))) + kron(kron(kron(Y, I), Y), speye(2^(2*L))) + kron(kron(kron(Z, I), Z), speye(2^(2*L))));
H1 = H1 - (J/1^alpha) * (kron(speye(2^(2*L)), kron(kron(X, X), I)) + kron(speye(2^(2*L)), kron(kron(Y, Y), I)) + kron(speye(2^(2*L)), kron(kron(Z, Z), I)));
H1 = H1 - (J/1^alpha) * (kron(speye(2^(2*L)), kron(I, kron(X, X))) + kron(speye(2^(2*L)), kron(I, kron(Y, Y))) + kron(speye(2^(2*L)), kron(I, kron(Z, Z))));
H1 = H1 - (J/(2^alpha)) * (kron(speye(2^(2*L)), kron(kron(X, I), X)) + kron(speye(2^(2*L)), kron(kron(Y, I), Y)) + kron(speye(2^(2*L)), kron(kron(Z, I), Z)));
H1 = H1 - (J/1^alpha) * (kron(speye(2^L), kron(kron(kron(X, X), I), speye(2^L))) + kron(speye(2^L), kron(kron(kron(Y, Y), I), speye(2^L))) + kron(speye(2^L), kron(kron(kron(Z, Z), I), speye(2^L))));
H1 = H1 - (J/1^alpha) * (kron(speye(2^L), kron(kron(I, kron(X, X)), speye(2^L))) + kron(speye(2^L), kron(kron(I, kron(Y, Y)), speye(2^L))) + kron(speye(2^L), kron(kron(I, kron(Z, Z)), speye(2^L))));
H1 = H1 - (J/(2)^alpha) * (kron(speye(2^L), kron(kron(kron(X, I), X), speye(2^L))) + kron(speye(2^L), kron(kron(kron(Y, I), Y), speye(2^L))) + kron(speye(2^L), kron(kron(kron(Z, I), Z), speye(2^L))));

H2 = H2 - (J/1^alpha) * (kron(kron(kron(X, kron(I, I)), X), speye(2^(N-4))) + kron(kron(kron(Y, kron(I, I)), Y), speye(2^(N-4))) + kron(kron(kron(Z, kron(I, I)), Z), speye(2^(N-4))));
H2 = H2 - (J/1^alpha) * (kron(kron(kron(kron(I, X), speye(2^(L-1))), X), speye(2^(L+1))) + kron(kron(kron(kron(I, Y), speye(2^(L-1))), Y), speye(2^(L+1))) + kron(kron(kron(kron(I, Z), speye(2^(L-1))), Z), speye(2^(L+1))));
H2 = H2 - (J/1^alpha) * (kron(kron(kron(kron(speye(2^(L-1)), X), speye(2^(L-1))), X), speye(2^L)) + kron(kron(kron(kron(speye(2^(L-1)), Y), speye(2^(L-1))), Y), speye(2^L)) + kron(kron(kron(kron(speye(2^(L-1)), Z), speye(2^(L-1))), Z), speye(2^L)));
H2 = H2 - (J/1^alpha) * (kron(kron(kron(kron(speye(2^L), X), speye(2^(L-1))), X), speye(2^(L-1))) +  kron(kron(kron(kron(speye(2^L), Y), speye(2^(L-1))), Y), speye(2^(L-1))) +  kron(kron(kron(kron(speye(2^L), Z), speye(2^(L-1))), Z), speye(2^(L-1))));
H2 = H2 - (J/1^alpha) * (kron(kron(kron(kron(speye(2^(L+1)), X), speye(2^(L-1))), X), I) + kron(kron(kron(kron(speye(2^(L+1)), Y), speye(2^(L-1))), Y), I) + kron(kron(kron(kron(speye(2^(L+1)), Z), speye(2^(L-1))), Z), I));
H2 = H2 - (J/1^alpha) * (kron(kron(kron(speye(2^(L+2)), X), speye(2^(L-1))), X) + kron(kron(kron(speye(2^(L+2)), Y), speye(2^(L-1))), Y) + kron(kron(kron(speye(2^(L+2)), Z), speye(2^(L-1))), Z));
H2 = H2 - (J/(2^alpha)) * (kron(kron(kron(X, speye(2^(L+2))), X), speye(2^(L-1))) + kron(kron(kron(Y, speye(2^(L+2))), Y), speye(2^(L-1))) + kron(kron(kron(Z, speye(2^(L+2))), Z), speye(2^(L-1))));
H2 = H2 - (J/(2^alpha)) * (kron(kron(kron(kron(I, X), speye(2^(L+2))), X), I) + kron(kron(kron(kron(I, Y), speye(2^(L+2))), Y), I) + kron(kron(kron(kron(I, Z), speye(2^(L+2))), Z), I));
H2 = H2 - (J/(2^alpha)) * (kron(kron(kron(speye(2^(L-1)), X), speye(2^(L+2))), X) + kron(kron(kron(speye(2^(L-1)), Y), speye(2^(L+2))), Y) + kron(kron(kron(speye(2^(L-1)), Z), speye(2^(L+2))), Z));

H3 = H3 - (J/(sqrt(2)^alpha)) * (kron(kron(kron(X, speye(2^L)), X), speye(2^(L+1))) + kron(kron(kron(Y, speye(2^L)), Y), speye(2^(L+1))) + kron(kron(kron(Z, speye(2^L)), Z), speye(2^(L+1))));
H3 = H3 - (J/(sqrt(2)^alpha)) * (kron(kron(kron(kron(I, X), I), X), speye(2^(L+2))) + kron(kron(kron(kron(I, Y), I), Y), speye(2^(L+2))) + kron(kron(kron(kron(I, Z), I), Z), speye(2^(L+2))));
H3 = H3 - (J/(sqrt(2)^alpha)) * (kron(kron(kron(kron(I, X), speye(2^L)), X), speye(2^L)) + kron(kron(kron(kron(I, Y), speye(2^L)), Y), speye(2^L)) + kron(kron(kron(kron(I, Z), speye(2^L)), Z), speye(2^L)));
H3 = H3 - (J/(sqrt(2)^alpha)) * (kron(kron(kron(kron(speye(2^(L-1)), X), I), X), speye(2^(L+1))) + kron(kron(kron(kron(speye(2^(L-1)), Y), I), Y), speye(2^(L+1))) + kron(kron(kron(kron(speye(2^(L-1)), Z), I), Z), speye(2^(L+1))));
H3 = H3 - (J/(sqrt(2)^alpha)) * (kron(kron(kron(kron(speye(2^L), X), speye(2^L)), X), I) + kron(kron(kron(kron(speye(2^L), Y), speye(2^L)), Y), I) + kron(kron(kron(kron(speye(2^L), Z), speye(2^L)), Z), I));
H3 = H3 - (J/(sqrt(2)^alpha)) * (kron(kron(kron(kron(speye(2^(L+1)), X), I), X), speye(2^(L-1))) + kron(kron(kron(kron(speye(2^(L+1)), Y), I), Y), speye(2^(L-1))) + kron(kron(kron(kron(speye(2^(L+1)), Z), I), Z), speye(2^(L-1))));
H3 = H3 - (J/(sqrt(2)^alpha)) * (kron(kron(kron(speye(2^(L+1)), X), speye(2^L)), X) + kron(kron(kron(speye(2^(L+1)), Y), speye(2^L)), Y) + kron(kron(kron(speye(2^(L+1)), Z), speye(2^L)), Z));
H3 = H3 - (J/(sqrt(2)^alpha)) * (kron(kron(kron(kron(speye(2^(L+2)), X), I), X), I) + kron(kron(kron(kron(speye(2^(L+2)), Y), I), Y), I) + kron(kron(kron(kron(speye(2^(L+2)), Z), I), Z), I));

H3 = H3 - (J/(sqrt(5)^alpha)) * (kron(kron(kron(X, speye(2^(L+1))), X), speye(2^L)) + kron(kron(kron(Y, speye(2^(L+1))), Y), speye(2^L)) + kron(kron(kron(Z, speye(2^(L+1))), Z), speye(2^L)));
H3 = H3 - (J/(sqrt(5)^alpha)) * (kron(kron(kron(I, I), kron(X, X)), speye(2^(L+2))) + kron(kron(kron(I, I), kron(Y, Y)), speye(2^(L+2))) + kron(kron(kron(I, I), kron(Z, Z)), speye(2^(L+2))));
H3 = H3 - (J/(sqrt(5)^alpha)) * (kron(kron(kron(speye(2^L), kron(I, I)),kron(X, X)), kron(I, I)) + kron(kron(kron(speye(2^L), kron(I, I)),kron(Y, Y)), kron(I, I)) + kron(kron(kron(speye(2^L), kron(I, I)),kron(Z, Z)), kron(I, I)));
H3 = H3 - (J/(sqrt(5)^alpha)) * (kron(kron(kron(speye(2^L), X), speye(2^(L+1))), X) + kron(kron(kron(speye(2^L), Y), speye(2^(L+1))), Y) + kron(kron(kron(speye(2^L), Z), speye(2^(L+1))), Z));
H3 = H3 - (J/(sqrt(5)^alpha)) * (kron(kron(kron(X, speye(2^(2*L))), X), I) + kron(kron(kron(Y, speye(2^(2*L))), Y), I) + kron(kron(kron(Z, speye(2^(2*L))), Z), I));
H3 = H3 - (J/(sqrt(5)^alpha)) * (kron(kron(kron(kron(I, X), speye(2^(L+1))), X), kron(I, I)) + kron(kron(kron(kron(I, Y), speye(2^(L+1))), Y), kron(I, I)) + kron(kron(kron(kron(I, Z), speye(2^(L+1))), Z), kron(I, I)));
H3 = H3 - (J/(sqrt(5)^alpha)) * (kron(kron(kron(I, X), speye(2^(2*L))), X) + kron(kron(kron(I, Y), speye(2^(2*L))), Y) + kron(kron(kron(I, Z), speye(2^(2*L))), Z));
H3 = H3 - (J/(sqrt(5)^alpha)) * (kron(kron(kron(kron(kron(I, I), X), speye(2^(L+1))), X), I) + kron(kron(kron(kron(kron(I, I), Y), speye(2^(L+1))), Y), I) + kron(kron(kron(kron(kron(I, I), Z), speye(2^(L+1))), Z), I));

H3 = H3 - (J/(sqrt(8)^alpha)) * (kron(kron(X, speye(2^(2*L + 1))), X) + kron(kron(Y, speye(2^(2*L + 1))), Y) + kron(kron(Z, speye(2^(2*L + 1))), Z));
H3 = H3 - (J/(sqrt(8)^alpha)) * (kron(kron(kron(kron(kron(I, I), X), speye(2^L)), X), kron(I, I)) + kron(kron(kron(kron(kron(I, I), Y), speye(2^L)), Y), kron(I, I)) + kron(kron(kron(kron(kron(I, I), Z), speye(2^L)), Z), kron(I, I)));

% Shift the Hamiltonian by the ground state energy
energy1 = sort(eigs(H1, 2^N));
energy2 = sort(eigs(H2, 2^N));
energy3 = sort(eigs(H3, 2^N));
ground_energy1 = min(energy1);
ground_energy2 = min(energy2);
ground_energy3 = min(energy3);
H1_shift = H1 - floor(ground_energy1) * speye(2^N);
H2_shift = H2 - floor(ground_energy2) * speye(2^N);
H3_shift = H3 - floor(ground_energy3) * speye(2^N);

% Derive the global Hamiltonian
H = H1 + H2 + H3;
H_shift = H1_shift+H2_shift+H3_shift;
[states, energy] = eigs(H_shift, 2^N);
energy = diag(energy);