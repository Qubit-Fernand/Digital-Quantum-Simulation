% Calculate Norm of Each Local Hamiltonian
global L;
global H;
global H1;
global H2;
global H3;

A = eigs(H1, 1);
B = eigs(H2, 1);
C = eigs(H3, 1);
N = 10;
dt = 0.010;

sum = A + B + C;
probabilities = [A/sum, B/sum, C/sum];

% Randomly Generate
U_qDrift = speye(2^(2*L));

for i = 1:N
    choice = rand();
    if choice < probabilities(1)
        selected_option = 1;
        U_qDrift = U_qDrift * exp(-1i*H1*dt);
    elseif choice < probabilities(1) + probabilities(2)
        selected_option = 2;
        U_qDrift = U_qDrift * exp(-1i*H2*dt);
    else
        selected_option = 3;
        U_qDrift = U_qDrift * exp(-1i*H3*dt);
    end
end

qDRIFT_error = U_qDrift - exp(-1i*H*dt*N);