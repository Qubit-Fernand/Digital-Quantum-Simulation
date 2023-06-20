% Calculate Norm of Each Local Hamiltonian
global L;
global H_shift;
global H1_shift;
global H2_shift;
global H3_shift;

A = eigs(H1_shift, 1);
B = eigs(H2_shift, 1);
C = eigs(H3_shift, 1);
N = 100;
dt = [0.005 0.008 0.010 0.012 0.015 0.017 0.020];

sum = A + B + C;
probabilities = [A/sum, B/sum, C/sum];

% Randomly Generate
qDRIFT_shift_Error_N_50 = cell(1,length(dt));
qDRIFT_shift_Error_N_100 = cell(1,length(dt));

for k = 1:length(dt)
    U_qDrift = speye(2^(2*L));
    for j = 1:N
        choice = rand();
        if choice < probabilities(1)
            selected_option = 1;
            U_qDrift = U_qDrift * expm(-1i*H1_shift*dt(k));
        elseif choice < probabilities(1) + probabilities(2)
            selected_option = 2;
            U_qDrift = U_qDrift * expm(-1i*H2_shift*dt(k));
        else
            selected_option = 3;
            U_qDrift = U_qDrift * expm(-1i*H3_shift*dt(k));
        end
    end
    qDRIFT_shift_Error_N_100{k}  = U_qDrift - expm(-1i*H_shift*N*dt(k));
end
