%% Calculate Norm of Each Local Hamiltonian
global L;
global H;
global H1;
global H2;
global H3;
global H_shift;
global H1_shift;
global H2_shift;
global H3_shift;

% Lanczos algorithm will return the eigenvalue with max norm
A = eigs(H1, 1);
B = eigs(H2, 1);
C = eigs(H3, 1);
t = 10.0;
N = [1000,10000,100000,200000,500000,1000000];

sum = A + B + C;
probabilities = [A/sum, B/sum, C/sum];

% Randomly Generate
qDRIFT_Error = cell(1,length(N));

%% qDRIFT
for k = 1:length(N)
    U_qDrift = speye(2^(2*L));
    U_1 = expm(-1i*H1*(t/N(k))*sum/A);
    U_2 = expm(-1i*H2*(t/N(k))*sum/B);
    U_3 = expm(-1i*H3*(t/N(k))*sum/C);
    for j = 1:N(k)
        choice = rand();
        if choice < probabilities(1)
            selected_option = 1;
            U_qDrift = U_qDrift * U_1;
        elseif choice < probabilities(1) + probabilities(2)
            selected_option = 2;
            U_qDrift = U_qDrift * U_2;
        else
            selected_option = 3;
            U_qDrift = U_qDrift * U_3;
        end
    end
    qDRIFT_Error{k}  = U_qDrift - expm(-1i*H*t);
end