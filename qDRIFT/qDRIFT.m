%% Calculate Norm of Each Local Hamiltonian
global N;
global H;
global H1;
global H2;
global H3;
global H_shift;
global H1_shift;
global H2_shift;
global H3_shift;

% Lanczos algorithm will return the eigenvalue with max norm
A = eigs(H1_shift, 1);
B = eigs(H2_shift, 1);
C = eigs(H3_shift, 1);

t_list = [1.0, 3.0, 5.0, 10.0];
r_list = [10000,50000,100000,200000,500000,1000000];
c_list = cell(1,length(r_list));

sum = A + B + C;
probabilities = [A/sum, B/sum, C/sum];

for i = 1:length(r_list)
    c_list{i} = zeros(1, r_list(i));
    for j = 1:r_list(i)
        choice = rand();
        if choice < probabilities(1)
            c_list{i}(j) = 1;
        elseif choice < probabilities(1) + probabilities(2)
            c_list{i}(j) = 2;
        else
            c_list{i}(j) = 3;
        end
    end
end

% Please load('qDRIFT.mat') to fix the choice sequence for every system
load('qDRIFT.mat');
%% qDRIFT: t = 1
t = 1.0;
qDRIFT_Error_1 = cell(1,length(r_list));

for i = 1:length(r_list)
    U_qDrift = speye(2^N);
    U_1 = expm(-1i*H1_shift*(t/r_list(i))*sum/A);
    U_2 = expm(-1i*H2_shift*(t/r_list(i))*sum/B);
    U_3 = expm(-1i*H3_shift*(t/r_list(i))*sum/C);
    U_list = {U_1, U_2, U_3};
    for j = 1:r_list(i)
        U_qDrift = U_qDrift * U_list{c_list{i}(j)};
    end
    qDRIFT_Error_1{i}  = U_qDrift - expm(-1i*H_shift*t);
end

%% qDRIFT: t = 3
t = 3.0;
qDRIFT_Error_3 = cell(1,length(r_list));

for i = 1:length(r_list)
    U_qDrift = speye(2^N);
    U_1 = expm(-1i*H1_shift*(t/r_list(i))*sum/A);
    U_2 = expm(-1i*H2_shift*(t/r_list(i))*sum/B);
    U_3 = expm(-1i*H3_shift*(t/r_list(i))*sum/C);
    U_list = {U_1, U_2, U_3};
    for j = 1:r_list(i)
        U_qDrift = U_qDrift * U_list{c_list{i}(j)};
    end
    qDRIFT_Error_3{i}  = U_qDrift - expm(-1i*H_shift*t);
end

%% qDRIFT: t = 5
t = 5.0;
qDRIFT_Error_5 = cell(1,length(r_list));

for i = 1:length(r_list)
    U_qDrift = speye(2^N);
    U_1 = expm(-1i*H1_shift*(t/r_list(i))*sum/A);
    U_2 = expm(-1i*H2_shift*(t/r_list(i))*sum/B);
    U_3 = expm(-1i*H3_shift*(t/r_list(i))*sum/C);
    U_list = {U_1, U_2, U_3};
    for j = 1:r_list(i)
        U_qDrift = U_qDrift * U_list{c_list{i}(j)};
    end
    qDRIFT_Error_5{i}  = U_qDrift - expm(-1i*H_shift*t);
end

%% qDRIFT: t = 10
t = 10.0;
qDRIFT_Error_10 = cell(1,length(r_list));

for i = 1:length(r_list)
    U_qDrift = speye(2^N);
    U_1 = expm(-1i*H1_shift*(t/r_list(i))*sum/A);
    U_2 = expm(-1i*H2_shift*(t/r_list(i))*sum/B);
    U_3 = expm(-1i*H3_shift*(t/r_list(i))*sum/C);
    U_list = {U_1, U_2, U_3};
    for j = 1:r_list(i)
        U_qDrift = U_qDrift * U_list{c_list{i}(j)};
    end
    qDRIFT_Error_10{i}  = U_qDrift - expm(-1i*H_shift*t);
end