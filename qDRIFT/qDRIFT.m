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
num_list = [1000,10000,100000,200000,500000,1000000];

sum = A + B + C;
probabilities = [A/sum, B/sum, C/sum];

%% qDRIFT: t = 1
t = 1.0;
qDRIFT_Error_1 = cell(1,length(num_list));

for k = 1:length(num_list)
    U_qDrift = speye(2^N);
    U_1 = expm(-1i*H1_shift*(t/num_list(k))*sum/A);
    U_2 = expm(-1i*H2_shift*(t/num_list(k))*sum/B);
    %% 
    U_3 = expm(-1i*H3_shift*(t/num_list(k))*sum/C);
    for j = 1:num_list(k)
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
    qDRIFT_Error_1{k}  = U_qDrift - expm(-1i*H_shift*t);
end

%% qDRIFT: t = 3
t = 3.0;
qDRIFT_Error_3 = cell(1,length(num_list));

for k = 1:length(num_list)
    U_qDrift = speye(2^N);
    U_1 = expm(-1i*H1_shift*(t/num_list(k))*sum/A);
    U_2 = expm(-1i*H2_shift*(t/num_list(k))*sum/B);
    U_3 = expm(-1i*H3_shift*(t/num_list(k))*sum/C);
    for j = 1:num_list(k)
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
    qDRIFT_Error_3{k}  = U_qDrift - expm(-1i*H_shift*t);
end

%% qDRIFT: t = 5
t = 5.0;
qDRIFT_Error_5 = cell(1,length(num_list));

for k = 1:length(num_list)
    U_qDrift = speye(2^N);
    U_1 = expm(-1i*H1_shift*(t/num_list(k))*sum/A);
    U_2 = expm(-1i*H2_shift*(t/num_list(k))*sum/B);
    U_3 = expm(-1i*H3_shift*(t/num_list(k))*sum/C);
    for j = 1:num_list(k)
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
    qDRIFT_Error_5{k}  = U_qDrift - expm(-1i*H_shift*t);
end

%% qDRIFT: t = 10
t = 10.0;
qDRIFT_Error_10 = cell(1,length(num_list));

for k = 1:length(num_list)
    U_qDrift = speye(2^N);
    U_1 = expm(-1i*H1_shift*(t/num_list(k))*sum/A);
    U_2 = expm(-1i*H2_shift*(t/num_list(k))*sum/B);
    U_3 = expm(-1i*H3_shift*(t/num_list(k))*sum/C);
    for j = 1:num_list(k)
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
    qDRIFT_Error_10{k}  = U_qDrift - expm(-1i*H_shift*t);
end