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

addpath('../Hamiltonian/')

t_list = [1.0, 3.0, 5.0, 10.0];
num_list = [1000,10000,100000,200000,500000,1000000];
step_number = num_list(5);

%% N = 2x5
L = 5;
run('Heisenberg.m')
% clearvars -except N;
% Lanczos algorithm will return the eigenvalue with max norm
A = eigs(H1_shift, 1);
B = eigs(H2_shift, 1);
C = eigs(H3_shift, 1);
sum = A + B + C;
probabilities = [A/sum, B/sum, C/sum];

qDRIFT_Error_2x5 = cell(1,length(t_list));

for i = 1:length(t_list)
    t = t_list(i);
    U_qDrift = speye(2^N);
    for j = 1:step_number
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
    qDRIFT_Error_2x5{i} = norm(U_qDrift - expm(-1i*H_shift*t));
end

%% N = 2x4
L = 4;
run('Heisenberg.m')
% clearvars -except N;
% Lanczos algorithm will return the eigenvalue with max norm
A = eigs(H1_shift, 1);
B = eigs(H2_shift, 1);
C = eigs(H3_shift, 1);
sum = A + B + C;
probabilities = [A/sum, B/sum, C/sum];

qDRIFT_Error_2x4 = cell(1,length(t_list));

for i = 1:length(t_list)
    t = t_list(i);
    U_qDrift = speye(2^N);
    for j = 1:step_number
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
    qDRIFT_Error_2x4{i} = norm(U_qDrift - expm(-1i*H_shift*t));
end

%% N = 2x3
L = 3;
run('Heisenberg.m')
% clearvars -except N;
% Lanczos algorithm will return the eigenvalue with max norm
A = eigs(H1_shift, 1);
B = eigs(H2_shift, 1);
C = eigs(H3_shift, 1);
sum = A + B + C;
probabilities = [A/sum, B/sum, C/sum];

qDRIFT_Error_2x3 = cell(1,length(t_list));

for i = 1:length(t_list)
    t = t_list(i);
    U_qDrift = speye(2^N);
    for j = 1:step_number
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
    qDRIFT_Error_2x3{i} = norm(U_qDrift - expm(-1i*H_shift*t));
end

%% N = 2x2
L = 2;
run('Heisenberg.m')
% clearvars -except N;
% Lanczos algorithm will return the eigenvalue with max norm
A = eigs(H1_shift, 1);
B = eigs(H2_shift, 1);
C = eigs(H3_shift, 1);
sum = A + B + C;
probabilities = [A/sum, B/sum, C/sum];

qDRIFT_Error_2x2 = cell(1,length(t_list));

for i = 1:length(t_list)
    t = t_list(i);
    U_qDrift = speye(2^N);
    for j = 1:step_number
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
    qDRIFT_Error_2x2{i} = norm(U_qDrift - expm(-1i*H_shift*t));
end