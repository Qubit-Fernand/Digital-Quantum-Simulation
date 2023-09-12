%% Random Permutation 2nd-order Trotter-Suzuki
global N;
global H;
global H1;
global H2;
global H3;
global H_shift;
global H1_shift;
global H2_shift;
global H3_shift;

t_list = [1.0 2.0];
r_list = [10,50,100,200,500,1000];
c_list = cell(1,length(r_list));

for k = 1:length(r_list)
    c_list{k} = zeros(1, r_list(k));
    for j = 1:r_list(k)
        choice = randi([1, 6]);
        c_list{k}(j) = choice;
    end
end

Random_Trotter_Error_p_2_1 = cell(1,length(r_list));
Random_Trotter_Error_p_2_2 = cell(1,length(r_list));
Random_Trotter_Error_p_4_1 = cell(1,length(r_list));
Random_Trotter_Error_p_4_2 = cell(1,length(r_list));

% Please load('Permutation.mat') to fix the choice sequence for every system
load('Permutation.mat');
%% t = 1.0 p = 2;
t = 1.0;
for i = 1:length(r_list)
    U_permutation = speye(2^N);
    U_1 = expm(-1i*H1_shift*(t/r_list(i))/2) * expm(-1i*H2_shift*(t/r_list(i))/2) * expm(-1i*H3_shift*(t/r_list(i))) * expm(-1i*H2_shift*(t/r_list(i))/2) * expm(-1i*H1_shift*(t/r_list(i))/2);
    U_2 = expm(-1i*H1_shift*(t/r_list(i))/2) * expm(-1i*H3_shift*(t/r_list(i))/2) * expm(-1i*H2_shift*(t/r_list(i))) * expm(-1i*H3_shift*(t/r_list(i))/2) * expm(-1i*H1_shift*(t/r_list(i))/2);
    U_3 = expm(-1i*H2_shift*(t/r_list(i))/2) * expm(-1i*H1_shift*(t/r_list(i))/2) * expm(-1i*H3_shift*(t/r_list(i))) * expm(-1i*H1_shift*(t/r_list(i))/2) * expm(-1i*H2_shift*(t/r_list(i))/2);
    U_4 = expm(-1i*H2_shift*(t/r_list(i))/2) * expm(-1i*H3_shift*(t/r_list(i))/2) * expm(-1i*H1_shift*(t/r_list(i))) * expm(-1i*H3_shift*(t/r_list(i))/2) * expm(-1i*H2_shift*(t/r_list(i))/2);
    U_5 = expm(-1i*H3_shift*(t/r_list(i))/2) * expm(-1i*H2_shift*(t/r_list(i))/2) * expm(-1i*H1_shift*(t/r_list(i))) * expm(-1i*H2_shift*(t/r_list(i))/2) * expm(-1i*H3_shift*(t/r_list(i))/2);
    U_6 = expm(-1i*H3_shift*(t/r_list(i))/2) * expm(-1i*H1_shift*(t/r_list(i))/2) * expm(-1i*H2_shift*(t/r_list(i))) * expm(-1i*H1_shift*(t/r_list(i))/2) * expm(-1i*H3_shift*(t/r_list(i))/2);
    U_list = {U_1, U_2, U_3, U_4, U_5, U_6};
    for j = 1:r_list(i)
        U_permutation = U_permutation * U_list{c_list{i}(j)};
    end
    Random_Trotter_Error_p_2_1{i} = U_permutation - expm(-1i * H_shift * t);
end

%% t = 1.0 p = 4;
for i = 1:length(r_list)
    u = 1/(4-power(4,1/3));
    U_permutation = speye(2^N);
    U_1a = expm(-1i*H1_shift*u*(t/r_list(i))/2) * expm(-1i*H2_shift*u*(t/r_list(i))/2) * expm(-1i*H3_shift*u*(t/r_list(i))) * expm(-1i*H2_shift*u*(t/r_list(i))/2) * expm(-1i*H1_shift*u*(t/r_list(i))/2);
    U_1b = expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))) * expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))/2);
    U_1 = U_1a * U_1a * U_1b * U_1a * U_1a;
    U_2a = expm(-1i*H1_shift*u*(t/r_list(i))/2) * expm(-1i*H3_shift*u*(t/r_list(i))/2) * expm(-1i*H2_shift*u*(t/r_list(i))) * expm(-1i*H3_shift*u*(t/r_list(i))/2) * expm(-1i*H1_shift*u*(t/r_list(i))/2);
    U_2b = expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))) * expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))/2);
    U_2 = U_2a * U_2a * U_2b * U_2a * U_2a;
    U_3a = expm(-1i*H2_shift*u*(t/r_list(i))/2) * expm(-1i*H1_shift*u*(t/r_list(i))/2) * expm(-1i*H3_shift*u*(t/r_list(i))) * expm(-1i*H1_shift*u*(t/r_list(i))/2) * expm(-1i*H2_shift*u*(t/r_list(i))/2);
    U_3b = expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))) * expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))/2);
    U_3 = U_3a * U_3a * U_3b * U_3a * U_3a;
    U_4a = expm(-1i*H2_shift*u*(t/r_list(i))/2) * expm(-1i*H3_shift*u*(t/r_list(i))/2) * expm(-1i*H1_shift*u*(t/r_list(i))) * expm(-1i*H3_shift*u*(t/r_list(i))/2) * expm(-1i*H2_shift*u*(t/r_list(i))/2);
    U_4b = expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))) * expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))/2);
    U_4 = U_4a * U_4a * U_4b * U_4a * U_4a;
    U_5a = expm(-1i*H3_shift*u*(t/r_list(i))/2) * expm(-1i*H2_shift*u*(t/r_list(i))/2) * expm(-1i*H1_shift*u*(t/r_list(i))) * expm(-1i*H2_shift*u*(t/r_list(i))/2) * expm(-1i*H3_shift*u*(t/r_list(i))/2);
    U_5b = expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))) * expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))/2);
    U_5 = U_5a * U_5a * U_5b * U_5a * U_5a;
    U_6a = expm(-1i*H3_shift*u*(t/r_list(i))/2) * expm(-1i*H1_shift*u*(t/r_list(i))/2) * expm(-1i*H2_shift*u*(t/r_list(i))) * expm(-1i*H1_shift*u*(t/r_list(i))/2) * expm(-1i*H3_shift*u*(t/r_list(i))/2);
    U_6b = expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))) * expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))/2);
    U_6 = U_6a * U_6a * U_6b * U_6a * U_6a;
    U_list = {U_1, U_2, U_3, U_4, U_5, U_6};
    for j = 1:r_list(i)
        U_permutation = U_permutation * U_list{c_list{i}(j)};
    end
    Random_Trotter_Error_p_4_1{i} = U_permutation - expm(-1i * H_shift * t);
end

%% t = 2.0 p = 2;
t = 2.0;
for i = 1:length(r_list)
    U_permutation = speye(2^N);
    U_1 = expm(-1i*H1_shift*(t/r_list(i))/2) * expm(-1i*H2_shift*(t/r_list(i))/2) * expm(-1i*H3_shift*(t/r_list(i))) * expm(-1i*H2_shift*(t/r_list(i))/2) * expm(-1i*H1_shift*(t/r_list(i))/2);
    U_2 = expm(-1i*H1_shift*(t/r_list(i))/2) * expm(-1i*H3_shift*(t/r_list(i))/2) * expm(-1i*H2_shift*(t/r_list(i))) * expm(-1i*H3_shift*(t/r_list(i))/2) * expm(-1i*H1_shift*(t/r_list(i))/2);
    U_3 = expm(-1i*H2_shift*(t/r_list(i))/2) * expm(-1i*H1_shift*(t/r_list(i))/2) * expm(-1i*H3_shift*(t/r_list(i))) * expm(-1i*H1_shift*(t/r_list(i))/2) * expm(-1i*H2_shift*(t/r_list(i))/2);
    U_4 = expm(-1i*H2_shift*(t/r_list(i))/2) * expm(-1i*H3_shift*(t/r_list(i))/2) * expm(-1i*H1_shift*(t/r_list(i))) * expm(-1i*H3_shift*(t/r_list(i))/2) * expm(-1i*H2_shift*(t/r_list(i))/2);
    U_5 = expm(-1i*H3_shift*(t/r_list(i))/2) * expm(-1i*H2_shift*(t/r_list(i))/2) * expm(-1i*H1_shift*(t/r_list(i))) * expm(-1i*H2_shift*(t/r_list(i))/2) * expm(-1i*H3_shift*(t/r_list(i))/2);
    U_6 = expm(-1i*H3_shift*(t/r_list(i))/2) * expm(-1i*H1_shift*(t/r_list(i))/2) * expm(-1i*H2_shift*(t/r_list(i))) * expm(-1i*H1_shift*(t/r_list(i))/2) * expm(-1i*H3_shift*(t/r_list(i))/2);
    U_list = {U_1, U_2, U_3, U_4, U_5, U_6};
    for j = 1:r_list(i)
        U_permutation = U_permutation * U_list{c_list{i}(j)};
    end
    Random_Trotter_Error_p_2_2{i} = U_permutation - expm(-1i * H_shift * t);
end

%% t = 2.0 p = 4
for i = 1:length(r_list)
    u = 1/(4-power(4,1/3));
    U_permutation = speye(2^N);
    U_1a = expm(-1i*H1_shift*u*(t/r_list(i))/2) * expm(-1i*H2_shift*u*(t/r_list(i))/2) * expm(-1i*H3_shift*u*(t/r_list(i))) * expm(-1i*H2_shift*u*(t/r_list(i))/2) * expm(-1i*H1_shift*u*(t/r_list(i))/2);
    U_1b = expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))) * expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))/2);
    U_1 = U_1a * U_1a * U_1b * U_1a * U_1a;
    U_2a = expm(-1i*H1_shift*u*(t/r_list(i))/2) * expm(-1i*H3_shift*u*(t/r_list(i))/2) * expm(-1i*H2_shift*u*(t/r_list(i))) * expm(-1i*H3_shift*u*(t/r_list(i))/2) * expm(-1i*H1_shift*u*(t/r_list(i))/2);
    U_2b = expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))) * expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))/2);
    U_2 = U_2a * U_2a * U_2b * U_2a * U_2a;
    U_3a = expm(-1i*H2_shift*u*(t/r_list(i))/2) * expm(-1i*H1_shift*u*(t/r_list(i))/2) * expm(-1i*H3_shift*u*(t/r_list(i))) * expm(-1i*H1_shift*u*(t/r_list(i))/2) * expm(-1i*H2_shift*u*(t/r_list(i))/2);
    U_3b = expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))) * expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))/2);
    U_3 = U_3a * U_3a * U_3b * U_3a * U_3a;
    U_4a = expm(-1i*H2_shift*u*(t/r_list(i))/2) * expm(-1i*H3_shift*u*(t/r_list(i))/2) * expm(-1i*H1_shift*u*(t/r_list(i))) * expm(-1i*H3_shift*u*(t/r_list(i))/2) * expm(-1i*H2_shift*u*(t/r_list(i))/2);
    U_4b = expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))) * expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))/2);
    U_4 = U_4a * U_4a * U_4b * U_4a * U_4a;
    U_5a = expm(-1i*H3_shift*u*(t/r_list(i))/2) * expm(-1i*H2_shift*u*(t/r_list(i))/2) * expm(-1i*H1_shift*u*(t/r_list(i))) * expm(-1i*H2_shift*u*(t/r_list(i))/2) * expm(-1i*H3_shift*u*(t/r_list(i))/2);
    U_5b = expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))) * expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))/2);
    U_5 = U_5a * U_5a * U_5b * U_5a * U_5a;
    U_6a = expm(-1i*H3_shift*u*(t/r_list(i))/2) * expm(-1i*H1_shift*u*(t/r_list(i))/2) * expm(-1i*H2_shift*u*(t/r_list(i))) * expm(-1i*H1_shift*u*(t/r_list(i))/2) * expm(-1i*H3_shift*u*(t/r_list(i))/2);
    U_6b = expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H2_shift*(1-4*u)*(t/r_list(i))) * expm(-1i*H1_shift*(1-4*u)*(t/r_list(i))/2) * expm(-1i*H3_shift*(1-4*u)*(t/r_list(i))/2);
    U_6 = U_6a * U_6a * U_6b * U_6a * U_6a;
    U_list = {U_1, U_2, U_3, U_4, U_5, U_6};
    for j = 1:r_list(i)
        U_permutation = U_permutation * U_list{c_list{i}(j)};
    end
    Random_Trotter_Error_p_4_2{i} = U_permutation - expm(-1i * H_shift * t);
end