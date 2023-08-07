%% Random Permutation 2nd-order Trotter-Suzuki
global N;
global H1;
global H2;
global H3;
global H1_shift;
global H2_shift;
global H3_shift;

% addpath('../Trotter');
r_list = [50 100 200 500];
dt_list = [0.005 0.008 0.010 0.012 0.015 0.017 0.020];

Random_Trotter_Error_50 = cell(1,length(dt_list));
Random_Trotter_Error_100 = cell(1,length(dt_list));
Random_Trotter_Error_200 = cell(1,length(dt_list));
Random_Trotter_Error_500 = cell(1,length(dt_list));

%% r = 50
r = 50;
for i = 1:length(dt_list)
    U_random = speye(2^N);
    dt = dt_list(i);
    for j = 1:r
        [H1_shift, H2_shift, H3_shift] = Random(H1_shift, H2_shift, H3_shift);
        U_random = U_random * Trotter_2(dt);
    end
    Random_Trotter_Error_50{i} = U_random - expm(-1i*(H1_shift+H2_shift+H3_shift)*dt*r);
end

%% r = 100
r = 100;
for i = 1:length(dt_list)
    U_random = speye(2^N);
    dt = dt_list(i);
    for j = 1:r
        [H1_shift, H2_shift, H3_shift] = Random(H1_shift, H2_shift, H3_shift);
        U_random = U_random * Trotter_2(dt);
    end
    Random_Trotter_Error_100{i} = U_random - expm(-1i*(H1_shift+H2_shift+H3_shift)*dt*r);
end

%% r = 200
r = 200;
for i = 1:length(dt_list)
    U_random = speye(2^N);
    dt = dt_list(i);
    for j = 1:r
        [H1_shift, H2_shift, H3_shift] = Random(H1_shift, H2_shift, H3_shift);
        U_random = U_random * Trotter_2(dt);
    end
    Random_Trotter_Error_200{i} = U_random - expm(-1i*(H1_shift+H2_shift+H3_shift)*dt*r);
end

%% r = 500
r = 500;
for i = 1:length(dt_list)
    U_random = speye(2^N);
    dt = dt_list(i);
    for j = 1:r
        [H1_shift, H2_shift, H3_shift] = Random(H1_shift, H2_shift, H3_shift);
        U_random = U_random * Trotter_2(dt);
    end
    Random_Trotter_Error_500{i} = U_random - expm(-1i*(H1_shift+H2_shift+H3_shift)*dt*r);
end