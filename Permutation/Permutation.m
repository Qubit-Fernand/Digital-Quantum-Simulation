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

addpath('Trotter');
r_list = [50 100 200 500];
dt_list = [0.010 0.020];

Random_Trotter_Error_2_001 = cell(1,length(r_list));
Random_Trotter_Error_2_002 = cell(1,length(r_list));
Random_Trotter_Error_4_001 = cell(1,length(r_list));
Random_Trotter_Error_4_002 = cell(1,length(r_list));

%% dt = 0.01 p = 2;
dt = 0.01;
for i = 1:length(r_list)
    r = r_list(i);
    for j = 1:r
        U_random = speye(2^N);
        [H1_shift, H2_shift, H3_shift] = Random(H1_shift, H2_shift, H3_shift);
        U_random = U_random * Trotter_2(dt);
    end
    Random_Trotter_Error_2_001{i} = U_random - expm(-1i * H_shift * dt * r);
end

%% dt = 0.01 p = 4;
for i = 1:length(r_list)
    r = r_list(i);
    for j = 1:r
        U_random = speye(2^N);
        [H1_shift, H2_shift, H3_shift] = Random(H1_shift, H2_shift, H3_shift);
        U_random = U_random * Trotter_4(dt);
    end
    Random_Trotter_Error_4_001{i} = U_random - expm(-1i * H_shift * dt * r);
end

%% dt = 0.02 p = 2;
dt = 0.02;
for i = 1:length(r_list)
    r = r_list(i);
    for j = 1:r
        U_random = speye(2^N);
        [H1_shift, H2_shift, H3_shift] = Random(H1_shift, H2_shift, H3_shift);
        U_random = U_random * Trotter_2(dt);
    end
    Random_Trotter_Error_2_002{i} = U_random - expm(-1i * H_shift * dt * r);
end

%% dt = 0.04 p = 4
for i = 1:length(r_list)
    r = r_list(i);
    for j = 1:r
        U_random = speye(2^N);
        [H1_shift, H2_shift, H3_shift] = Random(H1_shift, H2_shift, H3_shift);
        U_random = U_random * Trotter_4(dt);
    end
    Random_Trotter_Error_4_002{i} = U_random - expm(-1i * H_shift * dt * r);
end