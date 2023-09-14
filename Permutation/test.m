%% Set the energy threshold for the low-energy subspace
global N;

delta = 14;

% Projector onto low-energy subspace
projector = zeros(2^N);
for j = 1:2^N
    state = states(:,j);
    if energy(j) <= delta
        projector = projector + state * state.';
    end
end

%% Calculate the norm of error norm
r_list = [10,50,100,200,500,1000];

p_2_1 = zeros(1,length(r_list));
p_2_2 = zeros(1,length(r_list));
p_4_1 = zeros(1,length(r_list));
p_4_2 = zeros(1,length(r_list));
p_2_1_Delta = zeros(1,length(r_list));
p_2_2_Delta = zeros(1,length(r_list));
p_4_1_Delta = zeros(1,length(r_list));
p_4_2_Delta = zeros(1,length(r_list));

for i = 1:length(r_list)
    p_2_1(i) = norm(Random_Trotter_Error_p_2_1{i});
    p_2_2(i) = norm(Random_Trotter_Error_p_2_2{i});
    p_4_1(i) = norm(Random_Trotter_Error_p_4_1{i});
    p_4_2(i) = norm(Random_Trotter_Error_p_4_2{i});
    p_2_1_Delta(i) = norm(projector * Random_Trotter_Error_p_2_1{i});
    p_2_2_Delta(i) = norm(projector * Random_Trotter_Error_p_2_2{i});
    p_4_1_Delta(i) = norm(projector * Random_Trotter_Error_p_4_1{i});
    p_4_2_Delta(i) = norm(projector * Random_Trotter_Error_p_4_2{i});
end