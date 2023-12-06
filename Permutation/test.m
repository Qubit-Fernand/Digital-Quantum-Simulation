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
Repeat = 10;
r_list = [10,50,100,200,500,1000];

p_2_1 = zeros(Repeat,length(r_list));
p_2_2 = zeros(Repeat,length(r_list));
p_4_1 = zeros(Repeat,length(r_list));
p_4_2 = zeros(Repeat,length(r_list));
p_2_1_Delta = zeros(Repeat,length(r_list));
p_2_2_Delta = zeros(Repeat,length(r_list));
p_4_1_Delta = zeros(Repeat,length(r_list));
p_4_2_Delta = zeros(Repeat,length(r_list));

for i = 1:Repeat
    for j = 1:length(r_list)
        p_2_1(i, j) = norm(Random_Trotter_Error_p_2_1{i, j});
        p_2_2(i, j) = norm(Random_Trotter_Error_p_2_2{i, j});
        p_4_1(i, j) = norm(Random_Trotter_Error_p_4_1{i, j});
        p_4_2(i, j) = norm(Random_Trotter_Error_p_4_2{i, j});
        p_2_1_Delta(i, j) = norm(projector * Random_Trotter_Error_p_2_1{i, j});
        p_2_2_Delta(i, j) = norm(projector * Random_Trotter_Error_p_2_2{i, j});
        p_4_1_Delta(i, j) = norm(projector * Random_Trotter_Error_p_4_1{i, j});
        p_4_2_Delta(i, j) = norm(projector * Random_Trotter_Error_p_4_2{i, j});
    end
end