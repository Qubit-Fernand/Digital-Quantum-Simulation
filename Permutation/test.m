%% Set the energy threshold for the low-energy subspace
global N;

delta = 8;

% Projector onto low-energy subspace
projector = zeros(2^N);
for j = 1:2^N
    state = states(:,j);
    if energy(j) <= delta
        projector = projector + state * state.';
    end
end

%% Calculate the norm of error norm
r_list = [1000 10000];

p_2_001 = zeros(1,length(r_list));
p_2_002 = zeros(1,length(r_list));
p_4_001 = zeros(1,length(r_list));
p_4_002 = zeros(1,length(r_list));
p_2_001_Delta = zeros(1,length(r_list));
p_2_002_Delta = zeros(1,length(r_list));
p_4_001_Delta = zeros(1,length(r_list));
p_4_002_Delta = zeros(1,length(r_list));

for i = 1:length(r_list)
    p_2_001(i) = norm(Random_Trotter_Error_2_001{i});
    p_2_002(i) = norm(Random_Trotter_Error_2_002{i});
    p_4_001(i) = norm(Random_Trotter_Error_4_001{i});
    p_4_002(i) = norm(Random_Trotter_Error_4_002{i});
    p_2_001_Delta(i) = norm(projector * Random_Trotter_Error_2_001{i});
    p_2_002_Delta(i) = norm(projector * Random_Trotter_Error_2_002{i});
    p_4_001_Delta(i) = norm(projector * Random_Trotter_Error_4_001{i});
    p_4_002_Delta(i) = norm(projector * Random_Trotter_Error_4_002{i});
end