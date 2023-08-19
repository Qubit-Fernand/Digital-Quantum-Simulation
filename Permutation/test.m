%% Set the energy threshold for the low-energy subspace
global N;

delta = 4;

% Projector onto low-energy subspace
projector = zeros(2^N);
for j = 1:2^N
    state = states(:,j);
    if energy(j) <= delta
        projector = projector + state * state.';
    end
end

%% Calculate the norm of error norm
dt_list = [0.005 0.008 0.010 0.012 0.015 0.017 0.020];

r_50 = zeros(1,length(dt_list));
r_100 = zeros(1,length(dt_list));
r_200 = zeros(1,length(dt_list));
r_500 = zeros(1,length(dt_list));
r_50_Delta = zeros(1,length(dt_list));
r_100_Delta = zeros(1,length(dt_list));
r_200_Delta = zeros(1,length(dt_list));
r_500_Delta = zeros(1,length(dt_list));

for i = 1:length(dt_list)
    r_50(i) = norm(Random_Trotter_Error_50{i});
    r_100(i) = norm(Random_Trotter_Error_100{i});
    r_200(i) = norm(Random_Trotter_Error_200{i});
    r_500(i) = norm(Random_Trotter_Error_500{i});
    r_50_Delta(i) = norm(projector * Random_Trotter_Error_50{i});
    r_100_Delta(i) = norm(projector * Random_Trotter_Error_100{i});
    r_200_Delta(i) = norm(projector * Random_Trotter_Error_200{i});
    r_500_Delta(i) = norm(projector * Random_Trotter_Error_500{i});
end