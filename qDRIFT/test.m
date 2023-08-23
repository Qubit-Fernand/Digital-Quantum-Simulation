%% Set the energy threshold for the low-energy subspace
global N;

delta = 12;

% Projector onto low-energy subspace
projector = zeros(2^N);
for j = 1:2^N
    state = states(:,j);
    if energy(j) <= delta
        projector = projector + state * state.';
    end
end

%% Calculate the norm of the error
t_list = [1.0, 3.0, 5.0, 10.0];
num_list = [1000,10000,100000,200000,500000];

t_1 = zeros(1, length(num_list));
t_1_Delta = zeros(1, length(num_list));
t_3 = zeros(1, length(num_list));
t_3_Delta = zeros(1, length(num_list));
t_5 = zeros(1, length(num_list));
t_5_Delta = zeros(1, length(num_list));
t_10 = zeros(1, length(num_list));
t_10_Delta = zeros(1, length(num_list));

for i = 1:length(num_list)
    t_1(i) = norm(qDRIFT_Error_1{i});
    t_1_Delta(i) = norm(projector * qDRIFT_Error_1{i});
    t_3(i) = norm(qDRIFT_Error_3{i});
    t_3_Delta(i) = norm(projector * qDRIFT_Error_3{i});
    t_5(i) = norm(qDRIFT_Error_5{i});
    t_5_Delta(i) = norm(projector * qDRIFT_Error_5{i});
    t_10(i) = norm(qDRIFT_Error_10{i});
    t_10_Delta(i) = norm(projector * qDRIFT_Error_10{i});
end