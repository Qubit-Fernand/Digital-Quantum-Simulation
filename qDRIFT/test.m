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

%% Calculate the norm of the error
Repeat = 10;
t_list = [1.0, 3.0, 5.0, 10.0];
r_list = [10000,50000,100000,200000,500000,1000000];

t_1 = zeros(Repeat, length(r_list));
t_1_Delta = zeros(Repeat, length(r_list));
t_3 = zeros(Repeat, length(r_list));
t_3_Delta = zeros(Repeat, length(r_list));
t_5 = zeros(Repeat, length(r_list));
t_5_Delta = zeros(Repeat, length(r_list));
t_10 = zeros(Repeat, length(r_list));
t_10_Delta = zeros(Repeat, length(r_list));

for i = 1:Repeat
    for j = 1:length(r_list)
        t_1(i, j) = norm(qDRIFT_Error_1{i, j});
        t_1_Delta(i, j) = norm(projector * qDRIFT_Error_1{i, j});
        t_3(i, j) = norm(qDRIFT_Error_3{i, j});
        t_3_Delta(i, j) = norm(projector * qDRIFT_Error_3{i, j});
        t_5(i, j) = norm(qDRIFT_Error_5{i, j});
        t_5_Delta(i, j) = norm(projector * qDRIFT_Error_5{i, j});
        t_10(i, j) = norm(qDRIFT_Error_10{i, j});
        t_10_Delta(i, j) = norm(projector * qDRIFT_Error_10{i, j});
    end
end