%% Initailize
delta = 2;
% energy = diag(energy);

% Projector onto low-energy subspace
projector = zeros(2^L);
for j = 1:2^L
    state = states(:,j);
    if energy(j) < delta
        projector = projector + state * state.';
    end
end

N = [1000,10000,100000,200000,500000,1000000];

t_10 = zeros(1, length(N));
t_10_Delta_5 = zeros(1, length(N));

%% Calculate Norm
for i = 1:length(N)
    t_10(i) = norm(qDRIFT_Error{i});
    t_10_Delta_5(i) = norm(projector * qDRIFT_Error{i});
end

%% Save Data
save('/Users/AntiEntropy/Documents/Research/Quantum-Simulation/Numerical/Figure/t_10.mat', 't_10')
save('/Users/AntiEntropy/Documents/Research/Quantum-Simulation/Numerical/Figure/t_10_Delta_5.mat', 't_10_Delta_5')

