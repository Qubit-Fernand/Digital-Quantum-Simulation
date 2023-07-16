%% Set the energy threshold for the low-energy subspace
delta = 10;

energy = diag(energy);

% Projector onto low-energy subspace
projector = zeros(2^N);
for j = 1:2^N
    state = states(:,j);
    if energy(j) < delta
        projector = projector + state * state.';
    end
end

dt_list = [0.005 0.008 0.010 0.012 0.015 0.017 0.020];

r_50 = zeros(1,length(dt_list));
r_100 = zeros(1,length(dt_list));
r_200 = zeros(1,length(dt_list));
r_500 = zeros(1,length(dt_list));
r_50_Delta_5 = zeros(1,length(dt_list));
r_100_Delta_5 = zeros(1,length(dt_list));
r_200_Delta_5 = zeros(1,length(dt_list));
r_500_Delta_5 = zeros(1,length(dt_list));

for i = 1:length(dt_list)
    r_50(i) = norm(Trotter_Error_50{i});
    r_100(i) = norm(Trotter_Error_100{i});
    r_200(i) = norm(Trotter_Error_200{i});
    r_500(i) = norm(Trotter_Error_500{i});
    r_50_Delta_5(i) = norm(projector * Trotter_Error_50{i});
    r_100_Delta_5(i) = norm(projector * Trotter_Error_100{i});
    r_200_Delta_5(i) = norm(projector * Trotter_Error_200{i});
    r_500_Delta_5(i) = norm(projector * Trotter_Error_500{i});
end

%% save
save('/Users/AntiEntropy/Documents/Research/Quantum-Simulation/Numerical/Trotter_Suzuki/r_50.mat', 'r_50')
save('/Users/AntiEntropy/Documents/Research/Quantum-Simulation/Numerical/Trotter_Suzuki/r_100.mat', 'r_100')
save('/Users/AntiEntropy/Documents/Research/Quantum-Simulation/Numerical/Trotter_Suzuki/r_200.mat', 'r_200')
save('/Users/AntiEntropy/Documents/Research/Quantum-Simulation/Numerical/Trotter_Suzuki/r_500.mat', 'r_500')
save('/Users/AntiEntropy/Documents/Research/Quantum-Simulation/Numerical/Trotter_Suzuki/r_50_Delta_5.mat', 'r_50_Delta_5')
save('/Users/AntiEntropy/Documents/Research/Quantum-Simulation/Numerical/Trotter_Suzuki/r_100_Delta_5.mat', 'r_100_Delta_5')
save('/Users/AntiEntropy/Documents/Research/Quantum-Simulation/Numerical/Trotter_Suzuki/r_200_Delta_5.mat', 'r_200_Delta_5')
save('/Users/AntiEntropy/Documents/Research/Quantum-Simulation/Numerical/Trotter_Suzuki/r_500_Delta_5.mat', 'r_500_Delta_5')