%% Set the energy threshold for the low-energy subspace
delta = 5;

energy = diag(energy);

% Projector onto low-energy subspace
projector = zeros(2^(2*L));
for j = 1:2^(2*L)
    state = states(:,j);
    if energy(j) < delta
        projector = projector + state * state.';
    end
end

dt = [0.005 0.008 0.010 0.012 0.015 0.017 0.020];

r_50 = zeros(1,length(dt));
r_100 = zeros(1,length(dt));
r_200 = zeros(1,length(dt));
r_500 = zeros(1,length(dt));
r_50_Delta_5 = zeros(1,length(dt));
r_100_Delta_5 = zeros(1,length(dt));
r_200_Delta_5 = zeros(1,length(dt));
r_500_Delta_5 = zeros(1,length(dt));

for i = 1:length(dt)
    r_50(i) = norm(Error_r_50{i});
    r_100(i) = norm(Error_r_100{i});
    r_200(i) = norm(Error_r_200{i});
    r_500(i) = norm(Error_r_500{i});
    r_50_Delta_5(i) = norm(projector * Error_r_50{i});
    r_100_Delta_5(i) = norm(projector * Error_r_100{i});
    r_200_Delta_5(i) = norm(projector * Error_r_200{i});
    r_500_Delta_5(i) = norm(projector * Error_r_500{i});
end

%% save
save('/Users/AntiEntropy/Documents/Research/Quantum-Simulation/Numerical/Figure/r_50.mat', 'r_50')
save('/Users/AntiEntropy/Documents/Research/Quantum-Simulation/Numerical/Figure/r_100.mat', 'r_100')
save('/Users/AntiEntropy/Documents/Research/Quantum-Simulation/Numerical/Figure/r_200.mat', 'r_200')
save('/Users/AntiEntropy/Documents/Research/Quantum-Simulation/Numerical/Figure/r_500.mat', 'r_500')
save('/Users/AntiEntropy/Documents/Research/Quantum-Simulation/Numerical/Figure/r_50_Delta_5.mat', 'r_50_Delta_5')
save('/Users/AntiEntropy/Documents/Research/Quantum-Simulation/Numerical/Figure/r_100_Delta_5.mat', 'r_100_Delta_5')
save('/Users/AntiEntropy/Documents/Research/Quantum-Simulation/Numerical/Figure/r_200_Delta_5.mat', 'r_200_Delta_5')
save('/Users/AntiEntropy/Documents/Research/Quantum-Simulation/Numerical/Figure/r_500_Delta_5.mat', 'r_500_Delta_5')