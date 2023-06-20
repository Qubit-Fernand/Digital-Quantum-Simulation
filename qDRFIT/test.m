dt = [0.005 0.008 0.010 0.012 0.015 0.017 0.020];

N_50 = zeros(1,length(dt));
N_100 = zeros(1,length(dt));
N_50_Delta_5 = zeros(1,length(dt));
N_100_Delta_5 = zeros(1,length(dt));

for i = 1:length(dt)
    N_50(i) = norm(qDRIFT_shift_Error_N_50{i});
    N_50_Delta_5(i) = norm(projector * qDRIFT_shift_Error_N_50{i});
end

save('/Users/AntiEntropy/Documents/Physics/Research/low-energy/Simulation/Figure/N_50.mat', 'N_50')
save('/Users/AntiEntropy/Documents/Physics/Research/low-energy/Simulation/Figure/N_50_Delta_5.mat', 'N_50_Delta_5')
save('/Users/AntiEntropy/Documents/Physics/Research/low-energy/Simulation/Figure/N_100.mat', 'N_100')
save('/Users/AntiEntropy/Documents/Physics/Research/low-energy/Simulation/Figure/N_100_Delta_5.mat', 'N_100_Delta_5')