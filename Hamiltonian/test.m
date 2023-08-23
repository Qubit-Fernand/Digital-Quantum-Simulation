%% Act on a Random Generated State
% Set seed to ensure the same initial state is outputted each time
rng(0); % Set seed value to 0
v = randn(2^N, 1) + 1i * randn(2^N, 1);
rng('default'); % Restore default seed value

% Initialize quantum register
psi = v/norm(v);

% Project onto low-energy subspace
psi = projector * psi;