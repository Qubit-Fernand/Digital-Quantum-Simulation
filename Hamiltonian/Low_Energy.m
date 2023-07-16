%% Set the energy threshold for the low-energy subspace
delta = 5;
% energy = diag(energy);

% Projector onto low-energy subspace
projector = zeros(2^N);
for j = 1:2^N
    state = states(:,j);
    if energy(j) < delta
        projector = projector + state * state.';
    end
end

%% Act on a Random Generated State
% Set seed to ensure the same initial state is outputted each time
rng(0); % Set seed value to 0
v = randn(2^(2*L), 1) + 1i * randn(2^(2*L), 1);
rng('default'); % Restore default seed value

% Initialize quantum register
psi = v/norm(v);

% Project onto low-energy subspace
psi = projector * psi;