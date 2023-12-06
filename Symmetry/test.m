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

%% Calculate the norm of error norm
r_list = [10,50,100,200,500,1000];

% denote symmetry protection with postfix
Repeat = 10;
Raw = zeros(1,length(r_list));
SP = zeros(1,length(r_list));
Random = zeros(Repeat,length(r_list));
Raw_Delta = zeros(1,length(r_list));
SP_Delta = zeros(1,length(r_list));
Random_Delta = zeros(Repeat,length(r_list));

for i = 1:length(r_list)
    Raw(i) = norm(Raw_Error{i});
    SP(i) = norm(SP_Error{i});
    Raw_Delta(i) = norm(projector * Raw_Error{i});
    SP_Delta(i) = norm(projector * SP_Error{i});
    for j = 1:Repeat
        Random(j, i) = norm(Random_Error{j, i});
        Random_Delta(j, i) = norm(projector * Random_Error{j, i});
    end
end