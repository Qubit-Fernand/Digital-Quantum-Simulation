% Set the parameters
N = 4; % Number of lattice sites
alpha = 2; % Exponent

% Initialize the Hamiltonian matrices and coordinate matrix
H_1 = sparse(2^N, 2^N); % 左上-右下
H_2 = sparse(2^N, 2^N); % 水平-方向
H_3 = sparse(2^N, 2^N); % 右上-左下
coords = zeros(N, 2);

% Construct the Hamiltonian parts
for i = 1:N
    coords(i,:) = [ceil(i/N), mod(i-1,N)+1]; % Calculate the coordinates of the i-th qubit on the lattice
end
for i = 1:N
    for j = 1:N
        if i < j && (coords(i,1) < coords(j,1) || (coords(i,1) == coords(j,1) && coords(i,2) > coords(j,2))) % Check if i and j form a left-down pair
            dist = norm(coords(i,:)-coords(j,:))^alpha; % Calculate the Euclidean distance to the power of alpha
            Hij = -1/dist;
            op_i = speye(2^(i-1));
            op_j = speye(2^(j-i-1));
            op_k = speye(2^(N-j));
            H_ij = kron(kron(op_i,op_j),kron(Z,Hij)); % Use Kronecker product to construct H_{DL,ij}
            row = find(H_ij); % Find the vector indices of non-zero elements in H_{DL,ij}
            col = row+2^(i-1)+2^(j-1); % Calculate the vector indices of non-zero elements in H_{DL,ij}
            H_1(row, col) = H_ij(H_ij~=0); % Assign the non-zero elements of H_{DL,ij} to the corresponding positions
        elseif j == i+1 % Check if i and j form a horizontal pair
            dist = norm(coords(i,:)-coords(j,:))^alpha; % Calculate the Euclidean distance to the power of alpha
            Hij = -1/dist;
            op_i = speye(2^(i-1));
            op_j = speye(2^(j-i-1));
            op_k = speye(2^(N-j));
            H_ij = kron(kron(op_i,op_j),kron(Z,Hij)); % Use Kronecker product to construct H_{H,ij}
            row = find(H_ij); % Find the vector indices of non-zero elements in H_{H,ij}
            col = row+2^(i-1)+2^(j-1); % Calculate the vector indices of non-zero elements in H_{H,ij};
            H_2(row, col) = H_ij(H_ij~=0); % Assign the non-zero elements of H_{H,ij} to the corresponding positions
        elseif i < j && (coords(i,1) > coords(j,1) || (coords(i,1) == coords(j,1) && coords(i,2) < coords(j,2))) % Check if i and j form a left-up pair
            dist = norm(coords(i,:)-coords(j,:))^alpha; % Calculate the Euclidean distance to the power of alpha
            Hij = -1/dist;
            op_i = speye(2^(i-1));
            op_j = speye(2^(j-i-1));
            op_k = speye(2^(N-j));
            H_ij = kron(kron(op_i,op_j),kron(Z,Hij)); % Use Kronecker product to construct H_{UL,ij}
            row = find(H_ij); % Find the vector indices of non-zero elements in H_{UL,ij}
            col = row+2^(i-1)+2^(j-1); % Calculate the vector indices of non-zero elements in H_{UL,ij}
            H_3(row, col) = H_ij(H_ij~=0); % Assign the non-zero elements of H_{UL,ij} to the corresponding positions
        end
    end
end

% Shift the Hamiltonian by the ground state energy
energy1 = eigs(H1);
energy2 = eigs(H2);
energy3 = eigs(H3);
ground_energy1 = min(energy1);
ground_energy2 = min(energy2);
ground_energy3 = min(energy3);
H1_shift = H1 - floor(ground_energy1) * speye(2^N);
H2_shift = H2 - floor(ground_energy2) * speye(2^N);
H3_shift = H3 - floor(ground_energy3) * speye(2^N);

H = H1_shift+H2_shift+H3_shift;