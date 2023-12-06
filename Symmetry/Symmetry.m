global N;
global H;
global H1;
global H2;
global H3;
global H_shift;
global H1_shift;
global H2_shift;
global H3_shift;

t = 1.0;
r_list = [10,50,100,200,500,1000];

Hadamard = hadamard(2)/sqrt(2);
for i = 2:N
        Hadamard = kron(Hadamard, hadamard(2)/sqrt(2));
end

Repeat = 10;
Raw_Error = cell(1,length(r_list));
SP_Error = cell(1,length(r_list));
Random_Error = cell(Repeat,length(r_list));

%% Raw;
for i = 1:length(r_list)
    U = Trotter_1(t/r_list(i));
    Raw_Error{i} = U^r_list(i) - expm(-1i * H_shift * t);
end

%% Optimal Symmetry Protection;
for i = 1:length(r_list)
    U = Trotter_1(t/r_list(i));
    U_symmetry = (Hadamard')^(r_list(i))*(U * Hadamard)^(r_list(i));
    SP_Error{i} = U_symmetry - expm(-1i * H_shift * t);
end

%% Random Symmetry Protection

for k = 1:Repeat
    for i = 1:length(r_list)
        U_symmetry = speye(2^N);
        U = Trotter_1(t/r_list(i));
        for r = 1:r_list(i)
            realPart1 = randn();
            imagPart1 = randn();
            realPart2 = randn();
            imagPart2 = randn();
            
            alpha = complex(realPart1, imagPart1);
            beta = complex(realPart2, imagPart2);
            
            % Scaling s.t. |alpha|^2 + |beta|^2 = 1
            normFactor = sqrt(abs(alpha)^2 + abs(beta)^2);
            alpha = alpha / normFactor;
            beta = beta / normFactor;
            
            C = [alpha, -conj(beta); beta, conj(alpha)];
            for j = 2:N
                C = kron(C, [alpha, -conj(beta); beta, conj(alpha)]);
            end
            U_symmetry = U_symmetry * (C'*U*C);
        end
        Random_Error{k, i} = U_symmetry - expm(-1i * H_shift * t);
    end
end