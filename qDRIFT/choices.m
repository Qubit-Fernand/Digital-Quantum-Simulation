global H;
global H1;
global H2;
global H3;
global H_shift;
global H1_shift;
global H2_shift;
global H3_shift;

% Lanczos algorithm will return the eigenvalue with max norm
A = eigs(H1_shift, 1);
B = eigs(H2_shift, 1);
C = eigs(H3_shift, 1);

Repeat = 10;
t_list = [1.0, 3.0, 5.0, 10.0];
r_list = [10000,50000,100000,200000,500000,1000000];
c_list = cell(Repeat,length(r_list));

sum = A + B + C;
probabilities = [A/sum, B/sum, C/sum];

for i = 1:Repeat
    for j = 1:length(r_list)
        c_list{i, j} = zeros(1, r_list(j));
        for k = 1:r_list(j)
            choice = rand();
            if choice < probabilities(1)
                c_list{i, j}(k) = 1;
            elseif choice < probabilities(1) + probabilities(2)
                c_list{i, j}(k) = 2;
            else
                c_list{i, j}(k) = 3;
            end
        end
    end
end