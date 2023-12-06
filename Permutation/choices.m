global H;
global H1;
global H2;
global H3;
global H_shift;
global H1_shift;
global H2_shift;
global H3_shift;

Repeat = 10;
t_list = [1.0 2.0];
r_list = [10,50,100,200,500,1000];
c_list = cell(Repeat,length(r_list));

for i = 1:Repeat
    for j = 1:length(r_list)
        c_list{i, j} = zeros(1, r_list(j));
        for k = 1:r_list(j)
            choice = randi([1, 6]);
            c_list{i, j}(k) = choice;
        end
    end
end