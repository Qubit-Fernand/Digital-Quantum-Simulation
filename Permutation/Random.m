function [a, b, c] = Random(a, b, c)
% Generate a random integer between 1 and 6 representing different permutations
random_order = randi([1, 6]);

% Define all possible permutations
permutations = {
    a, b, c;
    a, c, b;
    b, a, c;
    b, c, a;
    c, a, b;
    c, b, a
    };

% Select the permutation corresponding to the random integer
a = permutations{random_order, 1};
b = permutations{random_order, 2};
c = permutations{random_order, 3};
end