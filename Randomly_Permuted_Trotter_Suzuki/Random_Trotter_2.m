% randomly permuted 2-nd order Trotter-Suzuki
function U_random = Random_Trotter_2(dt, r)
global L
global H1
global H2
global H3

U_random_1 = speye(2^(2*L));
for step = 1:r
    U_random_1 = U_random_1 * expm(-1i * H1 * dt /2);
    U_random_1 = U_random_1 * expm(-1i * H2 * dt /2);
    U_random_1 = U_random_1 * expm(-1i * H3 * dt);
    U_random_1 = U_random_1 * expm(-1i * H2 * dt /2);
    U_random_1 = U_random_1 * expm(-1i * H1 * dt /2);
end

U_random_2 = speye(2^(2*L));
for step = 1:r
    U_random_2 = U_random_2 * expm(-1i * H1 * dt /2);
    U_random_2 = U_random_2 * expm(-1i * H3 * dt /2);
    U_random_2 = U_random_2 * expm(-1i * H2 * dt);
    U_random_2 = U_random_2 * expm(-1i * H3 * dt /2);
    U_random_2 = U_random_2 * expm(-1i * H1 * dt /2);
end

U_random_3 = speye(2^(2*L));
for step = 1:r
    U_random_3 = U_random_3 * expm(-1i * H2 * dt /2);
    U_random_3 = U_random_3 * expm(-1i * H1 * dt /2);
    U_random_3 = U_random_3 * expm(-1i * H3 * dt);
    U_random_3 = U_random_3 * expm(-1i * H1 * dt /2);
    U_random_3 = U_random_3 * expm(-1i * H2 * dt /2);
end

U_random_4 = speye(2^(2*L));
for step = 1:r
    U_random_4 = U_random_4 * expm(-1i * H2 * dt /2);
    U_random_4 = U_random_4 * expm(-1i * H3 * dt /2);
    U_random_4 = U_random_4 * expm(-1i * H1 * dt);
    U_random_4 = U_random_4 * expm(-1i * H3 * dt /2);
    U_random_4 = U_random_4 * expm(-1i * H2 * dt /2);
end

U_random_5 = speye(2^(2*L));
for step = 1:r
    U_random_5 = U_random_5 * expm(-1i * H3 * dt /2);
    U_random_5 = U_random_5 * expm(-1i * H1 * dt /2);
    U_random_5 = U_random_5 * expm(-1i * H2 * dt);
    U_random_5 = U_random_5 * expm(-1i * H1 * dt /2);
    U_random_5 = U_random_5 * expm(-1i * H3 * dt /2);
end

U_random_6 = speye(2^(2*L));
for step = 1:r
    U_random_6 = U_random_6 * expm(-1i * H3 * dt /2);
    U_random_6 = U_random_6 * expm(-1i * H2 * dt /2);
    U_random_6 = U_random_6 * expm(-1i * H1 * dt);
    U_random_6 = U_random_6 * expm(-1i * H2 * dt /2);
    U_random_6 = U_random_6 * expm(-1i * H3 * dt /2);
end

U_random = (U_random_1 + U_random_2 + U_random_3 + U_random_4 + U_random_5 + U_random_6)/6;
end