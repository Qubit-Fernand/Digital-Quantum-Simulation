% randomly permuted 2-nd order Trotter-Suzuki
function U_random = Random_Trotter_2(dt)
global L;
global H1;
global H2;
global H3;
global H1_shift;
global H2_shift;
global H3_shift;

U_random_1 = speye(2^(2*L));
U_random_1 = U_random_1 * expm(-1i * H1_shift * dt /2);
U_random_1 = U_random_1 * expm(-1i * H2_shift * dt /2);
U_random_1 = U_random_1 * expm(-1i * H3_shift * dt);
U_random_1 = U_random_1 * expm(-1i * H2_shift * dt /2);
U_random_1 = U_random_1 * expm(-1i * H1_shift * dt /2);

U_random_2 = speye(2^(2*L));
U_random_2 = U_random_2 * expm(-1i * H2_shift * dt /2);
U_random_2 = U_random_2 * expm(-1i * H1_shift * dt /2);
U_random_2 = U_random_2 * expm(-1i * H3_shift * dt);
U_random_2 = U_random_2 * expm(-1i * H1_shift * dt /2);
U_random_2 = U_random_2 * expm(-1i * H2_shift * dt /2);

U_random_3 = speye(2^(2*L));
U_random_3 = U_random_3 * expm(-1i * H1_shift * dt /2);
U_random_3 = U_random_3 * expm(-1i * H3_shift * dt /2);
U_random_3 = U_random_3 * expm(-1i * H2_shift * dt);
U_random_3 = U_random_3 * expm(-1i * H3_shift * dt /2);
U_random_3 = U_random_3 * expm(-1i * H1_shift * dt /2);

U_random_4 = speye(2^(2*L));
U_random_4 = U_random_4 * expm(-1i * H3_shift * dt /2);
U_random_4 = U_random_4 * expm(-1i * H1_shift * dt /2);
U_random_4 = U_random_4 * expm(-1i * H2_shift * dt);
U_random_4 = U_random_4 * expm(-1i * H1_shift * dt /2);
U_random_4 = U_random_4 * expm(-1i * H3_shift * dt /2);

U_random_5 = speye(2^(2*L));
U_random_5 = U_random_5 * expm(-1i * H2_shift * dt /2);
U_random_5 = U_random_5 * expm(-1i * H3_shift * dt /2);
U_random_5 = U_random_5 * expm(-1i * H1_shift * dt);
U_random_5 = U_random_5 * expm(-1i * H3_shift * dt /2);
U_random_5 = U_random_5 * expm(-1i * H2_shift * dt /2);

U_random_6 = speye(2^(2*L));
U_random_6 = U_random_6 * expm(-1i * H3_shift * dt /2);
U_random_6 = U_random_6 * expm(-1i * H2_shift * dt /2);
U_random_6 = U_random_6 * expm(-1i * H1_shift * dt);
U_random_6 = U_random_6 * expm(-1i * H2_shift * dt /2);
U_random_6 = U_random_6 * expm(-1i * H3_shift * dt /2);

U_random = (U_random_1 + U_random_2 + U_random_3 + U_random_4 + U_random_5 + U_random_6)/6;
end