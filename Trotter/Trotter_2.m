% 2-nd order Trotter-Suzuki
function U2 = Trotter_2(dt)
global L;
global H1;
global H2;
global H3;
global H1_shift;
global H2_shift;
global H3_shift;

U2 = speye(2^(2*L));
U2 = U2 * expm(-1i * H1_shift * dt /2);
U2 = U2 * expm(-1i * H2_shift * dt /2);
U2 = U2 * expm(-1i * H3_shift * dt);
U2 = U2 * expm(-1i * H2_shift * dt /2);
U2 = U2 * expm(-1i * H1_shift * dt /2);
end
