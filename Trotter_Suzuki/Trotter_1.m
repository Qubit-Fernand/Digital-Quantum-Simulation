% 1-st order Trotter-Suzuki
function U1 = Trotter_1(dt)
global L;
global H1;
global H2;
global H3;
global H1_shift;
global H2_shift;
global H3_shift;

U1 = speye(2^(2*L));
U1 = U1 * expm(-1i * H1_shift * dt);
U1 = U1 * expm(-1i * H2_shift * dt);
U1 = U1 * expm(-1i * H3_shift * dt);

end