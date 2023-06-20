% 2-nd order Trotter-Suzuki
function U2 = Trotter_2(dt, r)
global L;
global H1_shift;
global H2_shift;
global H3_shift;

U2 = speye(2^(2*L));
for step = 1:r
    U2 = U2 * expm(-1i * H1_shift * dt /2);
    U2 = U2 * expm(-1i * H2_shift * dt /2);
    U2 = U2 * expm(-1i * H3_shift * dt);
    U2 = U2 * expm(-1i * H2_shift * dt /2);
    U2 = U2 * expm(-1i * H1_shift * dt /2);
end
end