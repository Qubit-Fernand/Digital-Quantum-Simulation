% 1-st order Trotter-Suzuki
function U1 = Trotter_1(dt, r)
global L;
global H_shift;

U1 = speye(2^(2*L));
for time_step = 1:r
    U1 = U1 * expm(-1i * H_shift * dt);
end
end