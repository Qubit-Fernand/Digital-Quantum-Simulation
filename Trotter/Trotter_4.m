% 4-th order Trotter-Suzuki
function U4 = Trotter_4(dt)
global N;
global H1;
global H2;
global H3;
global H1_shift;
global H2_shift;
global H3_shift;

U2 = speye(2^N);
U2 = U2 * expm(-1i * H1_shift * dt /2);
U2 = U2 * expm(-1i * H2_shift * dt /2);
U2 = U2 * expm(-1i * H3_shift * dt);
U2 = U2 * expm(-1i * H2_shift * dt /2);
U2 = U2 * expm(-1i * H1_shift * dt /2);

u = 1/(4-power(4,1/3));
U4 = Trotter_2(u*dt)*Trotter_2(u*dt)*Trotter_2((1-4*u)*dt)*Trotter_2(u*dt)*Trotter_2(u*dt);
end
