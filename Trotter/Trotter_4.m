% 4-th order Trotter-Suzuki
function U4 = Trotter_4(dt)
global N;
global H1;
global H2;
global H3;
global H1_shift;
global H2_shift;
global H3_shift;

u = 1/(4-power(4,1/3));
U4 = Trotter_2(u*dt)*Trotter_2(u*dt)*Trotter_2((1-4*u)*dt)*Trotter_2(u*dt)*Trotter_2(u*dt);
end
