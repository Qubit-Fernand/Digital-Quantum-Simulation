% p-th order Trotter-Suzuki: Recursion
function U_p = Trotter_p(dt, p)
global N;
global H1;
global H2;
global H3;
global H1_shift;
global H2_shift;
global H3_shift;

if p == 1
    U_p = Trotter_1(dt);
elseif p == 2
    U_p = Trotter_2(dt);
else
    u = 1/(4-4^(1/p-1));
    U_p = Trotter_p(u * dt, p-2)^2 * Trotter_p((1-4*u)*dt,p-2) * Trotter_p(u * dt, p-2)^2;
end
end