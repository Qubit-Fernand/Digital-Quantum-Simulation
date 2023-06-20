% p-th order Trotter-Suzuki: Recursion
function Up = Trotter_p(dt, r, p)
if p == 1
    Up = Trotter_1(dt, r);
elseif p == 2
    Up = Trotter_2(dt, r);
else
    u = 1/(4-4^(1/p-1));
    Up = Trotter_p(u * dt, r, p-2)^2 * Trotter_p((1-4*u)*dt, r,p-2) * Trotter_p(u * dt, r, p-2)^2;
end
end