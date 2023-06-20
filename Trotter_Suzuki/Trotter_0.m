% Time Evolution Operator
function U0 = Trotter_0(dt, r)
global H_shift;

U0 = expm(-1i * H_shift * r * dt);
end