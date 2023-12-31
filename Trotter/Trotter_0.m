% Time Evolution Operator
function U0 = Trotter_0(dt)
global N;
global H;
global H_shift;

U0 = expm(-1i * H_shift * dt);
end