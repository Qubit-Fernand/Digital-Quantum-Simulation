% Initialize time step
r = [50 100 200 500];
dt = [0.005 0.008 0.010 0.012 0.015 0.017 0.020];

Trotter_Error_50 = cell(1,length(dt));
Trotter_Error_100 = cell(1,length(dt));
Trotter_Error_200 = cell(1,length(dt));
Trotter_Error_500 = cell(1,length(dt));

for i = 1:length(dt)
    U = Trotter_2(dt(i));
    Trotter_Error_50{i} = U^r(1) - expm(-1i * H_shift * dt(i) * 50);
    Trotter_Error_100{i} = U^r(2) - expm(-1i * H_shift * dt(i) * 100);
    Trotter_Error_200{i} = U^r(3) - expm(-1i * H_shift * dt(i) * 200);
    Trotter_Error_500{i} = U^r(4) - expm(-1i * H_shift * dt(i) * 500);
end 