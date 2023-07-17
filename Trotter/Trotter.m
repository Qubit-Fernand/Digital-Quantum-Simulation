% Initialize time step
r_list = [50 100 200 500];
dt_list = [0.005 0.008 0.010 0.012 0.015 0.017 0.020];

Trotter_Error_50 = cell(1,length(dt_list));
Trotter_Error_100 = cell(1,length(dt_list));
Trotter_Error_200 = cell(1,length(dt_list));
Trotter_Error_500 = cell(1,length(dt_list));

for i = 1:length(dt_list)
    U = Trotter_2(dt_list(i));
    Trotter_Error_50{i} = U^r_list(1) - expm(-1i * H_shift * dt_list(i) * 50);
    Trotter_Error_100{i} = U^r_list(2) - expm(-1i * H_shift * dt_list(i) * 100);
    Trotter_Error_200{i} = U^r_list(3) - expm(-1i * H_shift * dt_list(i) * 200);
    Trotter_Error_500{i} = U^r_list(4) - expm(-1i * H_shift * dt_list(i) * 500);
end 