% Initialize time step
r = [50 100 200 500];
dt = [0.005 0.008 0.010 0.012 0.015 0.017 0.020];

Error_50 = cell(1,length(dt));
Error_100 = cell(1,length(dt));
Error_200 = cell(1,length(dt));
Error_500 = cell(1,length(dt));

for i = 1:length(dt)
    Error_100{i} = Trotter_2(dt(i), 100) - Trotter_0(dt(i), 100);
    Error_200{i} = Trotter_2(dt(i), 200) - Trotter_0(dt(i), 200);
    Error_500{i} = Trotter_2(dt(i), 500) - Trotter_0(dt(i), 500);
end 

% plot(dt, Error_50, 'red-', dt, Error_100, 'blue-', dt, Error_150, 'green-', dt, Error_200, 'cyan-')