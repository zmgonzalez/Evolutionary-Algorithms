delimeterIn=',';
data = importdata("variable_bot_energy.txt", delimeterIn);

x = data(1:22,1);
p = data(1:22,2);
k = data(1:22,3);
t = data(1:22,4);

hold on

plot(x, p, 'DisplayName', 'Potential Energy');
plot(x, k, 'DisplayName', 'Kinetic Energy');
plot(x, t, 'DisplayName', 'Total Energy');

legend('AutoUpdate', 'off');

xlabel('Time (s)');
ylabel('Energy (Joules, kg*(m/s)^2)');
title('Robot Energy');