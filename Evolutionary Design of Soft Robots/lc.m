%delimeterIn='';
rand = importdata("random_search_lc.txt");%, delimeterIn);
disp(rand)
hill = importdata("hill_climber_lc.txt");%, delimeterIn);
ev = importdata("ev_search_lc.txt");%, delimeterIn);

hold on

plot(rand(:,1), rand(:,2), 'DisplayName', 'Random Search');
plot(hill(:,1), hill(:,2), 'DisplayName', 'Parallel Hillclimber');
plot(ev(:,1), ev(:,2), 'DisplayName', 'Evolutionary Algorithm');

legend('AutoUpdate', 'off');

xlabel('Evaluations');
ylabel('Speed: m/s');
title('Learning Curve');