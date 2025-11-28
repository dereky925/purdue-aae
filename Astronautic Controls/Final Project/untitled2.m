% Data
disturbance_km = 1:5;
sim_time_sec   = [38.270, 41.231, 49.068, 50.611, 52.589];

% Plot
figure('Color','white');
plot(disturbance_km, sim_time_sec, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Disturbance Magnitude [km]', 'FontSize', 20);
ylabel('Simulation Time [s]', 'FontSize', 20);
title('Simulation Time vs. Disturbance Magnitude', 'FontSize', 20);
set(gca, 'FontSize', 20);

