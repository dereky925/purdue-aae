%% Simulate and plot increments of 2D Brownian motion for Problem Set 9(b)
clear;clc;close all

dt = 1e-3;              % time step
T = 1;                  % total time
N = T/dt;               % number of steps
M = 20;                 % number of sample paths

time = (0:dt:T-dt)';    % time vector for plotting

% Generate increments: each dw is N(0, dt)
dw = sqrt(dt) * randn(M, N);  % size M-by-N

% Compute 3-sigma bound for a single increment
sigma = sqrt(dt);
bound = 3 * sigma;

% Plot
figure('Color','w', 'Position', [100 100 1200 900]);
hold on;
plot(time, dw', 'LineWidth', 1);
plot(time, bound * ones(size(time)), 'k--', 'LineWidth', 2);
plot(time, -bound * ones(size(time)), 'k--', 'LineWidth', 2);
hold off;

title('Increments \Delta w_k of Standard Brownian Motion (M = 20)', 'FontSize', 20);
xlabel('Time t (s)', 'FontSize', 20);
ylabel('Increment \Delta w_k', 'FontSize', 20);
set(gca, 'FontSize', 20);
grid on;

%%
clear;clc;close all

dt = 1e-3;
T = 1;
N = T/dt;
time = 0:dt:T;
M = 20;

% System parameters
sigma = [2; 3];      % sigma1 = 2, sigma2 = 3
G = diag(sigma);

% Preallocate state trajectories: x(dim, time steps, sample index)
x = zeros(2, N+1, M);

% Simulate M sample paths
for m = 1:M
    % generate increments of Wiener process
    dw = sqrt(dt) * randn(2, N); % randn is delta w_k 
    % discrete update: x_{k+1} = x_k + G * dw_k
    for k = 1:N
        x(:, k+1, m) = x(:, k, m) + G * dw(:, k);
    end
end

% 3-sigma analytic bounds for each state component
bound1 = 3 * sigma(1) * sqrt(time);
bound2 = 3 * sigma(2) * sqrt(time);

% Plotting
figure('Color','w', 'Position',[200 200 1200 900]);

% Plot x1
subplot(2,1,1);
hold on;
plot(time, squeeze(x(1,:,:)), 'LineWidth', 1);
plot(time, bound1, 'k--', 'LineWidth', 2);
plot(time, -bound1, 'k--', 'LineWidth', 2);
hold off;
title('State x_1 Trajectories with \pm3\sigma Bounds','FontSize',20);
ylabel('x_1','FontSize',20);
set(gca,'FontSize',20);

grid on;

% Plot x2
subplot(2,1,2);
hold on;
plot(time, squeeze(x(2,:,:)), 'LineWidth', 1);
plot(time, bound2, 'k--', 'LineWidth', 2);
plot(time, -bound2, 'k--', 'LineWidth', 2);
hold off;
title('State x_2 Trajectories with \pm3\sigma Bounds','FontSize',20);
xlabel('Time t (s)','FontSize',20);
ylabel('x_2','FontSize',20);
set(gca,'FontSize',20);

grid on;

