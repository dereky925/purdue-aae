% =========================================================================
% 
% Filename:       HW1.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE590 - Fundamentals of State Estimation
% Professor:      Dr. Keith LeGrand
% Contact:        klegrand@purdue.edu
% Assignment:     HW 3
% Semester:       Fall 2025
% 
% Description:
%
% =========================================================================

%%

% Part a: Simulate constellation size and monthly budget with just process noise

clear; close all; clc; 
rng(1);  

K        = 36;          % months
p0       = 200;         % sats
b0       = 50;          % M$
sigma_b  = 1.5;         % M$, budget random-walk std

% Create p and b vectors
p = zeros(1, K+1);
b = zeros(1, K+1);
p(1) = p0;
b(1) = b0;

% Propogate
for k = 1:K
    % Budget random walk
    b(k+1) = b(k) + sigma_b*randn;

    % Constellation dynamics
    p(k+1) = 0.9*p(k) + b(k);
end

t = 0:K;

figure('Color','w','Position',[100 100 1400 800]);

subplot(1,2,1);
plot(t, p, 'LineWidth',1.8);
grid on; xlabel('Month'); ylabel('Satellites');
title('Constellation size p_k');
ax = gca;
ax.FontSize = 16;

subplot(1,2,2);
plot(t, b, 'LineWidth',1.8);
grid on; xlabel('Month'); ylabel('Budget (M$)');
title('Monthly budget b_k');

sgtitle('Derek Yu - Part (a): Simulation of Constellation and Budget','FontSize',25);
ax = gca;
ax.FontSize = 16;

%% Part c

clear; close all; clc; 
rng(1);              

% system dynamics
A = [0.9 1.0;   
     0.0 1.0];  

sigma_b = 1.5;              % $ 1 sig uncertainty
Qw = sigma_b^2;             % process noise 
G  = [0; 1];                % process noise for only b_k
Q  = G*Qw*G.';              % state space process noise covariance

H  = [1 0];                 % mapping, measure constellation size only
R  = 7^2;                   % meas noise satellites^2

K  = 36;  % months

% Init
x_true = zeros(2, K+1);
x_true(:,1) = [200; 50];    % p0=200 sats, b0=50 M$

% Generate a true trajectory and measurements
w = [zeros(1,K); sigma_b*randn(1,K)];   % only budget gets noise
z = zeros(1,K);                         % measurements z_k

for k = 1:K
    x_true(:,k+1) = A*x_true(:,k) + G*(sigma_b*randn);  % same as adding w(:,k)
    z(k) = H*x_true(:,k+1) + 7*randn;                   % measurement at month k
end

p_true = x_true(1,2:end);
b_true = x_true(2,2:end);

%  Kalman filter setup 
xhat = zeros(2, K+1);
xhat(:,1) = [170; 80];                 % init estimate

P = zeros(2,2,K+1);
P(:,:,1) = diag([50^2, 20^2]);         % init covariance

% storage for posterior (after-update) values
xhat_post = zeros(2,K);
P_post    = zeros(2,2,K);

for k = 1:K
    % Predict
    x_pred = A*xhat(:,k);
    P_pred = A*P(:,:,k)*A.' + Q;

    % Update
    S  = H*P_pred*H.' + R;        % innovation covariance
    Kk = (P_pred*H.')/S;          % Kalman gain (2x1)

    innov = z(k) - H*x_pred;      % residual
    x_upd = x_pred + Kk*innov;
    P_upd = (eye(2) - Kk*H)*P_pred;

    % Store posterior
    xhat_post(:,k) = x_upd;
    P_post(:,:,k)  = P_upd;

    % Roll to next step
    xhat(:,k+1) = x_upd;
    P(:,:,k+1)  = P_upd;
end

p_hat = xhat_post(1,:);
b_hat = xhat_post(2,:);
sig_p = sqrt(squeeze(P_post(1,1,:))).';
sig_b = sqrt(squeeze(P_post(2,2,:))).';

t = 1:K;


figure('Color','w','Position',[100 100 1400 900]);

subplot(2,2,1);
plot(t, p_true, 'LineWidth',2); hold on;
plot(t, p_hat,  '--', 'LineWidth',2);
grid on; xlabel('Month'); ylabel('Satellites');
title('i) Constellation size: true vs posterior estimate');
legend('True p_k','Estimated p_k','Location','best');
ax = gca;
ax.FontSize = 16;

subplot(2,2,2);
plot(t, b_true, 'LineWidth',2); hold on;
plot(t, b_hat,  '--', 'LineWidth',2);
grid on; xlabel('Month'); ylabel('Budget (M$)');
title('ii) Monthly budget: true vs posterior estimate');
legend('True b_k','Estimated b_k','Location','best');
ax = gca;
ax.FontSize = 16;

subplot(2,2,3);
err_p = p_true - p_hat;
plot(t, err_p, 'LineWidth',2); hold on;
plot(t,  sig_p, ':', 'LineWidth',3);
plot(t, -sig_p, ':', 'LineWidth',3);
grid on; xlabel('Month'); ylabel('Satellites');
title('iii) p-state error with \pm1\sigma bounds');
legend('p_k - p_{est,k}','+\sigma_p','-\sigma_p','Location','best');
ax = gca;
ax.FontSize = 16;

subplot(2,2,4);
err_b = b_true - b_hat;
plot(t, err_b, 'LineWidth',2); hold on;
plot(t,  sig_b, ':', 'LineWidth',3);
plot(t, -sig_b, ':', 'LineWidth',3);
grid on; xlabel('Month'); ylabel('M$');
title('iv) b-state error with \pm1\sigma bounds');
legend('b_k - b_{est,k}','+\sigma_b','-\sigma_b','Location','best');
ax = gca;
ax.FontSize = 16;

sgtitle('Kalman Filter for Constellation Size and Budget (36 months)');


%% Problem 2, Orbit EKF

clear; close all; clc;

% Given
mu = 3.986004418e5;           % [km^3/s^2]
Qr = 1e-8; Qth = 1e-16;       % process noise covariances
Qc = diag([Qr, Qth]);          % continuous-time process noise
R = 1e-4;                 % measurement noise variance (rad^2)

% True initial conditions
x_true0 = [6753.137; 0; 0; 0.00113975]; % [r, rdot, theta, thetadot]

% Estimated initial mean and covariance
x_hat = [6.75681367e3;
         -1.37015864e-1;
          4.45127583e-2;
          1.08593276e-3];
P = diag([40, 1e-2, 1e-3, 1e-8]);

% Orbital period
a = x_true0(1);                      % roughly semi-major axis
P_orbit = 2*pi*sqrt(a^3/mu);         % seconds

% Simulation parameters
dt_meas = 30;                        % measurement period [s]
dt_int = 1;                          % integration time step
tspan = 0:dt_int:P_orbit;

% Noise matrices
G = [0 0; 1 0; 0 0; 0 1];

% Storage
nx = length(x_true0);
x_true = zeros(nx, length(tspan));
x_hat_store = zeros(nx, length(tspan));
P_store = zeros(nx, nx, length(tspan));

x_true(:,1) = x_true0;
x_hat_store(:,1) = x_hat;
P_store(:,:,1) = P;

% Integrate true dynamics with process noise
for k = 2:length(tspan)
    t = tspan(k-1);
    x = x_true(:,k-1);
    
    % w = [sqrt(Qr)*randn; sqrt(Qth)*randn];   % process noise
    % dx = f_orbit(x, mu) + G*w;
    % dx = f_orbit(x, mu)*w;
    
    dx = f_orbit(x, mu);

    x_true(:,k) = x + dx*dt_int;
    x_true(3,k) = wrapTo2Pi(x_true(3,k));    % wrap angle
end

% Generate noisy measurements (every 30 s)
meas_times = 0:dt_meas:P_orbit;
z_meas = interp1(tspan, x_true(3,:), meas_times) + sqrt(R)*randn(size(meas_times));

% EKF
k_meas = 2;
for k = 2:length(tspan)

    % Predict
    x_pred = x_hat + f_orbit(x_hat, mu)*dt_int;          % Euler predict
    F = F_jacobian(x_hat, mu);                           % linearize and evaluate
    Qd = G*Qc*G.';                            
    P_pred = P + (F*P + P*F.' + Qd)*dt_int;              % simple Euler approx
    
    % Measurement update
    if abs(tspan(k) - meas_times(k_meas)) < dt_int/2  % Check if within time interval to get measurement

        H = [0 0 1 0];
        z_pred = x_pred(3);
        y = z_meas(k_meas) - z_pred;  % innovation
        S = H*P_pred*H.' + R;
        K = P_pred*H.'/S; % Kalman gain
        x_hat = x_pred + K*y; % Update estimate
        P = (eye(nx) - K*H)*P_pred;
        x_hat(3) = wrapTo2Pi(x_hat(3));
        k_meas = min(k_meas+1, length(meas_times));

    else % Otherwise, just propagate
        x_hat = x_pred; 
        P = P_pred;
    end

    % store
    x_hat_store(:,k) = x_hat;
    P_store(:,:,k) = P;
end


% Extract
r_true     = x_true(1,:);
theta_true = x_true(3,:);
r_hat      = x_hat_store(1,:);
theta_hat  = x_hat_store(3,:);
t          = tspan;

% Plot true vs estimated 2D trajectory
figure('Color','w','Position',[100 100 1000 900]);
plot(r_true.*cos(theta_true), r_true.*sin(theta_true), 'b', 'LineWidth', 2); hold on;
plot(r_hat.*cos(theta_hat), r_hat.*sin(theta_hat), 'r--', 'LineWidth', 2);
axis equal; grid minor; box on;
xlabel('x [km]', 'FontSize', 16);
ylabel('y [km]', 'FontSize', 16);
title('True vs EKF-Estimated Orbit', 'FontSize', 18, 'FontWeight', 'bold');
legend('True orbit','Estimated orbit','FontSize',14,'Location','best');
set(gca, 'FontSize', 14, 'LineWidth', 1.2);

% Estimation errors with 3 sigma bounds
state_names = {'r [km]','ṙ [km/s]','θ [rad]','θ̇ [rad/s]'};
state =  {'r','ṙ','θ','θ̇'};
figure('Color','w','Position',[1000 50 1200 900]);

for i = 1:4
    subplot(2,2,i);
    err   = x_true(i,:) - x_hat_store(i,:);
    sigma = sqrt(squeeze(P_store(i,i,:)))';
    
    plot(t, err, 'b', 'LineWidth', 1.8); hold on;
    plot(t,  3*sigma, 'r--', 'LineWidth', 1.6);
    plot(t, -3*sigma, 'r--', 'LineWidth', 1.6);
    
    grid minor; box on;
    xlabel('Time [s]', 'FontSize', 14);
    ylabel(state_names{i}, 'FontSize', 14);
    title([ state{i} ' Estimate Error ±3σ'], 'FontSize', 16, 'FontWeight', 'bold');
    set(gca, 'FontSize', 14, 'LineWidth', 1.2);
    legend('Error','+3σ','−3σ','FontSize',12,'Location','best');
end

sgtitle('EKF Estimation Errors and ±3σ Bounds', 'FontSize', 20, 'FontWeight', 'bold');


% Helper functions

% Orbit dynamics
function dx = f_orbit(x, mu)
    r = x(1); rd = x(2); th = x(3); thd = x(4);
    dx = zeros(4,1);
    dx(1) = rd;
    dx(2) = r*thd^2 - mu/r^2;
    dx(3) = thd;
    dx(4) = -2*rd*thd/r;
end

% Dynaqmics Jacobian
function F = F_jacobian(x, mu)
    r = x(1); rd = x(2); th = x(3); thd = x(4);
    F = [0 1 0 0;
         thd^2 + 2*mu/r^3 0 0 2*r*thd;
         0 0 0 1;
         2*rd*thd/r^2 -2*thd/r 0 -2*rd/r];
end



