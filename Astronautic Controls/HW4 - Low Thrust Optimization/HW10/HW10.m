% =========================================================================
% 
% Filename:       HW10.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE590 - Applied Control in Astronautics
% Professor:      Dr. Kenshiro Oguri
% Assignment:     HW 10
% Semester:       Spring 2025
% 
% Description: Homework 10
%
% =========================================================================

%%
clear; clc; close all;
format long

% PARAMETERS
mu = 1;              % Non-dimensional gravitational parameter
umax = 0.1;          % Maximum control magnitude
tf = 8;            % Fixed final time (non-dimensional)

% Smoothing parameters for continuation
rhos = [1e-3];

% Earth initial conditions (circular orbit at 1 AU)
% State vector: [position (2); velocity (2)]
x0 = [1; 0; 0; 1];

% Mars orbit parameters (circular orbit)
a_M = 1.524;         % Mars semi-major axis (AU)
nu_M0 = pi;          % Mars initial true anomaly (radians)
omega_M = sqrt(mu/(a_M^3));  % Mars angular rate

% Anonymous function for Mars' state (position and velocity) at time t
marsState = @(t) [ a_M*cos(nu_M0 + omega_M*t);
                   a_M*sin(nu_M0 + omega_M*t);
                  -a_M*omega_M*sin(nu_M0 + omega_M*t);
                   a_M*omega_M*cos(nu_M0 + omega_M*t) ];

% Terminal condition: spacecraft must rendezvous with Mars at tf.
x_target = marsState(tf);

% Plot: Approximated Thrust vs. Switching Function
% The switching function is S = ||lambda_v|| - 1, and the smoothed thrust is:
% Gamma = umax/2*(1 + tanh((||lambda_v|| - 1)/rho))
s_vec = linspace(-1,1,200);
figure('Color','white','Position',[2000 0 1000 800]); hold on;
colors = lines(length(rhos));
for i = 1:length(rhos)
    rho_val = rhos(i);
    Gamma = umax/2*(1 + tanh(s_vec./rho_val));
    plot(s_vec, Gamma, 'Color', colors(i,:), 'LineWidth', 3, 'DisplayName', sprintf('\\rho = %g', rho_val));
end
xlabel('Switching Function S = ||\lambda_v|| - 1');
ylabel('Approximated Thrust Magnitude \Gamma*');
title('Thrust Profile vs. Switching Function');
legend('Location','Best'); grid on;
ax = gca; ax.FontSize = 20;

% Solve the Two-Point Boundary Value Problem via Shooting
% The shooting function adjusts the initial costate lambda0 so that x(tf) = Mars' state.
lambda0_guess = [0.8; 0.1; 0.2; 1.1];  % Initial guess for costate (4x1)
options_fsolve = optimoptions('fsolve','Display','iter','TolFun',1e-8,'TolX',1e-8);
for i = 1:length(rhos)
    rho = rhos(i);
    shootingFunc = @(lambda0) shootFunc(lambda0, x0, tf, mu, umax, rho, x_target);
    [lambda0_sol, fval, exitflag, output] = fsolve(shootingFunc, lambda0_guess, options_fsolve);
    fprintf('For \\rho = %g, computed lambda0 = [%g; %g; %g; %g]\n', rho, lambda0_sol(1), lambda0_sol(2), lambda0_sol(3), lambda0_sol(4));
    lambda0_guess = lambda0_sol;  % Update guess for next iteration
end
rho_final = rhos(end);

% Integrate the Full ODE System (State and Costate)
% Full state: y = [r (2); v (2); lambda_r (2); lambda_v (2)]
options_ode = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_sol, y_sol] = ode45(@(t,y) odefun(t,y,mu,umax,rho_final), [0 tf], [x0; lambda0_sol], options_ode);
r_sol = y_sol(:,1:2);       % Spacecraft position
v_sol = y_sol(:,3:4);       % Spacecraft velocity
lambda_r_sol = y_sol(:,5:6);% Costate for r
lambda_v_sol = y_sol(:,7:8);% Costate for v

% Compute Control History
u_sol = zeros(length(t_sol),2);
for k = 1:length(t_sol)
    lam_v = lambda_v_sol(k,:)';
    norm_lam_v = norm(lam_v);
    if norm_lam_v > 1e-8
        Gamma = umax/2*(1 + tanh((norm_lam_v - 1)/rho_final));
        u_sol(k,:) = - Gamma * (lam_v/norm_lam_v);
    else
        u_sol(k,:) = [0; 0];
    end
end

% Compute Hamiltonian Variation Over Time
H_sol = zeros(length(t_sol),1);
for k = 1:length(t_sol)
    r = r_sol(k,:)';
    v = v_sol(k,:)';
    lam_r = lambda_r_sol(k,:)';
    lam_v = lambda_v_sol(k,:)';
    norm_r = norm(r);
    u = u_sol(k,:)';
    H_sol(k) = lam_r'*v + lam_v'*(-r/(norm_r^3)) + lam_v'*u + norm(u);
end
H0 = H_sol(1);
H_diff = H_sol - H0;

% Compute Polar Coordinates for State (for plotting states vs. time)
xSc = r_sol(:,1);
ySc = r_sol(:,2);
rPolar = sqrt(xSc.^2 + ySc.^2);
theta = atan2(ySc, xSc);
vx = v_sol(:,1);
vy = v_sol(:,2);
vr = (xSc.*vx + ySc.*vy)./rPolar;            % Radial velocity
vtheta = (xSc.*vy - ySc.*vx)./rPolar;          % Tangential velocity

% FIGURE 1: Trajectories in Cartesian (Optimal Transfer)
figure('Color','white','Position',[0 0 1000 1000]); hold on;
r0 = 1;              % Earth orbit radius
rMars = a_M;         % Mars orbit radius
t_circle = linspace(0,2*pi,200);
% Plot Earth orbit
plot(r0*cos(t_circle), r0*sin(t_circle), 'k--', 'DisplayName','Earth Orbit');
% Plot Mars orbit
plot(rMars*cos(t_circle), rMars*sin(t_circle), 'r--', 'DisplayName','Mars Orbit');
% Mark Earth start (blue marker)
plot(r0, 0, 'ko','MarkerFaceColor','b','MarkerSize',20, 'DisplayName','Earth Start');
% Mark Mars start (red marker at π)
plot(rMars*cos(pi), rMars*sin(pi), 'ro','MarkerFaceColor','r','MarkerSize',20, 'DisplayName','Mars Start');
% Plot spacecraft trajectory
plot(xSc, ySc, 'b-', 'LineWidth', 3, 'DisplayName', 'Spacecraft Trajectory');
% Mark intercept (final point)
scatter(xSc(end), ySc(end), 1200, 'bp', 'MarkerFaceColor','y','MarkerEdgeColor','k', 'DisplayName', 'Intercept');
axis equal; grid on;
xlabel('x (AU)'); ylabel('y (AU)');
title('Minimum-Fuel Trajectory from Earth to Mars');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

% FIGURE 2: States vs. Time (3 Subplots)
figure('Color','white','Position',[1000 0 1000 1000]);
% Subplot 1: Radial distance vs. time
subplot(3,1,1); hold on; grid on;
title('States vs. Time');
ylabel('r (AU)');
plot(t_sol, rPolar, 'b-', 'LineWidth', 3, 'DisplayName', 'r');
legend('Location','Best');
ax = gca; ax.FontSize = 20;
% Subplot 2: Radial and Tangential Velocities vs. time
subplot(3,1,2); hold on; grid on;
ylabel('Velocity (AU/yr)');
plot(t_sol, vr, 'b--', 'LineWidth', 3, 'DisplayName', 'v_r');
plot(t_sol, vtheta, 'b-', 'LineWidth', 3, 'DisplayName', 'v_\theta');
legend('Location','Best');
ax = gca; ax.FontSize = 20;
% Subplot 3: Angular Position vs. time
subplot(3,1,3); hold on; grid on;
xlabel('Time (yr)'); ylabel('\theta (rad)');
plot(t_sol, theta, 'b-', 'LineWidth', 3, 'DisplayName', '\theta');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

% FIGURE 3: Hamiltonian Difference vs. Time
figure('Color','white','Position',[0 1000 1000 1000]); hold on; grid on;
title('H(t) - H(0)');
xlabel('Time (yr)'); ylabel('H(t) - H(0)');
plot(t_sol, H_diff, 'b-', 'LineWidth', 3, 'DisplayName', 'H(t)-H(0)');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

% FIGURE 4: Switching Function and Thrust Profile vs. Time
% Compute the switching function S = ||lambda_v|| - 1 and thrust magnitude Gamma.
S_sol = zeros(length(t_sol),1);
Gamma_sol = zeros(length(t_sol),1);
for k = 1:length(t_sol)
    norm_lv = norm(lambda_v_sol(k,:)');
    S_sol(k) = norm_lv - 1;
    Gamma_sol(k) = umax/2*(1 + tanh((norm_lv - 1)/rho_final));
end
figure('Color','white','Position',[1000 1000 1000 1000]); hold on; grid on;
plot(t_sol, S_sol, 'b-', 'LineWidth', 3, 'DisplayName', 'Switching Function S = ||\lambda_v|| - 1');
plot(t_sol, Gamma_sol, 'r--', 'LineWidth', 3, 'DisplayName', 'Thrust Magnitude \Gamma*');
xlabel('Time (yr)');
ylabel('Value');
title('Switching Function and Thrust Profile vs. Time');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

% FIGURE 5: Control Time History from Final Iteration
figure('Color','white','Position',[1000 1000 1000 1000]); hold on; grid on;
plot(t_sol, u_sol(:,1), 'b-', 'LineWidth', 3, 'DisplayName', 'u_x');
plot(t_sol, u_sol(:,2), 'r-', 'LineWidth', 3, 'DisplayName', 'u_y');
plot(t_sol, vecnorm([u_sol(:,1) u_sol(:,2)],2,2), 'k-', 'LineWidth', 3, 'DisplayName', 'u_{mag}');
xlabel('Time (yr)');
ylabel('Control Acceleration');
title('Control Time History');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

% FIGURE 6: Costate Time History from Final Iteration
figure('Color','white','Position',[2000 1000 1000 1000]);
subplot(2,1,1); hold on; grid on;
plot(t_sol, lambda_r_sol(:,1), 'b-', 'LineWidth', 3, 'DisplayName', '\lambda_{r}');
plot(t_sol, lambda_r_sol(:,2), 'r-', 'LineWidth', 3, 'DisplayName', '\lambda_{\theta}');
xlabel('Time (yr)');
ylabel('Position');
title('Costate Time History: \lambda_r');
legend('Location','Best');
ax = gca; ax.FontSize = 20;
subplot(2,1,2); hold on; grid on;
plot(t_sol, lambda_v_sol(:,1), 'b-', 'LineWidth', 3, 'DisplayName', '\lambda_{vr}');
plot(t_sol, lambda_v_sol(:,2), 'r-', 'LineWidth', 3, 'DisplayName', '\lambda_{v\theta}');
xlabel('Time (yr)');
ylabel('Velocity');
title('Costate Time History: \lambda_v');
legend('Location','Best');
ax = gca; ax.FontSize = 20;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested Function: ODE Dynamics for State and Costate
function dydt = odefun(t,y,mu,umax,rho)
    % y = [r (2); v (2); lambda_r (2); lambda_v (2)]
    r = y(1:2);
    v = y(3:4);
    lambda_r = y(5:6);
    lambda_v = y(7:8);
    norm_r = norm(r);
    norm_lambda_v = norm(lambda_v);
    if norm_lambda_v > 1e-8
        Gamma = umax/2*(1 + tanh((norm_lambda_v - 1)/rho));
        u = - Gamma * (lambda_v/norm_lambda_v);
    else
        u = [0;0];
    end

    drdt = v;
    dvdt = - r/(norm_r^3) + u;
    dlambda_v_dt = - lambda_r;
    A = (1/(norm_r^3))*eye(2) - 3*(r*r')/(norm_r^5);
    dlambda_r_dt = A * lambda_v;
    dydt = [drdt; dvdt; dlambda_r_dt; dlambda_v_dt];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested Function: Shooting Function
function F = shootFunc(lambda0, x0, tf, mu, umax, rho, x_target)
    y0 = [x0; lambda0];
    options_ode_inner = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [~, y] = ode45(@(t,y) odefun(t,y,mu,umax,rho), [0 tf], y0, options_ode_inner);
    x_tf = y(end,1:4)';
    F = x_tf - x_target;
end

%% Monte Carlo simulation (Problem 10b)
close all;

M    = 20;           % number of runs
dt   = 1e-3;         % time step
t_MC = 0:dt:tf;      % time vector
N    = numel(t_MC)-1;

% dispersion & disturbance
sigma_pos  = 0.3e-4;
sigma_vel = sigma_pos;
P0         = diag([sigma_pos^2, sigma_pos^2, sigma_vel^2, sigma_vel^2]);
sigma_dist = 1e-3; % disturbance accel
G          = [zeros(2); sigma_dist*eye(2)];

% reference interpolation
Y_ref = interp1(t_sol, [r_sol, v_sol], t_MC);    % (N+1)x4
U_ref = interp1(t_sol,     u_sol,       t_MC);    % (N+1)x2

%-------------------------------------------------------------
% Simulate an ensemble of M stochastic trajectories (Monte Carlo)
%-------------------------------------------------------------

% Preallocate a 4×(N+1)×M array:
%   4 states per time step (x, y, vx, vy)
%   N+1 time steps from k=0 to k=N
%   M independent Monte Carlo runs
Xmc = zeros(4, N+1, M);

% Loop over each Monte Carlo sample
for i = 1:M
    %---------------------------------------------------------
    % 1) Sample the initial state with Gaussian uncertainty
    %    x0:     nominal initial state [x; y; vx; vy]
    %    P0:     covariance of initial‐condition dispersion
    %    chol(P0,'lower')*randn(4,1):
    %            draws a 4×1 zero‐mean Gaussian with cov P0
    %---------------------------------------------------------
    Xmc(:,1,i) = x0 + chol(P0,'lower')*randn(4,1);
    
    %---------------------------------------------------------
    % 2) Propagate forward in time using Euler–Maruyama
    %    dt:     time step size
    %    U_ref:  precomputed nominal control history [2×(N+1)]
    %    G:      noise‐injection matrix for disturbances
    %---------------------------------------------------------
    for k = 1:N
        % a) Extract the current state for run i at time step k
        %    xk = [ x; y; vx; vy ]
        xk = Xmc(:,k,i);
        
        % b) Compute the deterministic drift f(xk)
        %    dr/dt = v
        %    dv/dt = gravity + control
        %      gravity term    = -r / ||r||^3
        %      control term    = U_ref(k,:)'  (2×1)
        f = [ 
            xk(3:4);                                  % [vx; vy]
            - xk(1:2) / norm(xk(1:2))^3 + U_ref(k,:)'  % acceleration
        ];
        
        % c) Sample a Brownian increment in the acceleration channel
        %    dW = sqrt(dt) * N(0,I_2)
        dW = sqrt(dt) * randn(2,1);
        
        % d) Euler–Maruyama update:
        %    Xmc(:,k+1,i) = xk + f*dt + G*dW
        %    where G maps the 2×1 disturbance into the 4×1 state update
        Xmc(:,k+1,i) = xk + f*dt + G*dW;
    end
end

% compute polar states for MC runs
Rmc  = zeros(M, N+1);
Thmc = zeros(M, N+1);
for i = 1:M
    xi = squeeze(Xmc(1,:,i));
    yi = squeeze(Xmc(2,:,i));
    Rmc(i,:)  = sqrt(xi.^2 + yi.^2);
    Thmc(i,:) = atan2(yi, xi);
end

% compute polar velocities for MC runs
VRmc  = zeros(M, N+1);
VThmc = zeros(M, N+1);
for i = 1:M
    xi  = squeeze(Xmc(1,:,i));
    yi  = squeeze(Xmc(2,:,i));
    vxi = squeeze(Xmc(3,:,i));
    vyi = squeeze(Xmc(4,:,i));
    R   = Rmc(i,:);
    VRmc(i,:)  = (xi.*vxi + yi.*vyi) ./ R;
    VThmc(i,:) = (xi.*vyi - yi.*vxi) ./ R;
end

% compute polar states & velocities for reference
Rref   = sqrt(Y_ref(:,1).^2 + Y_ref(:,2).^2);
Thref  = atan2(Y_ref(:,2), Y_ref(:,1));
VRref  = (Y_ref(:,1).*Y_ref(:,3) + Y_ref(:,2).*Y_ref(:,4)) ./ Rref;
VThref = (Y_ref(:,1).*Y_ref(:,4) - Y_ref(:,2).*Y_ref(:,3)) ./ Rref;

%—— Figure: polar-coordinate states vs. time ——
figure('Color','white','Position',[100 100 900 1000]);

% 1) radial distance
subplot(4,1,1); hold on; grid on;
for i = 1:M
    plot(t_MC, Rmc(i,:), 'LineWidth',1, 'HandleVisibility','off');
end
plot(t_MC, Rref, 'k-', 'LineWidth',3);
ylabel('r (AU)','FontSize',20);
title('Radial Distance vs. Time','FontSize',20);
ax = gca; ax.FontSize = 20;

% 2) angular position
subplot(4,1,2); hold on; grid on;
for i = 1:M
    plot(t_MC, Thmc(i,:), 'LineWidth',1, 'HandleVisibility','off');
end
plot(t_MC, Thref, 'k-', 'LineWidth',3);
ylabel('\theta (rad)','FontSize',20);
title('Angular Position vs. Time','FontSize',20);
ax = gca; ax.FontSize = 20;

% 3) radial velocity
subplot(4,1,3); hold on; grid on;
for i = 1:M
    plot(t_MC, VRmc(i,:), 'LineWidth',1, 'HandleVisibility','off');
end
plot(t_MC, VRref, 'k-', 'LineWidth',3);
ylabel('v_r (AU/yr)','FontSize',20);
title('Radial Velocity vs. Time','FontSize',20);
ax = gca; ax.FontSize = 20;

% 4) tangential velocity
subplot(4,1,4); hold on; grid on;
for i = 1:M
    plot(t_MC, VThmc(i,:), 'LineWidth',1, 'HandleVisibility','off');
end
plot(t_MC, VThref, 'k-', 'LineWidth',3);
xlabel('Time (yr)','FontSize',20);
ylabel('v_\theta (AU/yr)','FontSize',20);
title('Tangential Velocity vs. Time','FontSize',20);
ax = gca; ax.FontSize = 20;


%—— Figure: Monte Carlo Transfer Trajectories in Cartesian Space ——
theta_circle = linspace(0,2*pi,200);

figure('Color','white','Position',[100 100 1000 1000]); hold on; grid on;

% Earth and Mars orbits
plot(cos(theta_circle),          sin(theta_circle),          'k--', 'LineWidth',1.5, 'DisplayName','Earth Orbit');
plot(a_M*cos(theta_circle), a_M*sin(theta_circle), 'r--', 'LineWidth',1.5, 'DisplayName','Mars Orbit');

% Monte Carlo trajectories (hidden from legend)
for i = 1:M
    xi = squeeze(Xmc(1,:,i));
    yi = squeeze(Xmc(2,:,i));
    plot(xi, yi, 'LineWidth',1, 'HandleVisibility','off');
end

% Nominal trajectory
plot(r_sol(:,1), r_sol(:,2), 'k-', 'LineWidth',3, 'DisplayName','Nominal Trajectory');

% Start & intercept markers
plot(1, 0,                           'ko', 'MarkerFaceColor','b', 'MarkerSize',12, 'DisplayName','Earth Start');
plot(a_M*cos(pi), a_M*sin(pi),      'ro', 'MarkerFaceColor','r', 'MarkerSize',12, 'DisplayName','Mars Start');
scatter(r_sol(end,1), r_sol(end,2), 250,  'p', 'MarkerFaceColor','y', 'MarkerEdgeColor','k', 'DisplayName','Intercept');

% Labels, title, legend, styling
xlabel('x (AU)','FontSize',20);
ylabel('y (AU)','FontSize',20);
title('Monte Carlo Transfer Trajectories vs. Nominal Optimal Path','FontSize',20);
legend('Location','Best','FontSize',16);

axis equal;
ax = gca; 
ax.FontSize = 20;

%—— Figure: Cartesian-coordinate states vs. time ——
figure('Color','white','Position',[700 100 900 1000]);

% 1) x position
subplot(4,1,1); hold on; grid on;
for i = 1:M
    plot(t_MC, squeeze(Xmc(1,:,i)), 'LineWidth',1, 'HandleVisibility','off');
end
plot(t_MC, Y_ref(:,1), 'k-', 'LineWidth',3);
ylabel('x (AU)','FontSize',20);
title('Cartesian State: x vs. Time','FontSize',20);
ax = gca; ax.FontSize = 20;

% 2) y position
subplot(4,1,2); hold on; grid on;
for i = 1:M
    plot(t_MC, squeeze(Xmc(2,:,i)), 'LineWidth',1, 'HandleVisibility','off');
end
plot(t_MC, Y_ref(:,2), 'k-', 'LineWidth',3);
ylabel('y (AU)','FontSize',20);
title('Cartesian State: y vs. Time','FontSize',20);
ax = gca; ax.FontSize = 20;

% 3) x velocity
subplot(4,1,3); hold on; grid on;
for i = 1:M
    plot(t_MC, squeeze(Xmc(3,:,i)), 'LineWidth',1, 'HandleVisibility','off');
end
plot(t_MC, Y_ref(:,3), 'k-', 'LineWidth',3);
ylabel('v_x (AU/yr)','FontSize',20);
title('Cartesian State: v_x vs. Time','FontSize',20);
ax = gca; ax.FontSize = 20;

% 4) y velocity
subplot(4,1,4); hold on; grid on;
for i = 1:M
    plot(t_MC, squeeze(Xmc(4,:,i)), 'LineWidth',1, 'HandleVisibility','off');
end
plot(t_MC, Y_ref(:,4), 'k-', 'LineWidth',3);
xlabel('Time (yr)','FontSize',20);
ylabel('v_y (AU/yr)','FontSize',20);
title('Cartesian State: v_y vs. Time','FontSize',20);
ax = gca; ax.FontSize = 20;

%
% HW10_cov_analysis.m
% Part (c) & (d): Linear covariance analysis vs. Monte Carlo for minimum-fuel transfer
% AAE590 HW10, Derek Yu

%-----------------------------------------------------
%--- Load or compute nominal trajectory into:
%    t_sol, r_sol, v_sol, u_sol, tf, mu
%  If running standalone, you can:
%    load('HW10_solution.mat','t_sol','r_sol','v_sol','u_sol','tf','mu');
%-----------------------------------------------------

%--- Common parameters ---
t_cov  = 0:dt_cov:tf;     % discrete times
dt_MC = dt;    % use the same dt you already defined for t_MC
N_MC  = N;     % and the same N
% — DOWN-SAMPLE MONTE CARLO FOR COVARIANCE COMPARISON —
dt_cov  = 0.1;
idx_cov = 1:round(dt_cov/dt_MC):N_MC+1;   % indices into t_MC/Y_ref/Xmc
time    = t_MC(idx_cov);
N_cov   = numel(time)-1;

% — PROPAGATE LINEAR COVARIANCE (re-using true G) —
Pcell = cell(1, N_cov+1);
Pcell{1} = P0;
for k = 1:N_cov
    r_k        = Y_ref(idx_cov(k),1:2)';       % nominal position
    [A_c, ~]   = linAB(r_k);                   % only the A matrix
    [Phi, Q]   = discretizeCov(A_c, G, dt_cov);% use G from MC (with σ_dist)
    Pcell{k+1} = Phi * Pcell{k} * Phi' + Q;
end

% — PLOT DEVIATIONS vs. ±3σ —  
figure('Color','white','Position',[1500 100 1200 800]);
states = {'x','y','v_x','v_y'};
for j = 1:4
    subplot(2,2,j); hold on; grid on;
      data_j = squeeze(Xmc(j, idx_cov, :));   % (N_cov+1)xM
      nom_j  = Y_ref(idx_cov, j);             % (N_cov+1)x1
      dev    = data_j - nom_j;                % true Cartesian deviations
      plot(time, dev); % all runs
      sigma_j = 3 * cellfun(@(Pmat) sqrt(Pmat(j,j)), Pcell);
      plot(time, +sigma_j, 'k:','LineWidth',4);
      plot(time, -sigma_j, 'k:','LineWidth',4);
    title(sprintf('Deviation of %s vs. \\pm3\\sigma', states{j}), 'FontSize',20);
    xlabel('Time','FontSize',20); ylabel('Deviation','FontSize',20);
    set(gca,'FontSize',20);
end


%-----------------------------------------------------
% Local function: continuous-time linearization
%-----------------------------------------------------
function [A_c,B_c] = linAB(r)
    % x = [r;v], f = [v; -r/|r|^3 + u]
    R    = norm(r);
    I2   = eye(2);
    dfdx = -(1/R^3)*I2 + 3*(r*r')/R^5;
    A_c  = [ zeros(2), I2;
            dfdx,    zeros(2) ];
    B_c  = [ zeros(2); I2 ];
end

%-----------------------------------------------------
% Local function: Van-Loan discretization of covariance
%-----------------------------------------------------
function [Phi,Qd] = discretizeCov(A_c,G_c,dt)
    n    = size(A_c,1);
    Mbig = [ -A_c,      G_c*G_c';
             zeros(n),   A_c'     ];
    E    = expm(Mbig*dt);
    Phi  = expm(A_c*dt);
    Qd   = Phi * E(1:n,n+1:end);
end