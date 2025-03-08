% =========================================================================
% 
% Filename:       HW4.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE590 - Applied Control in Astronautics
% Professor:      Dr. Kenshiro Oguri
% Assignment:     HW 4
% Semester:       Spring 2025
% 
% Description: Homework 4
%
% =========================================================================

%%
clear; clc; close all;
format long

% PARAMETERS
mu = 1;              % Non-dimensional gravitational parameter
umax = 0.1;          % Maximum control magnitude
tf = 6.5;            % Fixed final time (non-dimensional)

% Smoothing parameters for continuation
rhos = [1, 0.1, 1e-2, 1e-3];

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