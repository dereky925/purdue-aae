% =========================================================================
% 
% Filename:       Compare_HW4_Enhanced_ComparisonPlots_Polar_Control.m
% Author:         Your Name
% Institution:    Purdue University
% Course:         AAE590 - Applied Control in Astronautics
% Professor:      Dr. Kenshiro Oguri
% Assignment:     HW 4 Comparison (with Polar & Control Plots)
% Semester:       Spring 2025
% 
% Description: Compares the solutions of the constant-mass (4-state) and 
%              variable-mass (5-state) optimal transfer problems. In addition
%              to overlayed and difference plots in Cartesian coordinates,
%              the script also plots the state and costate histories in 
%              polar coordinates, and compares the control input time histories.
%
% =========================================================================

% Clear workspace and set formatting
clear; clc; close all;
format long

% PARAMETERS (common)
mu    = 1;              % Non-dimensional gravitational parameter
umax  = 0.1;            % Maximum control magnitude
tf    = 8;              % Fixed final time (non-dimensional)
rhos  = [1, 0.1, 1e-2, 1e-3];
rho_final = rhos(end);

% Mars orbit parameters (common)
a_M   = 1.524;          % Mars semi-major axis (AU)
nu_M0 = pi;             % Mars initial true anomaly (radians)
omega_M = sqrt(mu/(a_M^3));  % Mars angular rate
marsState = @(t) [ a_M*cos(nu_M0 + omega_M*t);
                   a_M*sin(nu_M0 + omega_M*t);
                  -a_M*omega_M*sin(nu_M0 + omega_M*t);
                   a_M*omega_M*cos(nu_M0 + omega_M*t) ];
x_target = marsState(tf);  % Terminal condition: [r; v]

% Specific impulse product for variable-mass problem
Ispg0 = 4000*9.8;

% Initial Conditions
% Constant mass (4-state): [rx; ry; vx; vy]
x0_const = [1; 0; 0; 1];
% Variable mass (5-state): [rx; ry; vx; vy; m]
x0_var   = [1; 0; 0; 1; 1];

% SOLVE CONSTANT MASS PROBLEM (4-state)
lambda0_guess_const = [0.8; 0.1; 0.2; 1.1];  % initial guess for costate
options_fsolve = optimoptions('fsolve','Display','iter','TolFun',1e-8,'TolX',1e-8);
for i = 1:length(rhos)
    rho = rhos(i);
    shootingFunc_const = @(lambda0) constantMass_shootFunc(lambda0, x0_const, tf, mu, umax, rho, x_target);
    [lambda0_sol_const, ~, ~, ~] = fsolve(shootingFunc_const, lambda0_guess_const, options_fsolve);
    fprintf('Constant Mass, for rho = %g, computed lambda0 = [%g; %g; %g; %g]\n', ...
        rho, lambda0_sol_const(1), lambda0_sol_const(2), lambda0_sol_const(3), lambda0_sol_const(4));
    lambda0_guess_const = lambda0_sol_const;  % update guess
end

% SOLVE VARIABLE MASS PROBLEM (5-state)
lambda0_guess_var = [0.849671; 0.101611; 0.201818; 1.08465; 1];  % user-specified guess
for i = 1:length(rhos)
    rho = rhos(i);
    shootingFunc_var = @(lambda0) variableMass_shootFunc(lambda0, x0_var, tf, mu, umax, rho, Ispg0, x_target);
    [lambda0_sol_var, ~, ~, ~] = fsolve(shootingFunc_var, lambda0_guess_var, options_fsolve);
    fprintf('Variable Mass, for rho = %g, computed lambda0 = [%g; %g; %g; %g; %g]\n', ...
        rho, lambda0_sol_var(1), lambda0_sol_var(2), lambda0_sol_var(3), lambda0_sol_var(4), lambda0_sol_var(5));
    lambda0_guess_var = lambda0_sol_var;  % update guess
end

% INTEGRATE ODEs

options_ode = odeset('RelTol',1e-12,'AbsTol',1e-12);

% Constant mass integration
[t_const, y_const] = ode45(@(t,y) constantMass_odefun(t,y,mu,umax,rho_final), [0 tf], [x0_const; lambda0_sol_const], options_ode);
r_const = y_const(:,1:2);
v_const = y_const(:,3:4);
lambda_r_const = y_const(:,5:6);
lambda_v_const = y_const(:,7:8);

% Variable mass integration
[t_var, y_var] = ode45(@(t,y) variableMass_odefun(t,y,mu,umax,rho_final,Ispg0), [0 tf], [x0_var; lambda0_sol_var], options_ode);
r_var = y_var(:,1:2);
v_var = y_var(:,3:4);
m_var = y_var(:,5);
lambda_r_var = y_var(:,6:7);
lambda_v_var = y_var(:,8:9);
lambda_m_var = y_var(:,10);

% COMPUTE CONTROL HISTORY

% Constant mass control (switching function: ||lambda_v|| - 1)
u_const = zeros(length(t_const),2);
for k = 1:length(t_const)
    lam_v = lambda_v_const(k,:)';
    norm_lv = norm(lam_v);
    if norm_lv > 1e-8
        Gamma = umax/2*(1 + tanh((norm_lv - 1)/rho_final));
        u_const(k,:) = - Gamma * (lam_v/norm_lv);
    else
        u_const(k,:) = [0 0];
    end
end

% Variable mass control (switching function: (||lambda_v||/m)+ (lambda_m/Ispg0) )
u_var = zeros(length(t_var),2);
for k = 1:length(t_var)
    lam_v = lambda_v_var(k,:)';
    norm_lv = norm(lam_v);
    if norm_lv > 1e-8
        S_val = (norm_lv/m_var(k)) + (lambda_m_var(k)/Ispg0);
        Gamma = umax/2*(1 + tanh((S_val - 1)/rho_final));
        u_var(k,:) = - Gamma * (lam_v/norm_lv);
    else
        u_var(k,:) = [0 0];
    end
end

% PLOT COMPARISONS (Cartesian)

% ----------------------------
% 1. Trajectory Comparison (with markers)
% ----------------------------
figure('Color','white','Position',[0 0 1000 1000]); hold on; grid on;
r0 = 1;          % Earth orbit radius
rMars = a_M;     % Mars orbit radius
t_circle = linspace(0,2*pi,200);
% Plot Earth and Mars orbits
plot(r0*cos(t_circle), r0*sin(t_circle), 'k--', 'LineWidth', 2, 'DisplayName','Earth Orbit');
plot(rMars*cos(t_circle), rMars*sin(t_circle), 'r--', 'LineWidth', 2, 'DisplayName','Mars Orbit');
% Mark starting points: Earth and Mars
plot(r0, 0, 'ko','MarkerFaceColor','b','MarkerSize',20, 'DisplayName','Earth Start');
plot(rMars*cos(pi), rMars*sin(pi), 'ro','MarkerFaceColor','r','MarkerSize',20, 'DisplayName','Mars Start');
% Plot trajectories
plot(r_const(:,1), r_const(:,2), 'b-', 'LineWidth', 3, 'DisplayName','Trajectory (Const Mass)');
plot(r_var(:,1), r_var(:,2), 'm--', 'LineWidth', 3, 'DisplayName','Trajectory (Var Mass)');
% Mark intercept points
scatter(r_const(end,1), r_const(end,2), 1200, 'bp','MarkerFaceColor','y','MarkerEdgeColor','k','DisplayName','Intercept (Const)');
scatter(r_var(end,1), r_var(end,2), 1200, 'mp','MarkerFaceColor','c','MarkerEdgeColor','k','DisplayName','Intercept (Var)');
xlabel('x (AU)'); ylabel('y (AU)');
title('Trajectory Comparison');
legend('Location','Best');
axis equal;

% ----------------------------
% 2. Trajectory Differences (Cartesian)
% ----------------------------
% Interpolate variable-mass solution onto constant-mass time grid.
t_common = t_const;
r_var_interp = interp1(t_var, r_var, t_common);
diff_r = r_const - r_var_interp;
diff_norm = vecnorm(diff_r,2,2);

figure('Color','white','Position',[1000 0 1000 1000]); hold on; grid on;
subplot(3,1,1);
plot(t_common, diff_r(:,1), 'k-', 'LineWidth', 2);
ylabel('\Delta x (AU)'); title('Trajectory Differences: Constant Mass - Variable Mass');
subplot(3,1,2);
plot(t_common, diff_r(:,2), 'k-', 'LineWidth', 2);
ylabel('\Delta y (AU)');
subplot(3,1,3);
plot(t_common, diff_norm, 'k-', 'LineWidth', 2);
xlabel('Time (yr)'); ylabel('||\Delta r|| (AU)');
legend('Difference'); grid on;

% ----------------------------
% 3. States Comparison (Cartesian)
% ----------------------------
% Interpolate variable mass states onto constant mass time grid.
y_var_interp = interp1(t_var, y_var, t_common);
r_const_interp = r_const;
v_const_interp = v_const;
r_var_interp = y_var_interp(:,1:2);
v_var_interp = y_var_interp(:,3:4);

figure('Color','white','Position',[0 1000 1000 1000]);
subplot(2,2,1); hold on; grid on;
plot(t_const, r_const_interp(:,1), 'b-', 'LineWidth', 2, 'DisplayName','r_x (Const)');
plot(t_common, r_var_interp(:,1), 'm--', 'LineWidth', 2, 'DisplayName','r_x (Var)');
xlabel('Time (yr)'); ylabel('r_x (AU)'); title('State Comparison: r_x');
legend('Location','Best');
subplot(2,2,2); hold on; grid on;
plot(t_const, r_const_interp(:,2), 'b-', 'LineWidth', 2, 'DisplayName','r_y (Const)');
plot(t_common, r_var_interp(:,2), 'm--', 'LineWidth', 2, 'DisplayName','r_y (Var)');
xlabel('Time (yr)'); ylabel('r_y (AU)'); title('State Comparison: r_y');
legend('Location','Best');
subplot(2,2,3); hold on; grid on;
plot(t_const, v_const_interp(:,1), 'b-', 'LineWidth', 2, 'DisplayName','v_x (Const)');
plot(t_common, v_var_interp(:,1), 'm--', 'LineWidth', 2, 'DisplayName','v_x (Var)');
xlabel('Time (yr)'); ylabel('v_x (AU/yr)'); title('State Comparison: v_x');
legend('Location','Best');
subplot(2,2,4); hold on; grid on;
plot(t_const, v_const_interp(:,2), 'b-', 'LineWidth', 2, 'DisplayName','v_y (Const)');
plot(t_common, v_var_interp(:,2), 'm--', 'LineWidth', 2, 'DisplayName','v_y (Var)');
xlabel('Time (yr)'); ylabel('v_y (AU/yr)'); title('State Comparison: v_y');
legend('Location','Best');

% ----------------------------
% 4. Costates Comparison (Cartesian)
% ----------------------------
% For constant mass, costates are directly available.
lambda_r_var_interp = interp1(t_var, lambda_r_var, t_common);
lambda_v_var_interp = interp1(t_var, lambda_v_var, t_common);

figure('Color','white','Position',[1000 1000 1000 1000]);
subplot(2,2,1); hold on; grid on;
plot(t_const, lambda_r_const(:,1), 'b-', 'LineWidth', 2, 'DisplayName','\lambda_{r_x} (Const)');
plot(t_common, lambda_r_var_interp(:,1), 'm--', 'LineWidth', 2, 'DisplayName','\lambda_{r_x} (Var)');
xlabel('Time (yr)'); ylabel('\lambda_{r_x}'); title('Costate Comparison: \lambda_{r_x}');
legend('Location','Best');
subplot(2,2,2); hold on; grid on;
plot(t_const, lambda_r_const(:,2), 'b-', 'LineWidth', 2, 'DisplayName','\lambda_{r_y} (Const)');
plot(t_common, lambda_r_var_interp(:,2), 'm--', 'LineWidth', 2, 'DisplayName','\lambda_{r_y} (Var)');
xlabel('Time (yr)'); ylabel('\lambda_{r_y}'); title('Costate Comparison: \lambda_{r_y}');
legend('Location','Best');
subplot(2,2,3); hold on; grid on;
plot(t_const, lambda_v_const(:,1), 'b-', 'LineWidth', 2, 'DisplayName','\lambda_{v_x} (Const)');
plot(t_common, lambda_v_var_interp(:,1), 'm--', 'LineWidth', 2, 'DisplayName','\lambda_{v_x} (Var)');
xlabel('Time (yr)'); ylabel('\lambda_{v_x}'); title('Costate Comparison: \lambda_{v_x}');
legend('Location','Best');
subplot(2,2,4); hold on; grid on;
plot(t_const, lambda_v_const(:,2), 'b-', 'LineWidth', 2, 'DisplayName','\lambda_{v_y} (Const)');
plot(t_common, lambda_v_var_interp(:,2), 'm--', 'LineWidth', 2, 'DisplayName','\lambda_{v_y} (Var)');
xlabel('Time (yr)'); ylabel('\lambda_{v_y}'); title('Costate Comparison: \lambda_{v_y}');
legend('Location','Best');

% POLAR COORDINATE PLOTS FOR STATES

% Interpolate variable mass states onto t_common:
y_var_interp = interp1(t_var, y_var, t_common);
% Constant Mass:
x_const = r_const(:,1); 
y_const = r_const(:,2);
rPolar_const = sqrt(x_const.^2 + y_const.^2);
theta_const = atan2(y_const, x_const);
vx_const = v_const(:,1); 
vy_const = v_const(:,2);
vr_const = (x_const .* vx_const + y_const .* vy_const) ./ rPolar_const;
vtheta_const = (x_const .* vy_const - y_const .* vx_const) ./ rPolar_const;
% Variable Mass:
r_var_interp = y_var_interp(:,1:2);
v_var_interp = y_var_interp(:,3:4);
x_var = r_var_interp(:,1); 
y_var = r_var_interp(:,2);
rPolar_var = sqrt(x_var.^2 + y_var.^2);
theta_var = atan2(y_var, x_var);
vx_var = v_var_interp(:,1); 
vy_var = v_var_interp(:,2);
vr_var = (x_var .* vx_var + y_var .* vy_var) ./ rPolar_var;
vtheta_var = (x_var .* vy_var - y_var .* vx_var) ./ rPolar_var;

figure('Color','white','Position',[0 0 1000 1000]);
subplot(2,2,1);
plot(t_common, rPolar_const, 'b-', 'LineWidth',2, 'DisplayName', 'r (Const)');
hold on;
plot(t_common, rPolar_var, 'm--', 'LineWidth',2, 'DisplayName', 'r (Var)');
xlabel('Time (yr)'); ylabel('Radial Distance (AU)');
title('r vs. Time');
legend('Location','Best');
subplot(2,2,2);
plot(t_common, theta_const, 'b-', 'LineWidth',2, 'DisplayName', '\theta (Const)');
hold on;
plot(t_common, theta_var, 'm--', 'LineWidth',2, 'DisplayName', '\theta (Var)');
xlabel('Time (yr)'); ylabel('Angle (rad)');
title('\theta vs. Time');
legend('Location','Best');
subplot(2,2,3);
plot(t_common, vr_const, 'b-', 'LineWidth',2, 'DisplayName', 'v_r (Const)');
hold on;
plot(t_common, vr_var, 'm--', 'LineWidth',2, 'DisplayName', 'v_r (Var)');
xlabel('Time (yr)'); ylabel('Radial Velocity (AU/yr)');
title('v_r vs. Time');
legend('Location','Best');
subplot(2,2,4);
plot(t_common, vtheta_const, 'b-', 'LineWidth',2, 'DisplayName', 'v_\theta (Const)');
hold on;
plot(t_common, vtheta_var, 'm--', 'LineWidth',2, 'DisplayName', 'v_\theta (Var)');
xlabel('Time (yr)'); ylabel('Tangential Velocity (AU/yr)');
title('v_\theta vs. Time');
legend('Location','Best');

% POLAR COORDINATE PLOTS FOR COSTATES

% Convert constant mass costates to polar components:
lambda_r_const_polar_r = (x_const .* lambda_r_const(:,1) + y_const .* lambda_r_const(:,2)) ./ rPolar_const;
lambda_r_const_polar_theta = (-y_const .* lambda_r_const(:,1) + x_const .* lambda_r_const(:,2)) ./ rPolar_const;
lambda_v_const_polar_r = (x_const .* lambda_v_const(:,1) + y_const .* lambda_v_const(:,2)) ./ rPolar_const;
lambda_v_const_polar_theta = (-y_const .* lambda_v_const(:,1) + x_const .* lambda_v_const(:,2)) ./ rPolar_const;
% Convert variable mass costates to polar components (interpolated onto t_common):
lambda_r_var_interp = interp1(t_var, lambda_r_var, t_common);
lambda_v_var_interp = interp1(t_var, lambda_v_var, t_common);
lambda_r_var_polar_r = (x_var .* lambda_r_var_interp(:,1) + y_var .* lambda_r_var_interp(:,2)) ./ rPolar_var;
lambda_r_var_polar_theta = (-y_var .* lambda_r_var_interp(:,1) + x_var .* lambda_r_var_interp(:,2)) ./ rPolar_var;
lambda_v_var_polar_r = (x_var .* lambda_v_var_interp(:,1) + y_var .* lambda_v_var_interp(:,2)) ./ rPolar_var;
lambda_v_var_polar_theta = (-y_var .* lambda_v_var_interp(:,1) + x_var .* lambda_v_var_interp(:,2)) ./ rPolar_var;

figure('Color','white','Position',[1000 0 1000 1000]);
subplot(2,2,1);
plot(t_common, lambda_r_const_polar_r, 'b-', 'LineWidth',2, 'DisplayName', '\lambda_{r_r} (Const)');
hold on;
plot(t_common, lambda_r_var_polar_r, 'm--', 'LineWidth',2, 'DisplayName', '\lambda_{r_r} (Var)');
xlabel('Time (yr)'); ylabel('\lambda_{r_r}');
title('Radial Component of \lambda_r');
legend('Location','Best');
subplot(2,2,2);
plot(t_common, lambda_r_const_polar_theta, 'b-', 'LineWidth',2, 'DisplayName', '\lambda_{r_\theta} (Const)');
hold on;
plot(t_common, lambda_r_var_polar_theta, 'm--', 'LineWidth',2, 'DisplayName', '\lambda_{r_\theta} (Var)');
xlabel('Time (yr)'); ylabel('\lambda_{r_\theta}');
title('Angular Component of \lambda_r');
legend('Location','Best');
subplot(2,2,3);
plot(t_common, lambda_v_const_polar_r, 'b-', 'LineWidth',2, 'DisplayName', '\lambda_{v_r} (Const)');
hold on;
plot(t_common, lambda_v_var_polar_r, 'm--', 'LineWidth',2, 'DisplayName', '\lambda_{v_r} (Var)');
xlabel('Time (yr)'); ylabel('\lambda_{v_r}');
title('Radial Component of \lambda_v');
legend('Location','Best');
subplot(2,2,4);
plot(t_common, lambda_v_const_polar_theta, 'b-', 'LineWidth',2, 'DisplayName', '\lambda_{v_\theta} (Const)');
hold on;
plot(t_common, lambda_v_var_polar_theta, 'm--', 'LineWidth',2, 'DisplayName', '\lambda_{v_\theta} (Var)');
xlabel('Time (yr)'); ylabel('\lambda_{v_\theta}');
title('Angular Component of \lambda_v');
legend('Location','Best');

% CONTROL INPUT TIME HISTORY COMPARISON

% Compute control norm for constant mass solution
u_norm_const = vecnorm(u_const, 2, 2);

% Compute control norm for variable mass solution
u_norm_var = vecnorm(u_var, 2, 2);

% Interpolate variable mass control norm onto the constant mass time grid (t_const)
u_norm_var_interp = interp1(t_var, u_norm_var, t_const);

% Plot the comparison of control input norm
figure('Color','white','Position',[2000 0 1000 1000]); hold on; grid on;
plot(t_const, u_norm_const, 'b-', 'LineWidth', 3, 'DisplayName', '||u|| (Constant Mass)');
plot(t_const, u_norm_var_interp, 'm--', 'LineWidth', 3, 'DisplayName', '||u|| (Variable Mass)');
xlabel('Time (yr)');
ylabel('Thrust Magnitude');
title('Comparison of Control Input Norms');
legend('Location','Best');

% LOCAL FUNCTIONS

% --- Constant Mass Dynamics and Shooting Function ---
function dydt = constantMass_odefun(~, y, mu, umax, rho)
    % y = [r (2); v (2); lambda_r (2); lambda_v (2)]
    r = y(1:2);
    v = y(3:4);
    lambda_r = y(5:6);
    lambda_v = y(7:8);
    norm_r = norm(r);
    norm_lv = norm(lambda_v);
    if norm_lv > 1e-8
        Gamma = umax/2*(1 + tanh((norm_lv - 1)/rho));
        u = - Gamma*(lambda_v/norm_lv);
    else
        u = [0; 0];
    end
    drdt = v;
    dvdt = -r/(norm_r^3) + u;
    dlam_v_dt = - lambda_r;
    A = (1/(norm_r^3))*eye(2) - 3*(r*r')/(norm_r^5);
    dlam_r_dt = A*lambda_v;
    dydt = [drdt; dvdt; dlam_r_dt; dlam_v_dt];
end

function F = constantMass_shootFunc(lambda0, x0, tf, mu, umax, rho, x_target)
    y0 = [x0; lambda0];
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [~, y] = ode45(@(t,y) constantMass_odefun(t,y,mu,umax,rho), [0 tf], y0, options);
    x_tf = y(end,1:4)';
    F = x_tf - x_target;
end

% --- Variable Mass Dynamics and Shooting Function ---
function dydt = variableMass_odefun(~, y, mu, umax, rho, Ispg0)
    % y = [r (2); v (2); m (1); lambda_r (2); lambda_v (2); lambda_m (1)]
    r = y(1:2);
    v = y(3:4);
    m = y(5);
    lambda_r = y(6:7);
    lambda_v = y(8:9);
    lambda_m = y(10);
    norm_r = norm(r);
    norm_lv = norm(lambda_v);
    if norm_lv > 1e-8
        S_val = (norm_lv/m) + (lambda_m/Ispg0);
        Gamma = umax/2*(1 + tanh((S_val - 1)/rho));
        u = - Gamma*(lambda_v/norm_lv);
    else
        u = [0; 0];
    end
    drdt = v;
    dvdt = -r/(norm_r^3) + (1/m)*u;
    dmdt = -(1/Ispg0)*norm(u);
    A = (1/(norm_r^3))*eye(2) - 3*(r*r')/(norm_r^5);
    dlam_r_dt = A*lambda_v;
    dlam_v_dt = - lambda_r;
    dlam_m_dt = lambda_v'*(1/m^2*u);
    dydt = [drdt; dvdt; dmdt; dlam_r_dt; dlam_v_dt; dlam_m_dt];
end

function F = variableMass_shootFunc(lambda0, x0, tf, mu, umax, rho, Ispg0, x_target)
    y0 = [x0; lambda0];
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [~, y] = ode45(@(t,y) variableMass_odefun(t,y,mu,umax,rho,Ispg0), [0 tf], y0, options);
    x_tf = y(end,1:4)';  % Compare only r and v
    F = x_tf - x_target;
end