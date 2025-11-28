% =========================================================================
%
% Filename:       HW4_mass_included.m
% Author:         <Your Name>
% Institution:    <Your Institution>
% Course:         AAE590 - Applied Control in Astronautics
% Professor:      Dr. Kenshiro Oguri
% Assignment:     HW4 (Modified for mass)
% Semester:       Spring 2025
%
% Description:    Example code incorporating spacecraft mass as a state
%
% =========================================================================

clear; clc; close all;
format long

% PARAMETERS
mu = 1;                 % (Non-dimensional) gravitational parameter
umax = 0.1;             % Max THRUST magnitude (force units if you prefer) 
                        % NOTE: If you were previously using 0.1 as max acceleration,
                        %       you should adjust accordingly if 'u' is now a force.
tf = 8;                 % Fixed final time (non-dimensional)  <--- Changed to 8

% -- New mass-related parameters --
m0 = 1500;              % Initial mass [kg]
Isp = 4000;             % [sec]
g0 = 9.81;              % [m/s^2] standard gravity for mass flow

% Smoothing parameters for continuation
rhos = [1, 0.1, 1e-2, 1e-3];

% Earth initial conditions (circular orbit at 1 AU)
% State vector: [position (2); velocity (2)]
x0 = [1; 0; 0; 1];

% Mars orbit parameters (circular orbit)
a_M = 1.524;            % Mars semi-major axis (AU)
nu_M0 = pi;             % Mars initial true anomaly (radians)
omega_M = sqrt(mu/(a_M^3));  % Mars angular rate

% Anonymous function for Mars' state (position and velocity) at time t
marsState = @(t) [ a_M*cos(nu_M0 + omega_M*t);
                   a_M*sin(nu_M0 + omega_M*t);
                  -a_M*omega_M*sin(nu_M0 + omega_M*t);
                   a_M*omega_M*cos(nu_M0 + omega_M*t) ];

% Terminal condition: spacecraft must rendezvous with Mars at t = tf
x_target = marsState(tf);

% -------------------------------------------------------------------------
% PLOT: Approximated Thrust vs. Switching Function
% The switching function is S = ||lambda_v|| - 1,
% and the smoothed thrust is: Gamma = umax/2*(1 + tanh((||lambda_v|| - 1)/rho))
% -------------------------------------------------------------------------
s_vec = linspace(-1,1,200);
figure('Color','white','Position',[2000 0 1000 800]); hold on;
colors = lines(length(rhos));
for i = 1:length(rhos)
    rho_val = rhos(i);
    Gamma = umax/2*(1 + tanh(s_vec./rho_val));
    plot(s_vec, Gamma, 'Color', colors(i,:), 'LineWidth', 3, ...
         'DisplayName', sprintf('\\rho = %g', rho_val));
end
xlabel('Switching Function S = ||\\lambda_v|| - 1');
ylabel('Approximated Thrust Magnitude \\Gamma^*');
title('Thrust Profile vs. Switching Function');
legend('Location','Best'); grid on;
ax = gca; ax.FontSize = 20;

% -------------------------------------------------------------------------
% Solve the Two-Point Boundary Value Problem via Shooting
%
% *** IMPORTANT: Now we have 5 unknown costates, not 4. ***
% We must also enforce the transversality condition lambda_m(tf) = 0 if
% final mass is free (typical for minimum-fuel).
% -------------------------------------------------------------------------

% *** CHANGED: New dimension of lambda0_guess = 5 *** 
lambda0_guess = [4.83864; 1.72265; 2.08651; 5.38362; 100];  % Provided

options_fsolve = optimoptions('fsolve','Display','iter','TolFun',1e-8,'TolX',1e-8);

for i = 1:length(rhos)
    rho = rhos(i);
    shootingFunc = @(lambda0) shootFunc(lambda0, x0, m0, tf, mu, umax, rho, Isp, g0, x_target);
    [lambda0_sol, fval, exitflag, output] = fsolve(shootingFunc, lambda0_guess, options_fsolve);
    fprintf('For \\rho = %g, computed lambda0 = [%g; %g; %g; %g; %g]\n', ...
            rho, lambda0_sol(1), lambda0_sol(2), lambda0_sol(3), ...
            lambda0_sol(4), lambda0_sol(5));
    lambda0_guess = lambda0_sol;  % Update guess for next iteration
end

rho_final = rhos(end);

% -------------------------------------------------------------------------
% Integrate the Full ODE System (State + Costate)
% Full state vector y = [ r(2); v(2); m;  lambda_r(2); lambda_v(2); lambda_m ]
% dimension = 5 + 5 = 10
% -------------------------------------------------------------------------
y0_full = [x0; m0; lambda0_sol];  % (4 + 1) + 5 costates = 10
options_ode = odeset('RelTol',1e-12,'AbsTol',1e-12);

[t_sol, y_sol] = ode45(@(t,y) odefun(t,y,mu,umax,rho_final,Isp,g0), ...
                       [0 tf], y0_full, options_ode);

% Unpack the solution
r_sol       = y_sol(:,1:2);
v_sol       = y_sol(:,3:4);
m_sol       = y_sol(:,  5);
lambda_r_sol= y_sol(:,6:7);
lambda_v_sol= y_sol(:,8:9);
lambda_m_sol= y_sol(:,10);

% -------------------------------------------------------------------------
% Compute Control History
% same smoothing approach, but note that dv/dt = -mu*r/r^3 + (1/m)*u
% so 'u' is a force (units of N if you prefer)
% -------------------------------------------------------------------------
u_sol = zeros(length(t_sol),2);
for k = 1:length(t_sol)
    lam_v = lambda_v_sol(k,:)';
    norm_lam_v = norm(lam_v);
    if norm_lam_v > 1e-8
        Gamma = umax/2 * (1 + tanh((norm_lam_v - 1)/rho_final));
        % Direction is negative of lambda_v
        u_sol(k,:) = -Gamma * (lam_v / norm_lam_v)';
    else
        u_sol(k,:) = [0, 0];
    end
end

% -------------------------------------------------------------------------
% Compute Hamiltonian Variation Over Time (Optional)
% Now the Hamiltonian includes the m-dot term:
%   H = lambda_r^T v + lambda_v^T[ -mu*r/r^3 + (1/m)*u ] + lambda_m[ -||u||/(Isp*g0) ] + ...
% -------------------------------------------------------------------------
H_sol = zeros(length(t_sol),1);
for k = 1:length(t_sol)
    r = r_sol(k,:)';
    v = v_sol(k,:)';
    m = m_sol(k);
    lam_r = lambda_r_sol(k,:)';
    lam_v = lambda_v_sol(k,:)';
    lam_m = lambda_m_sol(k);
    u = u_sol(k,:)';

    norm_r = norm(r);
    norm_u = norm(u);

    % Hamiltonian:
    H_sol(k) = lam_r' * v ...
               + lam_v' * ( -r/(norm_r^3) + (1/m)*u ) ...
               + lam_m  * ( -norm_u/(Isp*g0) );
end
H0 = H_sol(1);
H_diff = H_sol - H0;

% -------------------------------------------------------------------------
% Convert to Polar Coordinates for plotting
% -------------------------------------------------------------------------
xSc = r_sol(:,1);
ySc = r_sol(:,2);
rPolar = sqrt(xSc.^2 + ySc.^2);
theta = atan2(ySc, xSc);

vx = v_sol(:,1);
vy = v_sol(:,2);
vr    = (xSc.*vx + ySc.*vy) ./ rPolar;
vtheta= (xSc.*vy - ySc.*vx)./ rPolar;

% -------------------------------------------------------------------------
% FIGURE 1: Trajectories in Cartesian (Optimal Transfer)
% -------------------------------------------------------------------------
figure('Color','white','Position',[0 0 1000 1000]); hold on;
r0 = 1;       % Earth orbit radius
rMars = a_M;  % Mars orbit radius
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
title('Minimum-Fuel Trajectory from Earth to Mars (with Mass State)');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

% -------------------------------------------------------------------------
% FIGURE 2: States vs. Time (3 Subplots) + (Mass subplot if desired)
% -------------------------------------------------------------------------
figure('Color','white','Position',[1000 0 1000 1000]);

% Subplot 1: Radial distance vs. time
subplot(4,1,1); hold on; grid on;
title('States vs. Time');
ylabel('r (AU)');
plot(t_sol, rPolar, 'b-', 'LineWidth', 3, 'DisplayName', 'r');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

% Subplot 2: Radial and Tangential Velocities vs. time
subplot(4,1,2); hold on; grid on;
ylabel('Velocity (AU/yr)');
plot(t_sol, vr, 'b--', 'LineWidth', 3, 'DisplayName', 'v_r');
plot(t_sol, vtheta, 'b-', 'LineWidth', 3, 'DisplayName', 'v_\theta');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

% Subplot 3: Angular Position vs. time
subplot(4,1,3); hold on; grid on;
ylabel('\theta (rad)');
plot(t_sol, theta, 'b-', 'LineWidth', 3, 'DisplayName', '\theta');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

% Subplot 4: Mass vs. time
subplot(4,1,4); hold on; grid on;
xlabel('Time (yr)'); ylabel('Mass (kg)');
plot(t_sol, m_sol, 'r-', 'LineWidth', 3, 'DisplayName', 'm');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

% -------------------------------------------------------------------------
% FIGURE 3: Hamiltonian Difference vs. Time
% -------------------------------------------------------------------------
figure('Color','white','Position',[0 1000 1000 1000]); hold on; grid on;
title('H(t) - H(0)');
xlabel('Time (yr)'); ylabel('H(t) - H(0)');
plot(t_sol, H_diff, 'b-', 'LineWidth', 3, 'DisplayName', 'H(t)-H(0)');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

% -------------------------------------------------------------------------
% FIGURE 4: Switching Function and Thrust Profile vs. Time
% S = ||lambda_v|| - 1,  Gamma = umax/2*(1 + tanh((||lambda_v|| - 1)/rho))
% -------------------------------------------------------------------------
S_sol = zeros(length(t_sol),1);
Gamma_sol = zeros(length(t_sol),1);
for k = 1:length(t_sol)
    norm_lv = norm(lambda_v_sol(k,:)');
    S_sol(k) = norm_lv - 1;
    Gamma_sol(k) = umax/2*(1 + tanh((norm_lv - 1)/rho_final));
end
figure('Color','white','Position',[1000 1000 1000 1000]); hold on; grid on;
plot(t_sol, S_sol, 'b-', 'LineWidth', 3, ...
     'DisplayName', 'Switching Function S = ||\\lambda_v|| - 1');
plot(t_sol, Gamma_sol, 'r--', 'LineWidth', 3, ...
     'DisplayName', 'Thrust Magnitude \\Gamma^*');
xlabel('Time (yr)');
ylabel('Value');
title('Switching Function and Thrust Profile vs. Time');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

% -------------------------------------------------------------------------
% FIGURE 5: Control Time History
% -------------------------------------------------------------------------
figure('Color','white','Position',[1000 1000 1000 1000]); hold on; grid on;
plot(t_sol, u_sol(:,1), 'b-', 'LineWidth', 3, 'DisplayName', 'u_x');
plot(t_sol, u_sol(:,2), 'r-', 'LineWidth', 3, 'DisplayName', 'u_y');
plot(t_sol, vecnorm(u_sol,2,2), 'k-', 'LineWidth', 3, 'DisplayName', 'u_{mag}');
xlabel('Time (yr)');
ylabel('Control Force');
title('Control Time History');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

% -------------------------------------------------------------------------
% FIGURE 6: Costate Time History
% -------------------------------------------------------------------------
figure('Color','white','Position',[2000 1000 1000 1000]);
subplot(3,1,1); hold on; grid on;
plot(t_sol, lambda_r_sol(:,1), 'b-', 'LineWidth', 3, 'DisplayName', '\lambda_{r_x}');
plot(t_sol, lambda_r_sol(:,2), 'r-', 'LineWidth', 3, 'DisplayName', '\lambda_{r_y}');
xlabel('Time (yr)');
ylabel('Position Costates');
title('Costate Time History: \lambda_r');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

subplot(3,1,2); hold on; grid on;
plot(t_sol, lambda_v_sol(:,1), 'b-', 'LineWidth', 3, 'DisplayName', '\lambda_{v_x}');
plot(t_sol, lambda_v_sol(:,2), 'r-', 'LineWidth', 3, 'DisplayName', '\lambda_{v_y}');
xlabel('Time (yr)');
ylabel('Velocity Costates');
title('Costate Time History: \lambda_v');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

subplot(3,1,3); hold on; grid on;
plot(t_sol, lambda_m_sol, 'm-', 'LineWidth', 3, 'DisplayName', '\lambda_m');
xlabel('Time (yr)');
ylabel('Mass Costate');
title('Costate Time History: \lambda_m');
legend('Location','Best');
ax = gca; ax.FontSize = 20;


% =========================================================================
% NESTED FUNCTION: ODE dynamics for (r,v,m,lambda_r,lambda_v,lambda_m)
% =========================================================================
function dydt = odefun(t,y,mu,umax,rho,Isp,g0)
    %
    % y = [r(2); v(2); m;  lambda_r(2); lambda_v(2); lambda_m]
    %
    r        = y(1:2);
    v        = y(3:4);
    m        = y(5);
    lambda_r = y(6:7);
    lambda_v = y(8:9);
    lambda_m = y(10);

    norm_r = norm(r);
    norm_lam_v = norm(lambda_v);

    % Smoothing-based thrust magnitude
    if norm_lam_v > 1e-8
        % Gamma = umax/2 * (1 + tanh((norm_lam_v - 1)/rho));
        % u_dir = - (lambda_v / norm_lam_v);  % direction
        % u     = Gamma * u_dir;             % "u" is the force vector

        % S = (||lambda_v|| / m) + (lambda_m / (Isp*g0))
        S_val = (norm_lam_v / m) + (lambda_m/(Isp*g0));
        Gamma = umax/2 * (1 + tanh((S_val - 1)/rho));  % Force magnitude
        % Direction is opposite lambda_v
        u = -Gamma * (lambda_v / norm_lam_v);

    else
        u     = [0;0];
    end

    % ---------------------
    % State derivatives
    % ---------------------
    drdt    = v;
    dvdt    = - (mu/(norm_r^3)) * r + (1/m)*u;
    dmdt    = - (1/(Isp*g0)) * norm(u);

    % ---------------------
    % Costate derivatives
    %   dot(lambda_r) = -∂H/∂r
    %   dot(lambda_v) = -∂H/∂v
    %   dot(lambda_m) = -∂H/∂m
    % ---------------------
    % 1) lambda_r-dot
    %    The only r-dependence in the Hamiltonian is inside -mu*r/||r||^3
    %    as part of v-dot.  The same "A*lambda_v" form as original:
    A = (1/(norm_r^3))*eye(2) - 3*(r*r')/(norm_r^5);
    dlambda_r_dt = A * lambda_v;

    % 2) lambda_v-dot
    %    The only v-dependence is in (lambda_r^T * v). => partial wrt v => lambda_r
    dlambda_v_dt = - lambda_r;

    % 3) lambda_m-dot
    %    H = lambda_r^T v + lambda_v^T[ -mu*r/r^3 + (1/m)*u ] + lambda_m[ -||u||/(Isp*g0) ]
    %    The only m-dependence in the 2nd bracket: (1/m)*u => derivative wrt m is -1/m^2 * u
    %    => dot(lambda_m) = - dH/dm = + (lambda_v^T)*(1/m^2 * u)
    dlambda_m_dt = (1/m^2) * (lambda_v' * u);

    dydt = [drdt; dvdt; dmdt; dlambda_r_dt; dlambda_v_dt; dlambda_m_dt];
end

% =========================================================================
% NESTED FUNCTION: Shooting Function
% Now we have 5 unknown costates => we have 5 boundary conditions:
%   1-4) r(tf) - r_target, v(tf) - v_target
%   5)   lambda_m(tf) = 0  (free final mass -> transversality condition)
% =========================================================================
function F = shootFunc(lambda0, x0, m0, tf, mu, umax, rho, Isp, g0, x_target)
    %
    % lambda0 is 5x1 => [lambda_r(2); lambda_v(2); lambda_m]
    % x0 is 4x1 => [r(2); v(2)]
    %
    y0 = [x0; m0; lambda0];  % dimension 4 + 1 + 5 = 10
    options_ode_inner = odeset('RelTol',1e-12,'AbsTol',1e-12);

    [~, y] = ode45(@(t,Y) odefun(t,Y,mu,umax,rho,Isp,g0), [0 tf], y0, options_ode_inner);

    % final states
    x_tf = y(end,1:4)';    % r(tf), v(tf)
    lambda_m_tf = y(end,10);

    % We want: r(tf)=x_target(1:2), v(tf)=x_target(3:4), and lambda_m(tf)=0
    F = [ x_tf - x_target;   % 4 conditions
          lambda_m_tf ];    % 5th condition
end