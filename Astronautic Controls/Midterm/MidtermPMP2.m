clear; clc; close all;
format long

% PARAMETERS
mu = 1;               % Non-dimensional gravitational parameter
umax = 0.1;           % Maximum control magnitude
tf = 8;               % Fixed final time (non-dimensional)
Ispg0 = 4000*9.8;     % Given

% Smoothing parameters
rhos = [1, 0.1, 1e-2, 1e-3];

% Earth initial conditions (circular orbit at 1 AU)
% Now the state vector is 5D: [rx, ry, vx, vy, m].
% Assume initial mass = 1 in nondimensional units
x0 = [1; 0; 0; 1; 1];

% Mars orbit parameters (circular orbit)
a_M = 1.524;                % Mars semi-major axis (AU)
nu_M0 = pi;                 % Mars initial true anomaly (radians)
omega_M = sqrt(mu/(a_M^3)); % Mars angular rate

% Anonymous function for Mars' state (position and velocity) at time t
marsState = @(t) [ a_M*cos(nu_M0 + omega_M*t);
                   a_M*sin(nu_M0 + omega_M*t);
                  -a_M*omega_M*sin(nu_M0 + omega_M*t);
                   a_M*omega_M*cos(nu_M0 + omega_M*t) ];

% Terminal condition: spacecraft must rendezvous with Mars at tf.
x_target = marsState(tf);

% Plot Approximated Thrust vs. Switching Function
s_vec = linspace(-1,1,200);
figure('Color','white','Position',[2000 0 1000 800]); hold on;
colors = lines(length(rhos));
for i = 1:length(rhos)
    rho_val = rhos(i);
    Gamma = umax/2*(1 + tanh(s_vec./rho_val));
    plot(s_vec, Gamma, 'Color', colors(i,:), 'LineWidth', 3, ...
         'DisplayName', sprintf('\\rho = %g', rho_val));
end
xlabel('Switching Function S = ||\lambda_v|| - 1');
ylabel('Approximated Thrust Magnitude \Gamma^*');
title('Thrust Profile vs. Switching Function');
legend('Location','Best'); grid on;
ax = gca; ax.FontSize = 20;

% Solve the Two-Point Boundary Value Problem via Shooting
lambda0_guess = [0.849671; 0.101611; 0.201818; 1.08465; 1];

options_fsolve = optimoptions('fsolve','Display','iter','TolFun',1e-8,'TolX',1e-8);
for i = 1:length(rhos)
    rho = rhos(i);
    shootingFunc = @(lambda0) shootFunc(lambda0, x0, tf, mu, umax, rho, Ispg0, x_target);
    [lambda0_sol, fval, exitflag, output] = fsolve(shootingFunc, lambda0_guess, options_fsolve);
    fprintf('For \rho = %g, computed lambda0 = [%g; %g; %g; %g; %g]\n', ...
            rho, lambda0_sol(1), lambda0_sol(2), lambda0_sol(3), ...
            lambda0_sol(4), lambda0_sol(5));
    lambda0_guess = lambda0_sol;  % Update guess for next iteration
end
rho_final = rhos(end);

% Integrate the Full ODE System (State and Costate)
% Full state: y = [rx, ry, vx, vy, m, lam_rx, lam_ry, lam_vx, lam_vy, lam_m].
y0 = [x0; lambda0_sol];
options_ode = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_sol, y_sol] = ode45(@(t,y) odefun(t,y,mu,umax,rho_final,Ispg0), [0 tf], y0, options_ode);

% Extract state and costate
r_sol       = y_sol(:,1:2);
v_sol       = y_sol(:,3:4);
m_sol       = y_sol(:,5);
lambda_r_sol= y_sol(:,6:7);
lambda_v_sol= y_sol(:,8:9);
lambda_m_sol= y_sol(:,10);

% Compute Control History
u_sol = zeros(length(t_sol),2);
for k = 1:length(t_sol)
    lam_v = lambda_v_sol(k,:)';
    norm_lam_v = norm(lam_v);
    if norm_lam_v > 0
        % Updated switching function
        S_val = (norm_lam_v / m_sol(k)) + (lambda_m_sol(k)/Ispg0);
        Gamma = umax/2 * (1 + tanh((S_val - 1)/rho_final));  % Force magnitude
        u_sol(k,:) = - Gamma * (lam_v/norm_lam_v);
    else
        u_sol(k,:) = [0 0];
    end
end

% Compute Hamiltonian Variation Over Time
H_sol = zeros(length(t_sol),1);
for k = 1:length(t_sol)
    r_  = r_sol(k,:)';
    v_  = v_sol(k,:)';
    m_  = m_sol(k);
    lam_r= lambda_r_sol(k,:)';
    lam_v= lambda_v_sol(k,:)';
    lam_m= lambda_m_sol(k);
    norm_r = norm(r_);
    u_  = u_sol(k,:)';
    H_sol(k) = lam_r'*v_ ...
             + lam_v'*(-r_/norm_r^3 + (1/m_)*u_) ...
             + lam_m * ( -1/Ispg0 * norm(u_) ) ...
             + norm(u_);
end
H0 = H_sol(1);
H_diff = H_sol - H0;

% Compute Polar Coordinates for State (for plotting)
xSc = r_sol(:,1);
ySc = r_sol(:,2);
rPolar = sqrt(xSc.^2 + ySc.^2);
theta = atan2(ySc, xSc);
vx = v_sol(:,1);
vy = v_sol(:,2);
vr = (xSc.*vx + ySc.*vy)./rPolar;
vtheta = (xSc.*vy - ySc.*vx)./rPolar;

% Figure 1 - Trajectories in Cartesian 
figure('Color','white','Position',[0 0 1000 1000]); hold on;
t_circle = linspace(0,2*pi,200);
% Plot Earth orbit
plot(cos(t_circle), sin(t_circle), 'k--', 'DisplayName','Earth Orbit');
% Plot Mars orbit
plot(a_M*cos(t_circle), a_M*sin(t_circle), 'r--', 'DisplayName','Mars Orbit');
% Mark Earth start (blue marker)
plot(1, 0, 'ko','MarkerFaceColor','b','MarkerSize',20, 'DisplayName','Earth Start');
% Mark Mars start (red marker at pi)
plot(a_M*cos(pi), a_M*sin(pi), 'ro','MarkerFaceColor','r','MarkerSize',20, 'DisplayName','Mars Start');
% Plot spacecraft trajectory
plot(xSc, ySc, 'b-', 'LineWidth', 3, 'DisplayName', 'Spacecraft Trajectory');
% Mark intercept
scatter(xSc(end), ySc(end), 1200, 'bp', 'MarkerFaceColor','y','MarkerEdgeColor','k', 'DisplayName', 'Intercept');
axis equal; grid on;
xlabel('x (AU)'); ylabel('y (AU)');
title('Minimum-Fuel Trajectory from Earth to Mars (Variable Mass)');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

% Figure 2 - States vs. Time
figure('Color','white','Position',[1000 0 1000 1000]);
subplot(4,1,1); hold on; grid on;
title('States vs. Time');
ylabel('r (AU)');
plot(t_sol, rPolar, 'b-', 'LineWidth', 3, 'DisplayName', 'r');
legend('Location','Best'); ax = gca; ax.FontSize = 20;

subplot(4,1,2); hold on; grid on;
ylabel('Velocity (AU/yr)');
plot(t_sol, vr, 'b--', 'LineWidth', 3, 'DisplayName', 'v_r');
plot(t_sol, vtheta, 'b-', 'LineWidth', 3, 'DisplayName', 'v_\theta');
legend('Location','Best'); ax = gca; ax.FontSize = 20;

subplot(4,1,3); hold on; grid on;
ylabel('\theta (rad)');
plot(t_sol, theta, 'b-', 'LineWidth', 3, 'DisplayName', '\theta');
legend('Location','Best'); ax = gca; ax.FontSize = 20;

subplot(4,1,4); hold on; grid on;
xlabel('Time (yr)'); ylabel('Mass (nondim)');
plot(t_sol, m_sol*1500, 'r-', 'LineWidth', 3, 'DisplayName', 'm(t)');
legend('Location','Best'); ax = gca; ax.FontSize = 20;

% Figure 3 - Hamiltonian Difference vs. Time
figure('Color','white','Position',[0 1000 1000 1000]); hold on; grid on;
title('H(t) - H(0)');
xlabel('Time (yr)'); ylabel('H(t) - H(0)');
plot(t_sol, H_diff, 'b-', 'LineWidth', 3, 'DisplayName', 'H(t)-H(0)');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

% Switching Function and Thrust Profile vs. Time
S_sol = zeros(length(t_sol),1);
Gamma_sol = zeros(length(t_sol),1);
for k = 1:length(t_sol)
    norm_lv = norm(lambda_v_sol(k,:)');
    % Updated switching function 
    S_val = (norm_lv / m_sol(k)) + (lambda_m_sol(k)/Ispg0);
    S_sol(k) = S_val - 1;
    Gamma_sol(k) = umax/2 * (1 + tanh((S_val - 1)/rho_final));
end
% Figure 4 - Thrust profile
figure('Color','white','Position',[1000 1000 1000 1000]); hold on; grid on;
plot(t_sol, S_sol, 'b-', 'LineWidth', 3, 'DisplayName', 'S = (||\lambda_v||/m + \lambda_m/Ispg0 - 1)');
plot(t_sol, Gamma_sol, 'r--', 'LineWidth', 3, 'DisplayName', '\\Gamma^*');
xlabel('Time (yr)');
ylabel('Value');
title('Switching Function and Thrust Profile vs. Time');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

% Figure 5 - Control Time History
figure('Color','white','Position',[1000 1000 1000 1000]); hold on; grid on;
plot(t_sol, u_sol(:,1), 'b-', 'LineWidth', 3, 'DisplayName', 'u_x');
plot(t_sol, u_sol(:,2), 'r-', 'LineWidth', 3, 'DisplayName', 'u_y');
plot(t_sol, vecnorm(u_sol,2,2), 'k-', 'LineWidth', 3, 'DisplayName', '||u||');
xlabel('Time (yr)');
ylabel('Control Acceleration');
title('Control Time History');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

% Figure 6 - Costate Time History
lambda_r_radial = lambda_r_sol(:,1).*cos(theta) + lambda_r_sol(:,2).*sin(theta);
lambda_r_transverse = -lambda_r_sol(:,1).*sin(theta) + lambda_r_sol(:,2).*cos(theta);
% For lambda_v:
lambda_v_radial = lambda_v_sol(:,1).*cos(theta) + lambda_v_sol(:,2).*sin(theta);
lambda_v_transverse = -lambda_v_sol(:,1).*sin(theta) + lambda_v_sol(:,2).*cos(theta);

figure('Color','white','Position',[2000 1000 1000 1000]);

subplot(3,1,1); hold on; grid on;
plot(t_sol, lambda_r_radial, 'b-', 'LineWidth', 3, 'DisplayName', '\lambda_{r,rad}');
plot(t_sol, lambda_r_transverse, 'r-', 'LineWidth', 3, 'DisplayName', '\lambda_{r,\theta}');
xlabel('Time [years]'); ylabel('Costate (units)');
title('Costate \lambda_r');
legend('Location','Best'); ax = gca; ax.FontSize = 20;

subplot(3,1,2); hold on; grid on;
plot(t_sol, lambda_v_radial, 'b-', 'LineWidth', 3, 'DisplayName', '\lambda_{v,rad}');
plot(t_sol, lambda_v_transverse, 'r-', 'LineWidth', 3, 'DisplayName', '\lambda_{v,\theta}');
xlabel('Time [years]'); ylabel('Costate (units)');
title('Costate \lambda_v');
legend('Location','Best'); ax = gca; ax.FontSize = 20;

subplot(3,1,3); hold on; grid on;
plot(t_sol, lambda_m_sol, 'm-', 'LineWidth', 3, 'DisplayName', '\lambda_m');
xlabel('Time [years]'); ylabel('Mass Costate');
title('Costate \lambda_m');
legend('Location','Best'); ax = gca; ax.FontSize = 20;

% ODE Dynamics for State and Costate
function dydt = odefun(~, y, mu, umax, rho, Ispg0)
    % y = [r(2); v(2); m(1); lambda_r(2); lambda_v(2); lambda_m(1)]
    r       = y(1:2);
    v       = y(3:4);
    m       = y(5);
    lam_r   = y(6:7);
    lam_v   = y(8:9);
    lam_m   = y(10);

    norm_lam_v = norm(lam_v);
    if norm_lam_v > 0
        S_val = (norm_lam_v / m) + (lam_m/(Ispg0));
        Gamma = umax/2 * (1 + tanh((S_val - 1)/rho));  % Force magnitude
        u = -Gamma*(lam_v/norm_lam_v);
    else
        u = [0; 0];
    end

    % State equations
    norm_r = norm(r);
    drdt = v;
    dvdt = -r/(norm_r^3) + (1/m)*u;
    dmdt = -(1/Ispg0)*norm(u);

    % Costate equations
    A = (1/(norm_r^3))*eye(2) - 3*(r*r')/(norm_r^5);
    dlam_rdt = A*lam_v;
    dlam_vdt = - lam_r;
    dlam_mdt = lam_v'*(1/m^2 * u);

    dydt = [drdt; dvdt; dmdt; dlam_rdt; dlam_vdt; dlam_mdt];
end

function F = shootFunc(lambda0, x0, tf, mu, umax, rho, Ispg0, x_target)
    % x0 is [r(2), v(2), m], lambda0 is [lam_r(2), lam_v(2), lam_m]
    y0 = [x0; lambda0];
    options_ode_inner = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [~, y] = ode45(@(t,z) odefun(t,z,mu,umax,rho,Ispg0), [0 tf], y0, options_ode_inner);

    % Final state at tf
    x_tf = y(end,1:4)';  % Only compare r and v
    F = x_tf - x_target; % final mass is free => do not constrain y(end,5)
end