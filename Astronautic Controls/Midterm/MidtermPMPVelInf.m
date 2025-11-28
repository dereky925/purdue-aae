clear; clc; close all;
format long

mu = 1;               % Non-dimensional gravitational parameter
umax = 0.1;           % Maximum control magnitude
tf = 8;               % Fixed final time (non-dimensional)
Ispg0 = 4000*9.8;

% vInf_dept
% Earth velocity roughly 30 km/s, so 0.1 for vinf
vInf_dept = 0.1; 

% Earth initial conditions (circular orbit at 1 AU)
% State: [rx; ry; vx; vy; m]
rEarth0 = [1; 0];
vEarth0 = [0; 1];  % Earth’s circular orbital velocity (nondim)
m0 = 1;

% Mars orbit parameters (circular orbit)
a_M = 1.524;                
nu_M0 = pi;                 
omega_M = sqrt(mu/(a_M^3));
marsState = @(t) [ a_M*cos(nu_M0 + omega_M*t);
                   a_M*sin(nu_M0 + omega_M*t);
                  -a_M*omega_M*sin(nu_M0 + omega_M*t);
                   a_M*omega_M*cos(nu_M0 + omega_M*t) ];
x_target = marsState(tf);

% Smoothing parameters for continuation
rhos = [1, 0.1, 1e-2, 1e-3];

% Shooting Setup
% Shooting unknowns are:
%   X = [deltaVx0; deltaVy0; lam_r0x; lam_r0y; lam_v0x; lam_v0y; lam_m0]
% where deltaV0 is the departure velocity offset so that
% v(0) = vEarth0 + deltaV0, with ||deltaV0|| = vInf_dept
lambda0_guess = [ 0.0;           % deltaVx0 (guess)
                  vInf_dept;     % deltaVy0 (guess to yield proper magnitude)
                  0.8; 0.1;      % lam_r0 (guess)
                  0.2; 1.0;      % lam_v0 (guess)
                  1.0];          % lam_m0 (guess)

options_fsolve = optimoptions('fsolve','Display','iter','TolFun',1e-8,'TolX',1e-8);

% Continuation in smoothing parameter ρ
for i = 1:length(rhos)
    rho = rhos(i);
    fun = @(X) shootingFunc_fixedVinf(X, rEarth0, vEarth0, m0, tf, mu, umax, rho, Ispg0, ...
                                      vInf_dept, x_target, marsState);
    [solX, fval, exitflag, output] = fsolve(fun, lambda0_guess, options_fsolve);
    
    % Compute initial velocity:
    % v(t0) = vEarth0 + deltaV0, where deltaV0 is the first two components of solX.
    v0 = vEarth0 + solX(1:2);
    v0_mag = norm(v0);
    % Launch azimuth angle (in degrees) computed in the inertial frame:
    launchAzimuth = atan2(v0(2), v0(1)) * 180/pi;
    fprintf('For ρ = %g: v(t0) = %.6f, Launch Azimuth = %.2f degrees\n', rho, v0_mag, launchAzimuth);
    % Report lambda_0 for this iteration:
    fprintf('For ρ = %g: lambda0 = %s\n', rho, mat2str(solX,6));
    
    lambda0_guess = solX;  % update guess for next iteration
end

% Propagate the Final Solution
rho_final = rhos(end);
[xf, Y, t_sol, u_sol, H_sol] = propagateSolution(solX, rEarth0, vEarth0, m0, tf, mu, umax, rho_final, Ispg0, marsState);

% Extract state and costate histories from Y
r_sol         = Y(:,1:2);
v_sol         = Y(:,3:4);
m_sol         = Y(:,5);
lambda_r_sol  = Y(:,6:7);
lambda_v_sol  = Y(:,8:9);
lambda_m_sol  = Y(:,10);

% Compute Switching Function & Thrust Profile History
S_sol = zeros(length(t_sol),1);
Gamma_sol = zeros(length(t_sol),1);
for k = 1:length(t_sol)
    norm_lv = norm(lambda_v_sol(k,:)');
    if norm_lv > 1e-12
        S_val = (norm_lv / m_sol(k)) + (lambda_m_sol(k)/Ispg0);
        S_sol(k) = S_val - 1;
        Gamma_sol(k) = umax/2 * (1 + tanh((S_val - 1)/rho_final));
    else
        S_sol(k) = -1;
        Gamma_sol(k) = 0;
    end
end

% Compute Hamiltonian Variation
H0 = H_sol(1);
H_diff = H_sol - H0;

% Compute Polar Coordinates for the Spacecraft State
xSc = r_sol(:,1);
ySc = r_sol(:,2);
rPolar = sqrt(xSc.^2 + ySc.^2);
theta = atan2(ySc, xSc);
vx = v_sol(:,1);
vy = v_sol(:,2);
vr = (xSc.*vx + ySc.*vy)./rPolar;
vtheta = (xSc.*vy - ySc.*vx)./rPolar;

% PLOTS

% Plot 0: Approximated Thrust vs. Switching Function (for different ρ)
s_vec = linspace(-1,1,200);
figure('Color','white','Position',[2000 0 1000 800]); hold on;
colors = lines(length(rhos));
for i = 1:length(rhos)
    rho_val = rhos(i);
    Gamma_vec = umax/2*(1 + tanh(s_vec./rho_val));
    plot(s_vec, Gamma_vec, 'Color', colors(i,:), 'LineWidth', 3, 'DisplayName', sprintf('ρ = %g', rho_val));
end
xlabel('Switching Function S');
ylabel('Approximated Thrust Magnitude \Gamma^*');
title('Thrust Profile vs. Switching Function');
legend('Location','Best'); grid on;
set(gca, 'FontSize', 20);

% Figure 1 - Trajectory
figure('Color','white','Position',[0 0 1000 1000]); hold on;
t_circle = linspace(0,2*pi,200);
% Plot Earth orbit (unit circle)
plot(cos(t_circle), sin(t_circle), 'k--', 'DisplayName', 'Earth Orbit');
% Plot Mars orbit
plot(a_M*cos(t_circle), a_M*sin(t_circle), 'r--', 'DisplayName', 'Mars Orbit');
% Mark Earth start
plot(rEarth0(1), rEarth0(2), 'ko','MarkerFaceColor','b','MarkerSize',20, 'DisplayName', 'Earth Start');
% Mark Mars start (using marsState at t = 0 for illustration)
xMars0 = marsState(0);
plot(xMars0(1), xMars0(2), 'ro','MarkerFaceColor','r','MarkerSize',20, 'DisplayName', 'Mars Start');
% Plot spacecraft trajectory
plot(r_sol(:,1), r_sol(:,2), 'b-', 'LineWidth', 3, 'DisplayName', 'Spacecraft Trajectory');
% Mark intercept (final position)
scatter(r_sol(end,1), r_sol(end,2), 500, 'bp','MarkerFaceColor','y','MarkerEdgeColor','k', 'DisplayName', 'Intercept');
axis equal; grid on;
xlabel('x (AU)'); ylabel('y (AU)');
title('Minimum-Fuel Trajectory from Earth to Mars (Variable Mass + Vel Direction)');
legend('Location','Best');
set(gca, 'FontSize', 20);

% Figure 2 - States vs. Time
figure('Color','white','Position',[1000 0 1000 1000]);
subplot(4,1,1); hold on; grid on;
title('States vs. Time');
ylabel('r (AU)');
plot(t_sol, rPolar, 'b-', 'LineWidth', 3, 'DisplayName', 'r');
legend('Location','Best'); set(gca, 'FontSize', 20);

subplot(4,1,2); hold on; grid on;
ylabel('\theta (rad)');
plot(t_sol, theta, 'b-', 'LineWidth', 3, 'DisplayName', '\theta');
legend('Location','Best'); set(gca, 'FontSize', 20);

subplot(4,1,3); hold on; grid on;
ylabel('Velocity (AU/yr)');
plot(t_sol, vr, 'b--', 'LineWidth', 3, 'DisplayName', 'v_r');
plot(t_sol, vtheta, 'b-', 'LineWidth', 3, 'DisplayName', 'v_\theta');
legend('Location','Best'); set(gca, 'FontSize', 20);

subplot(4,1,4); hold on; grid on;
xlabel('Time (yr)'); ylabel('Mass (nondim)');
plot(t_sol, m_sol*1500, 'r-', 'LineWidth', 3, 'DisplayName', 'm(t)');
legend('Location','Best'); set(gca, 'FontSize', 20);

% Figure 3 - Hamiltonian Difference vs. Time
figure('Color','white','Position',[0 1000 1000 1000]); hold on; grid on;
title('Hamiltonian Difference H(t) - H(0)');
xlabel('Time (yr)'); ylabel('H(t) - H(0)');
plot(t_sol, H_diff, 'b-', 'LineWidth', 3, 'DisplayName', 'H(t)-H(0)');
legend('Location','Best'); set(gca, 'FontSize', 20);

% Figure 4 - Switching Function and Thrust Profile vs. Time
figure('Color','white','Position',[1000 1000 1000 1000]); hold on; grid on;
plot(t_sol, S_sol, 'b-', 'LineWidth', 3, 'DisplayName', 'S = (||\lambda_v||/m + \lambda_m/Ispg0 - 1)');
plot(t_sol, Gamma_sol, 'r--', 'LineWidth', 3, 'DisplayName', '\\Gamma^*');
xlabel('Time (yr)');
ylabel('Value');
title('Switching Function and Thrust Profile vs. Time');
legend('Location','Best'); set(gca, 'FontSize', 20);

% Figure 5 - Control Time History 
figure('Color','white','Position',[1000 1000 1000 1000]); hold on; grid on;
plot(t_sol, u_sol(:,1), 'b-', 'LineWidth', 3, 'DisplayName', 'u_x');
plot(t_sol, u_sol(:,2), 'r-', 'LineWidth', 3, 'DisplayName', 'u_y');
plot(t_sol, vecnorm(u_sol,2,2), 'k-', 'LineWidth', 3, 'DisplayName', '||u||');
xlabel('Time (yr)');
ylabel('Control Acceleration');
title('Control Time History');
legend('Location','Best'); set(gca, 'FontSize', 20);

% Figure 6 - Costate Time History
lambda_r_radial = lambda_r_sol(:,1).*cos(theta) + lambda_r_sol(:,2).*sin(theta);
lambda_r_transverse = -lambda_r_sol(:,1).*sin(theta) + lambda_r_sol(:,2).*cos(theta);
lambda_v_radial = lambda_v_sol(:,1).*cos(theta) + lambda_v_sol(:,2).*sin(theta);
lambda_v_transverse = -lambda_v_sol(:,1).*sin(theta) + lambda_v_sol(:,2).*cos(theta);

figure('Color','white','Position',[2000 1000 1000 1000]);
subplot(3,1,1); hold on; grid on;
plot(t_sol, lambda_r_radial, 'b-', 'LineWidth', 3, 'DisplayName', '\lambda_{r,rad}');
plot(t_sol, lambda_r_transverse, 'r-', 'LineWidth', 3, 'DisplayName', '\lambda_{r,\theta}');
xlabel('Time (yr)'); ylabel('Costate (units)');
title('Costate \lambda_r');
legend('Location','Best'); set(gca, 'FontSize', 20);

subplot(3,1,2); hold on; grid on;
plot(t_sol, lambda_v_radial, 'b-', 'LineWidth', 3, 'DisplayName', '\lambda_{v,rad}');
plot(t_sol, lambda_v_transverse, 'r-', 'LineWidth', 3, 'DisplayName', '\lambda_{v,\theta}');
xlabel('Time (yr)'); ylabel('Costate (units)');
title('Costate \lambda_v');
legend('Location','Best'); set(gca, 'FontSize', 20);

subplot(3,1,3); hold on; grid on;
plot(t_sol, lambda_m_sol, 'm-', 'LineWidth', 3, 'DisplayName', '\lambda_m');
xlabel('Time (yr)'); ylabel('Mass Costate');
title('Costate \lambda_m');
legend('Location','Best'); set(gca, 'FontSize', 20);


function F = shootingFunc_fixedVinf(X, rEarth0, vEarth0, m0, tf, mu, umax, rho, Ispg0, ...
                                    vInf_dept, x_target, marsState)
    % X = [deltaVx0; deltaVy0; lam_r0x; lam_r0y; lam_v0x; lam_v0y; lam_m0]
    deltaV0 = X(1:2);       % Departure velocity offset
    lam_r0  = X(3:4);
    lam_v0  = X(5:6);
    lam_m0  = X(7);
    
    % Enforce initial state: r(0)=rEarth0, v(0)=vEarth0+deltaV0, m(0)=m0.
    x0 = [rEarth0; vEarth0 + deltaV0; m0];
    
    % Construct initial full state (state + costate)
    y0 = [x0; lam_r0; lam_v0; lam_m0];
    
    % Integrate ODEs
    options_ode = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [~, Y] = ode45(@(t,z) odefun(t,z,mu,umax,rho,Ispg0), [0 tf], y0, options_ode);
    
    % Final state
    xf = Y(end,1:5).';    % [rx; ry; vx; vy; m]
    lam_mf = Y(end, 10);   % final costate for mass
    
    % Boundary conditions
    xMars_tf = marsState(tf);
    F = zeros(6,1);
    F(1:2) = xf(1:2) - xMars_tf(1:2);    % Position error
    F(3:4) = xf(3:4) - xMars_tf(3:4);    % Velocity error
    F(5)   = lam_mf;                    % Final costate for mass condition
    F(6)   = norm(deltaV0) - vInf_dept;   % Departure speed constraint
end

function dydt = odefun(~, y, mu, umax, rho, Ispg0)
    % y = [r(2); v(2); m(1); lam_r(2); lam_v(2); lam_m(1)]
    r       = y(1:2);
    v       = y(3:4);
    m       = y(5);
    lam_r   = y(6:7);
    lam_v   = y(8:9);
    lam_m   = y(10);
    
    % Compute control
    norm_lam_v = norm(lam_v);
    if norm_lam_v > 0
        S_val = (norm_lam_v / m) + lam_m/Ispg0;
        Gamma = umax/2 * (1 + tanh((S_val - 1)/rho));
        u = -Gamma*(lam_v/norm_lam_v);
    else
        u = [0; 0];
    end
    
    % State dynamics
    rnorm = norm(r);
    drdt  = v;
    dvdt  = -r/(rnorm^3) + (1/m)*u;
    dmdt  = - (1/Ispg0)*norm(u);
    
    % Costate dynamics
    A = (1/(rnorm^3))*eye(2) - 3*(r*r')/(rnorm^5);
    dlam_rdt = A*lam_v;
    dlam_vdt = -lam_r;
    dlam_mdt = lam_v.'*((1/m^2)*u);
    
    dydt = [drdt; dvdt; dmdt; dlam_rdt; dlam_vdt; dlam_mdt];
end

function [xf, Y, t_sol, u_sol, H_sol] = ...
    propagateSolution(solX, rEarth0, vEarth0, m0, tf, mu, umax, rho, Ispg0, marsState)
    % Unpack solution vector
    deltaV0 = solX(1:2);
    lam_r0  = solX(3:4);
    lam_v0  = solX(5:6);
    lam_m0  = solX(7);
    
    % Set initial state with departure velocity offset
    x0 = [rEarth0; vEarth0 + deltaV0; m0];
    y0 = [x0; lam_r0; lam_v0; lam_m0];
    
    options_ode = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [t_sol, Y] = ode45(@(t,z) odefun(t,z,mu,umax,rho,Ispg0), [0 tf], y0, options_ode);
    xf = Y(end,1:5).';
    
    % Preallocate
    u_sol = zeros(length(t_sol),2);
    H_sol = zeros(length(t_sol),1);
    
    for k = 1:length(t_sol)
        r_     = Y(k,1:2).';
        v_     = Y(k,3:4).';
        m_     = Y(k,5);
        lam_r_ = Y(k,6:7).';
        lam_v_ = Y(k,8:9).';
        lam_m_ = Y(k,10);
        
        norm_lam_v = norm(lam_v_);
        if norm_lam_v > 1e-12
            S_val = (norm_lam_v/m_) + lam_m_/Ispg0;
            Gamma = umax/2*(1 + tanh((S_val - 1)/rho));
            u_ = -Gamma*(lam_v_/norm_lam_v);
        else
            u_ = [0;0];
        end
        u_sol(k,:) = u_.';
        
        % Compute Hamiltonian
        rnorm_ = norm(r_);
        H_sol(k) = lam_r_.'*v_ + lam_v_.'*( -r_/rnorm_^3 + (1/m_)*u_ ) ...
                    + lam_m_*( -1/Ispg0*norm(u_) ) + norm(u_);
    end
end