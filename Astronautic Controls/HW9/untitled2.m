clear; close all; clc;

% Set flag for saturation of control input (true = apply saturation)
flag_saturate = true;
% Saturation threshold for control input in dimensional units (km/s^2);
% 1 m/s^2 = 0.001 km/s^2.
u_sat_dim = 0.001;

% Dimensional gravitational parameter [km^3/s^2]
mu_Earth = 3.9860e5;
day2sec  = 86400;

% Define the current (initial) orbital elements (from your table)
% x_slow = [a, e, i, Omega, omega] (angles in deg)
x_slow_current = [10000;    % a (km)
                  0.3;      % e
                  10;       % i (deg)
                  0;        % Omega (deg)
                  0];       % omega (deg)
M_current = 0;                % Mean anomaly at t=0 in deg
x_current = [x_slow_current; M_current];

% Define the target orbital elements:
x_slow_target = [60000;     % target a (km)
                 0.7;       % target e
                 130;       % target i (deg)
                 180;       % target Omega (deg)
                 270];      % target omega (deg)
M_target = 0;                 % target mean anomaly (deg)
x_target = [x_slow_target; M_target];

% Feedback gain on the slow variables: K = I_5
K = eye(5);

% Time span for propagation: 10 days (in seconds)
T_span = [0, 10*day2sec];

% ODE SOLUTION (tolerances set to 1E-12; simulation in dimensional units)
opts = odeset('AbsTol',1e-12,'RelTol',1e-12);
[tSol, xSol] = ode45(@(t,x) transferDynamics(t, x, x_target, mu_Earth, K, u_sat_dim, flag_saturate), T_span, x_current, opts);

% Extract orbital elements from solution:
aVals     = xSol(:,1);
eVals     = xSol(:,2);
iVals     = xSol(:,3);
OmegaVals = xSol(:,4);
omegaVals = xSol(:,5);
MVals     = xSol(:,6);

% Convert time to days:
tDays = tSol/day2sec;

figure('Color','w','Position',[50,50,1200,800]);
subplot(3,2,1); plot(tDays,aVals,'LineWidth',2); grid on;
    xlabel('Time (days)','FontSize',20); ylabel('a (km)','FontSize',20); title('Semi-major axis','FontSize',20);
subplot(3,2,2); plot(tDays,eVals,'LineWidth',2); grid on;
    xlabel('Time (days)','FontSize',20); ylabel('e','FontSize',20); title('Eccentricity','FontSize',20);
subplot(3,2,3); plot(tDays,iVals,'LineWidth',2); grid on;
    xlabel('Time (days)','FontSize',20); ylabel('i (deg)','FontSize',20); title('Inclination','FontSize',20);
subplot(3,2,4); plot(tDays,OmegaVals,'LineWidth',2); grid on;
    xlabel('Time (days)','FontSize',20); ylabel('\Omega (deg)','FontSize',20); title('RAAN','FontSize',20);
subplot(3,2,5); plot(tDays,omegaVals,'LineWidth',2); grid on;
    xlabel('Time (days)','FontSize',20); ylabel('\omega (deg)','FontSize',20); title('Argument of perigee','FontSize',20);
subplot(3,2,6); plot(tDays,MVals,'LineWidth',2); grid on;
    xlabel('Time (days)','FontSize',20); ylabel('M (deg)','FontSize',20); title('Mean anomaly','FontSize',20);

% ------------------ Recompute Control History ------------------------------
uSol = zeros(length(tSol),3);
for k = 1:length(tSol)
    xk = xSol(k,:)';
    f0k = f0Gauss_transfer(xk, mu_Earth);
    Bk = BmatrixGauss_transfer(xk, mu_Earth);
    B_slow = Bk(1:5,:);
    delta_x = xk(1:5) - x_target(1:5);
    u = -pinv(B_slow)*(K*delta_x);
    if flag_saturate && (norm(u) > u_sat_dim)
        u = u*(u_sat_dim/norm(u));
    end
    uSol(k,:) = u';
end

figure('Color','w','Position',[50,900,1200,800]);
plot(tDays, uSol*1000, 'LineWidth',2); % convert km/s^2 -> m/s^2
hold on;
yline(u_sat_dim*1000, 'k--', 'LineWidth',2);
xlabel('Time (days)','FontSize',20); ylabel('Control Input (m/s^2)','FontSize',20);
title('Control History','FontSize',20); legend('u_R','u_T','u_N','u_{max}','FontSize',20,'Location','best');
grid on;

% --------------------- 3D ORBIT PLOTS ---------------------------------------
% Propagate a high-resolution time vector for current and target orbits.
t_high = linspace(T_span(1), T_span(2), 5000);  % 5000-point time vector (in seconds)
x_high = interp1(tSol, xSol, t_high);

% Compute ECI positions for current orbit:
rECI = zeros(length(t_high),3);
for k = 1:length(t_high)
    xk = x_high(k,:)';
    rECI(k,:) = orbitalElements2ECI(xk(1:6), mu_Earth);
end

% Propagate target orbit using its two-body evolution:
% For target, a, e, i, Omega, omega are constant; only M changes.
a_target = x_target(1);
n_target = sqrt(mu_Earth/(a_target^3));  % rad/s (dimensional)
n_target_deg = n_target*180/pi;          % deg/s
rECI_target = zeros(length(t_high),3);
for k = 1:length(t_high)
    M_t = x_target(6) + n_target_deg*t_high(k);  % propagate target M
    x_target_prop = [x_target(1:5); M_t];
    rECI_target(k,:) = orbitalElements2ECI(x_target_prop, mu_Earth);
end

figure('Color','w','Position',[1300,50,1200,800]);
plot3(rECI(:,1), rECI(:,2), rECI(:,3), 'r', 'LineWidth',2); hold on;
plot3(rECI_target(:,1), rECI_target(:,2), rECI_target(:,3), 'b--', 'LineWidth',4);
xlabel('X (km)','FontSize',20); ylabel('Y (km)','FontSize',20); zlabel('Z (km)','FontSize',20);
title('3D Orbit Transfer: Current vs. Target Orbits','FontSize',20);
legend('Current Orbit','Target Orbit','FontSize',20,'Location','best');
grid on; axis equal;

% -------------------- δx PLOT (Error in Slow Orbital Elements) -------------
deltaX = xSol(:,1:5) - repmat(x_target(1:5)', length(tSol), 1);
deltaX_norm = sqrt(sum(deltaX.^2,2));
figure('Color','w','Position',[1300,900,1200,800]);
plot(tDays, deltaX(:,1), 'r', 'LineWidth',2); hold on;
plot(tDays, deltaX(:,2), 'g', 'LineWidth',2);
plot(tDays, deltaX(:,3), 'b', 'LineWidth',2);
plot(tDays, deltaX(:,4), 'm', 'LineWidth',2);
plot(tDays, deltaX(:,5), 'k', 'LineWidth',2);
plot(tDays, deltaX_norm, 'c--', 'LineWidth',2);
xlabel('Time (days)','FontSize',20); ylabel('\delta x components','FontSize',20);
title('\delta x (x_{slow} - x_{target}) vs. Time','FontSize',20);
legend('\delta a','\delta e','\delta i','\delta \Omega','\delta \omega','|\delta x|','FontSize',20,'Location','best');
grid on; set(gca,'FontSize',20);

% ------------------ 3D δr PLOT (Error in ECI Position) ---------------------
% Compute the 3D error: δr = rECI(current) - rECI(target)
delta_r = rECI - rECI_target;
figure('Color','w','Position',[2450,50,1200,800]);
plot3(delta_r(:,1), delta_r(:,2), delta_r(:,3), 'LineWidth',2); hold on;
% Mark initial and final error vectors:
plot3(delta_r(1,1), delta_r(1,2), delta_r(1,3), 'ko', 'MarkerSize',12, 'MarkerFaceColor', 'g'); 
plot3(delta_r(end,1), delta_r(end,2), delta_r(end,3), 'ks', 'MarkerSize',12, 'MarkerFaceColor', 'm');
xlabel('\delta r_1 (km)','FontSize',20); ylabel('\delta r_2 (km)','FontSize',20); zlabel('\delta r_3 (km)','FontSize',20);
title('3D \delta r: Error in ECI Position (Current - Target)','FontSize',20);
legend('Error Trajectory','Initial Error','Final Error','FontSize',20,'Location','best');
grid on; axis equal; set(gca,'FontSize',20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE Dynamics for Orbit Transfer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdot = transferDynamics(t, x, x_target, mu, K, u_sat, flag_sat)
    % x = [a; e; i; Omega; omega; M]
    f0 = f0Gauss_transfer(x, mu);
    B = BmatrixGauss_transfer(x, mu);
    delta_x = x(1:5) - x_target(1:5);
    B_slow = B(1:5,:);
    u = -pinv(B_slow)*(K*delta_x);
    if flag_sat && (norm(u) > u_sat)
        u = u*(u_sat/norm(u));
    end
    xdot = f0 + B*u;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f0Gauss_transfer: Natural rates (f0) for orbital elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f0 = f0Gauss_transfer(x, mu)
    a = x(1);
    n = sqrt(mu/a^3);  % rad/s, mean motion
    n_deg = n*180/pi;  % deg/s
    f0 = [0; 0; 0; 0; 0; n_deg];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BmatrixGauss_transfer: Mapping from thrust (RTN) to orbital element rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = BmatrixGauss_transfer(x, mu)
    a = x(1);
    e = x(2);
    i_deg = x(3);
    omega_deg = x(5);
    M_deg = x(6);
    % Clamp e to the physical range
    e = max(min(e, 0.99), eps);
    i_rad = deg2rad(i_deg);
    omega_rad = deg2rad(omega_deg);
    M_rad = deg2rad(M_deg);
    E = M_rad;
    for k = 1:50
        E_new = M_rad + e*sin(E);
        if abs(E_new - E) < 1e-12, break; end
        E = E_new;
    end
    nu = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
    r = a*(1 - e*cos(E));
    p = a*(1 - e^2);
    b = a*sqrt(1-e^2);
    argp = omega_rad;
    factor = 1 / sqrt(mu * p);
    row1 = [2*a^2*e*sin(nu),   2*a^2*p/r,      0];
    row2 = [p*sin(nu), (p + r)*cos(nu) + r*e,   0];
    row3 = [0, 0, r*cos(nu + argp)];
    if abs(sin(i_rad)) < 1e-12
        row4 = [0, 0, 0];
    else
        row4 = [0, 0, r*sin(nu + argp)/sin(i_rad)];
    end
    if abs(e) < 1e-12
        row5 = [0, 0, 0];
    else
        if abs(tan(i_rad)) < 1e-12
            row5 = [-p*cos(nu)/e, (p + r)*sin(nu)/e, 0];
        else
            row5 = [-p*cos(nu)/e, (p + r)*sin(nu)/e, -r*sin(nu+argp)/tan(i_rad)];
        end
    end
    row6 = [b*p*cos(nu)/(a*e)-2*b*r/a, -(b*(p + r)*sin(nu))/(a*e), 0];
    B = factor * [row1; row2; row3; row4; row5; row6];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% orbitalElements2ECI: Convert orbital elements to Cartesian ECI position.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rECI = orbitalElements2ECI(x, mu, L_scale)
    % x = [a, e, i, Omega, omega, M] with a in km and angles in deg.
    % In this simulation, we keep a dimensional so no non-dimensionalization here.
    a = x(1);
    e = x(2);
    % Clamp e to ensure real values
    e = max(min(e, 0.99), 1e-6);
    i_rad = deg2rad(x(3));
    Omega_rad = deg2rad(x(4));
    omega_rad = deg2rad(x(5));
    M_rad = deg2rad(x(6));
    E = M_rad;
    for k = 1:50
        E_new = M_rad + e*sin(E);
        if abs(E_new - E) < 1e-12, break; end
        E = E_new;
    end
    E = real(E); % ensure real E
    nu = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
    r = a*(1 - e*cos(E));
    r_pf = [r*cos(nu); r*sin(nu); 0];
    R = [ cos(Omega_rad)*cos(omega_rad)-sin(Omega_rad)*sin(omega_rad)*cos(i_rad), -cos(Omega_rad)*sin(omega_rad)-sin(Omega_rad)*cos(omega_rad)*cos(i_rad),  sin(Omega_rad)*sin(i_rad);
          sin(Omega_rad)*cos(omega_rad)+cos(Omega_rad)*sin(omega_rad)*cos(i_rad), -sin(Omega_rad)*sin(omega_rad)+cos(Omega_rad)*cos(omega_rad)*cos(i_rad), -cos(Omega_rad)*sin(i_rad);
          sin(omega_rad)*sin(i_rad),                                          cos(omega_rad)*sin(i_rad),                                          cos(i_rad)];
    rECI = R * r_pf;
end