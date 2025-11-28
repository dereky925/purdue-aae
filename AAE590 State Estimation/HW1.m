% =========================================================================
% 
% Filename:       HW1.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE590 - Fundamentals of State Estimation
% Professor:      Dr. Keith LeGrand
% Contact:        klegrand@purdue.edu
% Assignment:     HW 1
% Semester:       Fall 2025
% 
% Description:
%
%
%
% =========================================================================

%% 1b
set(groot, 'DefaultAxesFontSize', 20);   
set(groot, 'DefaultLegendFontSize', 20); 

clear;clc;close all


syms t Omega real
A = [ 0  0   1     0
      0  0   0     1
      0  0   0    Omega
      0  0   -Omega  0  ];

Phi = expm(A*t);        
disp('Phi(t):'); pretty(Phi)

%% 
clear;clc;close all

% Speed of jet
v0_mag   = 230;  
x0       = [0; 0; v0_mag; 0];

t1       = 120;    % s - begin turn
t2       = 240;    % s - end turn
t3       = 440;    % s - to the end of sim

Omega_deg_per_min = 60;         % given
Omega_rps = deg2rad(Omega_deg_per_min/60);  % rad/s  - Converted

dt = 0.1;

% Create A function to pass in omegas
A_func = @(Om) [0 0 1 0;
               0 0 0 1;
               0 0 0 -Om;
               0 0 Om 0];

% STM
Phi = @(Om, t) expm(A_func(Om)*t);

% Pre-allocate
t = 0:dt:t3;
N = numel(t);
x = zeros(N,4);                       

% Initial condition at t=0 
x(1,:) = [0, 0, v0_mag, 0];


for k = 1:N-1
    tk   = t(k);
    tk1  = t(k+1);
    del_t  = tk1 - tk;

    % Piecewise-constant Omega depending on time
    if     tk < t1
        Om = 0;                          % straight
    elseif tk < t2
        Om = Omega_rps;            % coordinated turn (CCW)
    else
        Om = 0;                          % straight again
    end

    x(k+1,:)= (Phi(Om, del_t) * x(k,:)')';
end

% Extract
xi  = x(:,1);
eta = x(:,2);
u   = x(:,3);
v   = x(:,4);


figure('Position',[500 500 1300 1000],'Color','white'); 
hold on; grid minor; axis equal
plot(xi, eta, 'LineWidth', 2);
xlabel('\xi [m]'); ylabel('\eta [m]');
title('Coordinated-Turn Trajectory');

% Get indices nearest to key times
[~,i1] = min(abs(t - t1));
[~,i2] = min(abs(t - t2));
[~,i3] = min(abs(t - (t2 + t3)));

% Markers for segment boundaries
scatter(xi(1),  eta(1), 150, 'ko', 'MarkerFaceColor','k');   text(xi(1),eta(1),'  t=0', 'FontSize', 20);
scatter(xi(i1), eta(i1), 150,'ks', 'MarkerFaceColor','y');   text(xi(i1),eta(i1),'  t=120 s', 'FontSize', 20);
scatter(xi(i2), eta(i2), 150, 'kd', 'MarkerFaceColor','c');   text(xi(i2),eta(i2),'  t=240 s', 'FontSize', 20);
scatter(xi(i3), eta(i3), 150,'k^', 'MarkerFaceColor','g');   text(xi(i3),eta(i3),'  t=440 s', 'FontSize', 20);

legend({'Trajectory','t=0','t=120 s','t=240 s','t=440 s'}, 'Location','bestoutside');




% Plot over NYC

lat0_deg = 40.7128;      % [deg]
lon0_deg = -74.0060;     % [deg]

R_earth = 6378137;       % [m] WGS-84 mean equatorial radius

% Small-angle ENU -> geodetic approximation
lat_deg = lat0_deg + (eta ./ R_earth) * 180/pi;
lon_deg = lon0_deg + (xi  ./ (R_earth * cosd(lat0_deg))) * 180/pi;

% Geographic plot (uses base MATLAB's geographic axes)
figure('Position',[800 500 1300 1000],'Color','white'); 
geoplot(lat_deg, lon_deg, 'LineWidth', 2); hold on
geoplot(lat_deg(1),  lon_deg(1),  'o', 'MarkerSize', 15, 'MarkerFaceColor','k', 'MarkerEdgeColor','k')
geoplot(lat_deg(i1), lon_deg(i1), 's', 'MarkerSize', 15,'MarkerFaceColor','y', 'MarkerEdgeColor','k')
geoplot(lat_deg(i2), lon_deg(i2), 'd', 'MarkerSize', 15,'MarkerFaceColor','c', 'MarkerEdgeColor','k')
geoplot(lat_deg(i3), lon_deg(i3), '^', 'MarkerSize', 15,'MarkerFaceColor','g', 'MarkerEdgeColor','k')
text(lat_deg(1),  lon_deg(1),  '  t=0', 'FontSize', 20)
text(lat_deg(i1), lon_deg(i1), '  t=120 s', 'FontSize', 20)
text(lat_deg(i2), lon_deg(i2), '  t=240 s', 'FontSize', 20)
text(lat_deg(i3), lon_deg(i3), '  t=440 s', 'FontSize', 20)

title('Trajectory in Latitude/Longitude', 'FontSize', 25)
geobasemap streets

gx = gca;       
gx.FontSize = 20;  % increase tick label size


%% Problem 4

O = [0 1; 0 1];
rank(O)

%%

clear; clc; close all

par.rho0 = 3.4e-3;   % [lb·s^2/ft^4]
par.g    = 32.2;     % [ft/s^2]
par.kp   = 22000;    % [ft]

x0_nom = [60000; 0; 2000];  % [ft; ft/s; lb/ft^2]
tspan = [0 200];    

opts_event = odeset('RelTol',1e-10,'AbsTol',1e-10, ...
                    'Events',@(t,x) impactEvent(t,x));

% Nonlinear simulation to impact (nominal ballistic coeff, beta)
[t_nom, x_nom] = ode45(@(t,x) F_nonlinear(t,x,par), tspan, x0_nom, opts_event);
alt_nom = x_nom(:,1);
vel_nom = x_nom(:,2);
tf_nom  = t_nom(end);        % nominal impact time
x_nom_tf = x_nom(end,:).';   % nominal state at impact (x1 ~ 0)

% Part (a) altitude & velocity
fig1 = figure('Position',[80 100 1200 800],'Color','w');
tiledlayout(fig1,2,1,'TileSpacing','compact','Padding','compact');

nexttile; hold on; grid minor;
plot(t_nom, alt_nom, 'LineWidth', 1.8);
xlabel('Time, t [s]'); ylabel('Altitude, x_1 [ft]');
title('Altitude vs Time (ODE45)');

nexttile; hold on; grid minor;
plot(t_nom, vel_nom, 'LineWidth', 1.8);
xlabel('Time, t [s]'); ylabel('Vertical Velocity, x_2 [ft/s]');
title('Velocity vs Time (ODE45)');

% Part b - Sensitivity to ballistic coefficient
deltas = [0, -5, -100, -1000];    % delta x3 [lb/ft^2]
colors = lines(numel(deltas));
labels   = cell(1,numel(deltas));
tf_b   = zeros(1,numel(deltas));
vi_b   = zeros(1,numel(deltas));

% Altitude plot w/ different betas
fig2 = figure('Position',[140 140 1200 800],'Color','w'); hold on; grid minor;
for i = 1:numel(deltas)
    beta0 = 2000 + deltas(i);
    x0 = [60000; 0; beta0];
    [t_i, x_i] = ode45(@(t,x) F_nonlinear(t,x,par), tspan, x0, opts_event);

    plot(t_i, x_i(:,1), 'LineWidth', 1.8, 'Color', colors(i,:));
    labels{i} = sprintf('\\beta_0 = %d lb/ft^2', beta0);
    tf_b(i) = t_i(end);
    vi_b(i) = x_i(end,2);
end
xlabel('Time, t [s]'); ylabel('Altitude, x_1 [ft]');
title('Altitude Trajectories - Ballistic coefficient deviations');
legend(labels, 'Location', 'best');

fprintf('\nPart b impact results (nonlinear truth):\n');
for i = 1:numel(deltas)
    fprintf('  beta0 = %5d lb/ft^2:  tf ≈ %7.1f s,  v_impact ≈ %8.1f ft/s\n', ...
            2000 + deltas(i), tf_b(i), vi_b(i));
end

% Part c Analytic Jacobian F(x) 

syms x_1 x_2 x_3 rho0 kp g real
d = 0.5*rho0*exp(-x_1/kp)*x_2^2/x_3;
A = [x_2; d - g; 0];
STM = jacobian(A, [x_1 x_2 x_3]);

fprintf('\nPart c Analytical Jacobian:\n');
pretty(STM)

% Part d Linearized model
Phi0 = eye(3);    % initial STM
z0   = [x0_nom; Phi0(:)];  % augmented initial state w/ Phi verticalized
tspan_nom = [0 tf_nom];  % integrate to nominal impact time

% Integrate nominal x*(t) and STM simultaneously
[t_phi, z_phi] = ode45(@(t,z) dyn_with_stm(t,z,par), tspan_nom, z0);
x_star   = z_phi(:,1:3);                     % nominal trajectory x*(t)
Phi_all  = z_phi(:,4:end);                   % STM entries (rows hold 9 elements)
Phi_tf   = reshape(Phi_all(end,:).', 3,3);   

fprintf('\nPart (d): Linear (STM) prediction at t_f = %.1f s vs nonlinear truth:\n', tf_nom);
for i = 1:numel(deltas)
    delta0   = [0; 0; deltas(i)];
    x_lin_tf = x_nom_tf + Phi_tf * delta0;                 % linearized prediction at t_f
    [~, x_chk] = ode45(@(t,x) F_nonlinear(t,x,par), [0 tf_nom], [60000; 0; 2000+deltas(i)]);
    x_true_tf = x_chk(end,:).';
    err = x_true_tf - x_lin_tf;                            % error at t_f
    fprintf('  beta0 = %5d:  |err| = [%.2e, %.2e, %.2e]^T   (x1,x2,x3)\n', ...
            2000 + deltas(i), err(1), err(2), err(3));
end

fig3 = figure('Position',[200 200 1200 800],'Color','w');
tiledlayout(fig3,2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile; hold(ax1,'on'); grid(ax1,'minor');
title(ax1,'Part (d): Altitude Error (Linear − Nonlinear Truth)');
xlabel(ax1,'Time, t [s]'); ylabel(ax1,'Altitude error [ft]');

ax2 = nexttile; hold(ax2,'on'); grid(ax2,'minor');
title(ax2,'Part (d): Velocity Error (Linear − Nonlinear Truth)');
xlabel(ax2,'Time, t [s]'); ylabel(ax2,'Velocity error [ft/s]');

fprintf('\nPart (d) error summaries over 0..t_f (Linear − Truth):\n');
for i = 1:numel(deltas)
    delta0 = [0;0;deltas(i)];

    % Build linearized state over t_phi
    x_lin = zeros(numel(t_phi),3);
    for k = 1:numel(t_phi)
        Phi_k = reshape(Phi_all(k,:).', 3,3);
        x_lin(k,:) = (x_star(k,:).' + Phi_k * delta0).';
    end

    % Nonlinear truth sampled at the SAME times (no event stop)
    [~, x_truth] = ode45(@(t,x) F_nonlinear(t,x,par), t_phi(:).', [60000;0;2000+deltas(i)]);

    % Errors (linear − truth)
    e_alt = x_lin(:,1) - x_truth(:,1);
    e_vel = x_lin(:,2) - x_truth(:,2);

    plot(ax1, t_phi, e_alt, 'LineWidth', 1.5, 'Color', colors(i,:));
    plot(ax2, t_phi, e_vel, 'LineWidth', 1.5, 'Color', colors(i,:));

    % Final and RMS errors
    e_alt_final = e_alt(end);   e_vel_final = e_vel(end);
    e_alt_rms   = sqrt(mean(e_alt.^2));
    e_vel_rms   = sqrt(mean(e_vel.^2));
    fprintf('  Δβ=%6d:  final [alt, vel] = [%9.3f ft, %9.3f ft/s],  RMS = [%9.3f ft, %9.3f ft/s]\n', ...
            deltas(i), e_alt_final, e_vel_final, e_alt_rms, e_vel_rms);
end

legend(ax1, arrayfun(@(b)sprintf('\\Delta\\beta=%d',b),deltas,'UniformOutput',false), 'Location','best');
legend(ax2, arrayfun(@(b)sprintf('\\Delta\\beta=%d',b),deltas,'UniformOutput',false), 'Location','best');


function dx = F_nonlinear(~, x, par)

    rho     = par.rho0 * exp(-x(1)/par.kp);     
    drag_up = 0.5 * rho * (x(2)^2) / x(3);     
    dx = [ x(2);
           drag_up - par.g;
           0 ];
end

function [value, isterminal, direction] = impactEvent(~, x)
    value      = x(1);   % altitude
    isterminal = 1;      % stop integration
    direction  = -1;     % only when crossing downward
end

function [A] = F_jacobian(x, par)
    x1 = x(1); x2 = x(2); x3 = x(3);
    
    rho = par.rho0 * exp(-x1/par.kp);

    S   = 0.5 * rho * (x(2)^2) / x3;

    A = [ 0,            1,           0;
         -S/par.kp,  rho*x2/x3,   -S/x3;
          0,            0,           0 ];
end

function dz = dyn_with_stm(~, z, par)
    x   = z(1:3);
    Phi = reshape(z(4:end), 3,3);

    fx  = F_nonlinear([], x, par);
    F   = F_jacobian(x, par);

    Phi_dot = F * Phi;
    dz = [fx; Phi_dot(:)];
end