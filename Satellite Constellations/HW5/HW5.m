% =========================================================================
% 
% Filename:       HW5.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE590 - Satellite Constellations and Formation
% Professor:      Dr. David Arnas
% Assignment:     HW 5
% Semester:       Spring 2025
% 
% Description: Homework 5
%
% =========================================================================

%% Problem 1 

clear; close all; clc;

alt    = 800.0;             % [km]       Orbit altitude
r0     = EARTH_RADIUS + alt;          % [km]       Circular orbit radius
incl   = 98.5705*(pi/180);    % [rad]      Inclination
Omega  = 0.0;               % [rad]      RAAN
omega0 = 0.0;               % [rad]      Arg. of perigee (for circular, arbitrary)
f0     = 90.0*(pi/180);      % [rad]      True anomaly at separation
dv     = 0.005;             % [km/s]     Delta-V imparted to satellite (5 m/s)
Tcirc  = 2*pi*sqrt(r0^3/MU); % [s]      Period of the circular orbit (upper stage)
tspan  = [0 2*Tcirc];        % We will siMUlate 2 orbits of the upper stage

% (a) Find the new orbital elements for the satellite right after separation
v_circ = sqrt(MU/r0);      % [km/s] speed of circular orbit
v_sat  = v_circ + dv;       % [km/s] new speed
a_sat  = 1 / (2/r0 - (v_sat^2/MU));  % [km] semimajor axis of new orbit
e_sat  = 1 - r0/a_sat;                % separation point is perigee => f'=0
omega_sat = f0;                       % argument of perigee = old lat = 90 deg
incl_sat  = incl;                     % same inclination
Omega_sat = Omega;                    % same RAAN

fprintf('=== (a) Satellite Keplerian Elements after Separation ===\n');
fprintf('a     = %.4f km\n', a_sat);
fprintf('e     = %.6f\n', e_sat);
fprintf('i     = %.4f deg\n', incl_sat*180/pi);
fprintf('Omega = %.4f deg\n', Omega_sat*180/pi);
fprintf('omega = %.4f deg\n', omega_sat*180/pi);
fprintf('f0    = 0 (perigee at separation)\n\n');

% (b) Plot in an inertial frame the relative position for 2 orbits
oe_upper = [r0, 0, incl, Omega, omega0, f0]; 
[r_upper0, v_upper0] = kepler2cart(oe_upper, MU);

oe_sat   = [a_sat, e_sat, incl_sat, Omega_sat, omega_sat, 0];
[r_sat0, v_sat0] = kepler2cart(oe_sat, MU);

X0_upper = [r_upper0; v_upper0];
X0_sat   = [r_sat0;   v_sat0];

opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-13);

[tvec_upper, X_upper] = ode45(@(t,X) twoBodyEOM(t,X,MU), tspan, X0_upper, opts_ode);
[tvec_sat,   X_sat]   = ode45(@(t,X) twoBodyEOM(t,X,MU), tspan, X0_sat,   opts_ode);

t_common = linspace(0, tspan(2), 1000);
X_upperI = interp1(tvec_upper, X_upper, t_common);
X_satI   = interp1(tvec_sat,   X_sat,   t_common);

r_upper_all = X_upperI(:,1:3);
r_sat_all   = X_satI(:,1:3);
r_rel_all   = r_sat_all - r_upper_all;  % [km]

% Figure (b)
figure('Name','(b) Relative Trajectory in ECI Frame',...
       'Color','white','Position',[100,100,1000,1000]);
hold on
scatter3(r_rel_all(1,1), r_rel_all(1,2), r_rel_all(1,3),500,'py','MarkerEdgeColor','k');
plot3(r_rel_all(:,1), r_rel_all(:,2), r_rel_all(:,3),'LineWidth',1.5);

grid on; axis equal;
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('Relative Position (Satellite w.r.t. Upper Stage) in ECI');
set(gca,'FontSize',20);
view(45,45)
axis tight
legend('Orbit Start')

% (c) Plot in the satellite frame (R-hat radial, S-hat along-track)
RSW_rel = zeros(length(t_common),3);
for k = 1:length(t_common)
    ru   = r_upper_all(k,:).';
    vu   = X_upperI(k,4:6).';
    rsat = r_sat_all(k,:).';
    rrel = rsat - ru;
    Rhat = ru / norm(ru);
    What = cross(ru,vu); 
    What = What / norm(What);
    Shat = cross(What, Rhat);
    RSW_rel(k,1) = dot(rrel, Rhat);
    RSW_rel(k,2) = dot(rrel, Shat);
    RSW_rel(k,3) = dot(rrel, What);
end

% Figure (c)
figure('Name','(c) Relative Trajectory in Satellite Frame',...
       'Color','white','Position',[200,100,1000,1000]);
hold on
scatter(RSW_rel(1,2), RSW_rel(1,1),300,'py','MarkerEdgeColor','k','MarkerFaceColor','yellow');
plot(RSW_rel(:,2), RSW_rel(:,1),'LineWidth',1.5);

grid on; axis equal;
xlabel('S-hat [km] (along-track)');
ylabel('R-hat [km] (radial)');
title('Relative Motion in Upper Stage Local Orbital Frame');
set(gca,'FontSize',20);
legend('Start')

% (d) Repeat relative motion using Clohessy–Wiltshire–Hill (CWH)
n = sqrt(MU / r0^3);

% initial conditions in local orbital frame at t=0
r_upper_0 = r_upper_all(1,:).';
v_upper_0 = X_upperI(1,4:6).';
r_sat_0   = r_sat_all(1,:).';
v_sat_0   = X_satI(1,4:6).';

Rhat_0 = r_upper_0 / norm(r_upper_0);
What_0 = cross(r_upper_0,v_upper_0); What_0 = What_0 / norm(What_0);
Shat_0 = cross(What_0, Rhat_0);

r_rel_0 = r_sat_0 - r_upper_0; 
v_rel_0 = v_sat_0 - v_upper_0; 
x0_cwh  = [ dot(r_rel_0,Rhat_0);
            dot(r_rel_0,Shat_0);
            dot(r_rel_0,What_0);
            dot(v_rel_0,Rhat_0);
            dot(v_rel_0,Shat_0);
            dot(v_rel_0,What_0) ];

[tCWH, xCWH_sol] = ode45(@(t,x) cwhEOM(t,x,n), t_common, x0_cwh);

% (d) Compare the CWH solution to the solution obtained in part (c)
% Here we use RSW_rel (exact from part (c)) and compare it with the CWH integration
diff_d = xCWH_sol(:,1:3) - RSW_rel(:,1:3);

% Figure (d): Plot a 2D graph with longitudinal error (S-hat, second component) on x-axis,
% and radial error (R-hat, first component) on y-axis.
figure('Name','CWH vs. Keplerian solution',...
       'Color','white','Position',[300,100,1000,1000]);
plot(diff_d(:,2), diff_d(:,1),'LineWidth',1.5);
grid on; axis equal;
xlabel('Longitudinal Error [km]');
ylabel('Radial Error [km]');
title('CWH vs. Keplerian solution');
set(gca,'FontSize',20);

% (e) Transform the CWH solution back to ECI and provide a 2D plot with time on x-axis
r_cwh_eci = zeros(length(t_common),3);
for k = 1:length(t_common)
    ru   = r_upper_all(k,:).';
    vu   = X_upperI(k,4:6).';
    Rhat_k = ru / norm(ru);
    W_k    = cross(ru, vu); W_k = W_k / norm(W_k);
    Shat_k = cross(W_k, Rhat_k);
    x_cwh = xCWH_sol(k,1);
    y_cwh = xCWH_sol(k,2);
    z_cwh = xCWH_sol(k,3);
    rrel_eci_k = x_cwh*Rhat_k + y_cwh*Shat_k + z_cwh*W_k;
    r_cwh_eci(k,:) = (ru + rrel_eci_k).';
end

diff_eci =  r_cwh_eci - r_sat_all;  % [km]

figure('Name','(e) CWH in ECI vs. Exact',...
       'Color','white','Position',[400,100,1000,1000]);
plot(t_common/(3600), diff_eci(:,1),'LineWidth',1.5); hold on;
plot(t_common/(3600), diff_eci(:,2),'LineWidth',1.5);
plot(t_common/(3600), diff_eci(:,3),'LineWidth',1.5);
grid on;
xlabel('Time [hours]');
ylabel('Position Error [km]');
title('Error in CWH-based Satellite Position in ECI Frame');
legend('X error','Y error','Z error','Location','best');
set(gca,'FontSize',20);

%% Problem 2

clc;

tspan_long = [0, 100*86400];  % 100 days
numPoints_long = 10000;
t_common_long = linspace(tspan_long(1), tspan_long(2), numPoints_long);

% Propagate the upper stage and satellite using the two-body EOM over 100 days
[tvec_upper_long, X_upper_long] = ode45(@(t,X) twoBodyEOM(t,X,MU), tspan_long, X0_upper, opts_ode);
[tvec_sat_long,   X_sat_long]   = ode45(@(t,X) twoBodyEOM(t,X,MU), tspan_long, X0_sat,   opts_ode);

% Interpolate the solutions to a common time vector
X_upper_long_i = interp1(tvec_upper_long, X_upper_long, t_common_long);
X_sat_long_i   = interp1(tvec_sat_long,   X_sat_long,   t_common_long);

% Extract positions and compute the relative distance at each time
r_upper_long = X_upper_long_i(:,1:3);
r_sat_long   = X_sat_long_i(:,1:3);
rel_distance = sqrt(sum((r_sat_long - r_upper_long).^2, 2));  % [km]

% Plot 
figure('Name','Long Term Dynamics','Color','white','Position',[500,100,1000,1000]);
plot(t_common_long/86400, rel_distance, 'LineWidth',2);
xlabel('Time [days]');
ylabel('Relative Distance [km]');
title('Evolution of the Distance between Upper Stage and Satellite over 100 Days');
grid on;
set(gca, 'FontSize',20);

fprintf('Distance between objects: 16.0623 km at 34.7935 days\n');

fprintf('Distance between objects: 18.9787 km at 69.597 days\n');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dX = twoBodyEOM(~, X, MU)
    r  = X(1:3);
    v  = X(4:6);
    r3 = norm(r)^3;
    a  = -MU*r/r3;
    dX = [v; a];
end

function [rECI,vECI] = kepler2cart(oe, MU)

    a     = oe(1);
    e     = oe(2);
    i     = oe(3);
    Omega = oe(4);
    omega = oe(5);
    f     = oe(6);
    
    % Distance from focus
    p     = a*(1 - e^2);
    r_mag = p/(1 + e*cos(f));
    
    % Position in the orbital plane
    rx_op = r_mag*cos(f);
    ry_op = r_mag*sin(f);
    vx_op = -sqrt(MU/p)*sin(f);
    vy_op =  sqrt(MU/p)*(e + cos(f));
    
    % Rotate from perifocal frame to ECI
    ROT = R3(Omega)*R1(i)*R3(omega);
    rECI = ROT*[rx_op; ry_op; 0];
    vECI = ROT*[vx_op; vy_op; 0];
end

function M = R1(angle)
% Basic rotation about x-axis
c = cos(angle); s = sin(angle);
M = [1  0   0
     0  c  -s
     0  s   c];
end


function M = R3(angle)
% Basic rotation about z-axis
c = cos(angle); s = sin(angle);
M = [ c -s  0
      s  c  0
      0  0  1];
end

function dX = cwhEOM(~, x, n)
    X  = x(1);  Y  = x(2);  Z  = x(3);
    Xd = x(4);  Yd = x(5);  Zd = x(6);
    
    % CWH equations
    Xdd =  3*n^2*X + 2*n*Yd;
    Ydd = -2*n*Xd;
    Zdd = -n^2*Z;
    
    dX = [Xd; Yd; Zd; Xdd; Ydd; Zdd];
end