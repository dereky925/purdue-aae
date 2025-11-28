% =========================================================================
% 
% Filename:       HW4.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE590 - Satellite Constellations and Formation
% Professor:      Dr. David Arnas
% Assignment:     HW 4
% Semester:       Spring 2025
% 
% Description: Homework 4
%
% =========================================================================

%% Problem 1

clear; clc; close all;

deg2rad       = pi/180;
rad2deg       = 180/pi;

% Orbit altitude [km]
h = 650;
a = EARTH_RADIUS + h;  % [km], circular orbit => semi-major axis

% Sun-sync RAAN regression rate
desiredRegression_degday = -0.9856;   

% RAAN Error ±15 min => 30 min * 15°/hr Earth rate
dRAAN_deg = 0.5*EARTH_RATE;
dRAAN_rad = dRAAN_deg * deg2rad;   % convert to radians


% Satellite's inclination drift rate: 0.05 deg/year
% Convert to [rad/s]
didt_deg_per_year  = 0.05; 
didt_deg_per_sec   = didt_deg_per_year / (365*24*3600);
didt_deg_per_day   = didt_deg_per_year / 365;  
didt_rad_per_sec   = didt_deg_per_sec * deg2rad;  



n_rad_s   = sqrt(MU / a^3);            % [rad/s]   mean motion
n_deg_day = n_rad_s * rad2deg * 86400; % [deg/day] convert rad/s -> deg/day

K = (3/2)*J2*(EARTH_RADIUS^2 / a^2)*n_deg_day;  % [deg/day]

i0_rad = acos(abs(desiredRegression_degday)/K);
i0_deg = i0_rad*rad2deg;

ecc = 0;
omega_sunsync = 1.991063853E-7; % rad/s

% Compute inclination necessary for sun sync orbit
incl = acosd(-2*a^(7/2)*omega_sunsync*(1-ecc^2)^2/(3*EARTH_RADIUS^2*J2*sqrt(MU)));

% For a sun-sync orbit, inclination is typically > 90°
if desiredRegression_degday < 0
    i0_deg = 180 - i0_deg;
    i0_rad = i0_deg * deg2rad;
end

disp('-------------------------------------------');
disp('Nominal sun-sync inclination:');
disp(['  i0 = ', num2str(i0_deg, '%.4f'), ' deg']);


raan_rate_rad_s = - (3/2)*J2*(EARTH_RADIUS^2 / a^2)*n_rad_s*cos(i0_rad);


Delta_i_rad = -2 * sqrt(   (-2 * dRAAN_rad * didt_rad_per_sec) / ...
                           (raan_rate_rad_s * tan(i0_rad))   );

Delta_i_deg = Delta_i_rad * rad2deg;

disp('-------------------------------------------');
disp('Max allowable inclination variation:');
disp(['  Δi_max = ', num2str(Delta_i_deg, '%.4f'), ' deg']);


% drift rate in rad/day
didt_rad_per_day = -didt_rad_per_sec * 86400; 
time_to_drift_days = Delta_i_rad / didt_rad_per_day;
time_to_drift_years = time_to_drift_days / 365;

disp('-------------------------------------------');
disp('Time between maneuvers:');
disp(['  T_maneuver = ', num2str(time_to_drift_days,'%.2f'), ' days (',...
      num2str(time_to_drift_years,'%.2f'), ' years)']);


v_circ = sqrt(MU / a);  % [km/s]
Delta_v = -2 * v_circ * sin(0.5*Delta_i_rad);  % [km/s]

disp('-------------------------------------------');
disp('Impulse for station-keeping:');
disp(['  Δv = ', num2str(Delta_v, '%.3f'), ' km/s']);


i_mean_deg = i0_deg + 0.5*(Delta_i_deg);

disp('-------------------------------------------');
disp('(d) Mean inclination after station-keeping:');
disp(['  i_mean = ', num2str(i_mean_deg, '%.4f'), ' deg']);


Omega0 = 1.5 * J2 * (EARTH_RADIUS^2 / a^2) * n_deg_day;  % [deg/day]

raan_rate_rad_s = - (3/2)*J2*(EARTH_RADIUS^2 / a^2)*n_rad_s*cos(i0_rad);
raan_rate_rad_day = raan_rate_rad_s * 86400;


t_days  = 0:time_to_drift_days;

didt_rad_per_day = didt_deg_per_day*deg2rad;
i_evol_deg = Delta_i_deg/2 + didt_deg_per_day .* t_days;
i_evol_rad = i_evol_deg*deg2rad;

dOmega_rad = dRAAN_rad/2 + -raan_rate_rad_day * tan(i0_rad) .* (Delta_i_rad/2.* t_days +  0.5*didt_rad_per_day.*(t_days.^2) );

dOmega_deg = dOmega_rad .* rad2deg;


figure('Position',[1000 100 1400 1000],'Color','white'); 

set(gca,'FontSize',20);
subplot(2,1,1);
plot(t_days, i_evol_deg, 'LineWidth',2, 'Color','b');
grid minor;
xlabel('Time [days]', 'FontSize',20);
ylabel('i(t) [deg]', 'FontSize',20);
title('Variation of Inclination from Nominal Sun-Sync i_0', 'FontSize',20);
axis('fill')
ax = gca;
ax.FontSize = 20;

subplot(2,1,2);
plot(t_days, dOmega_deg, 'LineWidth',2, 'Color','r');
grid minor;
xlabel('Time [days]', 'FontSize',20);
ylabel('\Omega(t) [deg]', 'FontSize',20);
title('Variation of RAAN from Nominal Sun-Sync \Omega', 'FontSize',20);
axis('tight')
ax = gca;
ax.FontSize = 20;

%% Problem 2

clear;clc;close all

h_min = 650;
h_max = 900;
dh    = 10; % altitude steps
aVec  = (EARTH_RADIUS + h_min) : dh : (EARTH_RADIUS + h_max);

i_min = 0;
i_max = 90;
di    = 1;
iVec  = i_min : di : i_max;


for iInd = 1:length(iVec)
    inc_deg = iVec(iInd);
    inc     = deg2rad(inc_deg);  % inclination in radians

    for aInd = 1:length(aVec)
        a = aVec(aInd);  % semimajor axis [km]

        % Mean motion [rad/s] (unperturbed)
        n = sqrt(MU / a^3);

        % J2-induced precession rates [rad/s]
        % circular orbit so p => a, e is 0. p = a(1-e^2)
        dotOmega = -3/2 * J2 * (EARTH_RADIUS/a)^2 * n * cos(inc);  % RAAN rate

        dotomega =  3/4 * J2 * (EARTH_RADIUS/a)^2 * n * (4-5*sin(inc)^2); % Apsidal rotation rate

        % Nodal period [s] - time between successive passes at ascending node
        P_n = 2*pi / (n - dotOmega);

        % Anomalistic period [s] - time between successive perigee passages
        P_a = 2*pi / (n - dotomega);

        % Difference in minutes
        diff_min = (P_a - P_n) / 60;

        Pa_minus_Pn(iInd, aInd) = diff_min;
    end
end


[Agrid, Igrid] = meshgrid(aVec, iVec);
figure('Color','white', 'Position',[100,100,1400,1000]);

numLevels = 15;
[cs, hCont] = contourf(Agrid, Igrid, Pa_minus_Pn, numLevels, ...
                       'LineColor','k');  % black contour lines
hold on;

cb = colorbar;
cb.Label.String = 'P_a - P_n [minutes]';
cb.FontSize = 20;

colormap('jet');  % color map
% set(cb, 'YDir', 'reverse');

xlabel('Semimajor axis [km]', 'FontSize', 20);
ylabel('Inclination [deg]',   'FontSize', 20);
title('\DeltaP = P_a - P_n [minutes]', 'FontSize', 20);

set(gca, 'FontSize', 20, 'XGrid','on', 'YGrid','on');


%% Problem 3

clear; clc; close all

format long

w_e     = 7.2921159e-5;   % [rad/s]    Earth sidereal rotation rate
i_deg   = 98.5705;        % [deg]      Inclination
i_rad   = i_deg*(pi/180); % [rad]

Np = 385;  % number of orbits for ground track repeat
Nd = 27;   % number of Earth rotations for ground track repeat

f = @(a) w_e - ( -1.5*J2*sqrt(MU/a^3)*(EARTH_RADIUS/a)^2*cos(i_rad) ) ...
         - (Nd/Np)*sqrt(MU/a^3);

a_lower = 6500;  % [km] ~100 km alt
a_upper = 8000;  % [km] ~1500 km alt

% Upper and lower bounds
f_low = f(a_lower);
f_high = f(a_upper);

% Bisection method parameters
tol = 1e-12;      % Tolerance for 'a' (in km)
maxIter = 200;    % Maximum number of iterations

aL = a_lower;
aU = a_upper;

for k = 1:maxIter
    aMid = 0.5*(aL + aU);
    fMid = f(aMid);

    % Check the sign to decide which half to keep
    if f(aL)*fMid < 0
        aU = aMid;
    else
        aL = aMid;
    end
    
    % Check if bounds are sufficiently small
    if abs(aU - aL) < tol
        break;
    end
end

% Compute another midpoint for final answer
a_sol = 0.5*(aL + aU);

disp('==================================================================')
disp('Solution for repeat-ground-track orbit (Bisection):')
fprintf('  Iterations used   = %d\n', k);
fprintf('  Semimajor axis a  = %.9f km\n', a_sol)
% fprintf('  Altitude h        = %.9f km', a_sol - EARTH_RADIUS)

fprintf('\n')
n_sol       = sqrt(MU/a_sol^3);
dotOmega_sol= -1.5*J2*n_sol*(EARTH_RADIUS/a_sol)^2*cos(i_rad);
lhs         = (w_e - dotOmega_sol)/n_sol;    % Should match Nd/Np
rhs         = Nd/Np;
disp('==================================================================')


%% Problem 4

clear;clc;close all
format long

% Given
i_deg = 56;    % inclination (degrees)
e = 0.02;      % eccentricity
dr = 10;       % ±5 km altitude tolerance, delta radial distance

% Ground track repeat data
N_orbits       = 17;  % The satellite completes 17 orbits
N_earthRot     = 10;  % while Earth rotates 10 times (sidereal)
T_earthSidereal= 86164; % [s] approximate sidereal day of Earth

% Argument of perigee drift due to Moon
wdot_moon = 2.63e-8; % [rad/s]

T_orb = (N_earthRot / N_orbits) * T_earthSidereal;  % [s]

a = ( (T_orb^2 * MU) / (4*pi^2) )^(1/3);   % [km]

% Maximum variation in argp to keep altitude within ±5 km
dw_max = sqrt( 2*dr/(a*e) * (1-e)/(1+e) );

% Compute time between maneuvers
t_maneuver = dw_max / wdot_moon;  % [s]

% Impulse (delta v) needed to rotate the ellipse back by ~dw_max
T = sqrt(a^3/MU);
n = 1/T;

% Compute required delta v (simplified equation)
dV = n*a*e*dw_max/2; 
dV_m_s = dV * 1000; % [m/s]

disp('==================================================================')
fprintf('(a) Compute semi major axis:\n');

fprintf('    - Semimajor axis (a)          = %.2f  km\n\n', a);

fprintf('(b) Max |delta_omega| to keep ±5 km altitude:\n');
fprintf('    - delta_omega_max = %.6f rad  (%.4f deg)\n\n', ...
         dw_max, dw_max*180/pi);

fprintf('(c) Time between maneuvers & impulse:\n');
fprintf('    - Time between maneuvers   = %.2f  s   (%.2f days)\n', ...
         t_maneuver, t_maneuver/86400);
fprintf('    - Required delta-v (1 burn)= %.4f km/s  (%.2f m/s)\n', ...
         dV, dV_m_s);
disp('==================================================================')




