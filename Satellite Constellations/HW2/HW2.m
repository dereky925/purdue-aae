% =========================================================================
% 
% Filename:       HW2.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE590 - Satellite Constellations and Formation
% Professor:      Dr. David Arnas
% Assignment:     HW 2
% Semester:       Spring 2025
% 
% Description: Homework 2
%
% =========================================================================

%% Problem 1
clear;clc;close all

mu = 398600.4418; % [km^3/s^2]

orbit_plot = OrbitPlotSetup();
orbit_plot.Position = [0 0 1000 1000];

altitude = 800; % [km] Choose altitude
a = altitude + EARTH_RADIUS; % Compute a
raan = 0;   % Choose raan
ecc = 0.01;  % Choose eccentricity
nu = 0;     % Choose nu
argp = 0;   % Choose argp

omega_sunsync = 1.991063853E-7; % rad/s

% Compute inclination necessary for sun sync orbit
incl = acosd(-2*a^(7/2)*omega_sunsync*(1-ecc^2)^2/(3*EARTH_RADIUS^2*J2*sqrt(MU)));

[r_vec,v_vec] = keplerian2eci(a, ecc, incl, raan, argp, nu);

num_orbits = 1;
T = num_orbits*2*pi*sqrt(a^3/mu); % Orbit period
t = 1:T;

utc_datetime = datetime(2024, 10, 19, 17, 54, 19, 'TimeZone', 'UTC');

options = odeset('RelTol',1E-12,'AbsTol',1E-12);
[t_out,eci_out] = ode45(@(t,x) two_body_J2_ode(t,x,utc_datetime),...
                t, [r_vec v_vec], options);

altitude = vecnorm(eci_out(:,1:3),2,2) - EARTH_RADIUS;

for i = 1:length(t)
    [a(i), ecc(i), incl(i), raan(i), argp(i), nu(i)] = eci2keplerian(eci_out(i,1:3),eci_out(i,4:6));
end

plot3(eci_out(:,1),eci_out(:,2),eci_out(:,3),'r','LineWidth',4)

i_fig = figure('Color','white');
hold on
grid minor
ax = gca;
ax.FontSize = 20;
i_fig.Position = [1100 0 1000 800];
plot(altitude,incl,'LineWidth',4)
title('Inclination vs. Altitude')
xlabel('Altitude [km]')
ylabel('Inclination [°]')
xlim([500 1200])



%% Problem 2
% Create J2 propagator

%% Problem 3
% 

%% Problem 4
clear;clc;close all

a = 6994.471731884392; % [km]
ex = 0;
ey = -4.976E-4;
incl = 49.981; % [°]
raan = 0;   % [°]
theta = 0;  % [°] Argument of latitude

T = 1.617608582630545; % [hours]
t = 1:T*3600;

[r_vec, v_vec] = keplerian2eci(a,ex,ey,incl,raan,theta,'alt');
% [a, ecc, incl, raan, ~, ~, ecc_x, ecc_y, arglat] = eci2keplerian(r_vec, rdot_vec);

% [r_vec, v_vec] = elements2cartesian(a, ex, ey, incl, raan, theta);

orbit_plot = OrbitPlotSetup();
orbit_plot.Position = [0 0 1000 1000];

utc_datetime = datetime(2024, 10, 19, 17, 54, 19, 'TimeZone', 'UTC');

options = odeset('RelTol',1E-12,'AbsTol',1E-12);
[t_out,eci_out] = ode45(@(t,x) two_body_J2_ode(t,x,utc_datetime),...
                t, [r_vec v_vec], options);

for i = 1:length(t)
    [a(i), ecc(i), incl(i), raan(i),...
        ~, ~, ecc_x(i), ecc_y(i), arglat(i)] = ...
            eci2keplerian(eci_out(i,1:3),eci_out(i,4:6));

    % [a(i), ecc_x(i), ecc_y(i), incl(i), raan(i), arglat(i)] ...
    %     = cartesian2elements(eci_out(i,1:3), eci_out(i,4:6));

end

% scatter3(eci_out(1,1),eci_out(1,2),eci_out(1,3),1000,'c.')
plot3(eci_out(:,1),eci_out(:,2),eci_out(:,3),'r','LineWidth',4)


figure(2);
set(gcf, 'Color', 'white', 'Position', [0 300 2000 1000]); 
sgtitle('Orbital Elements vs Time', 'FontSize', 20, 'FontWeight', 'bold');


% Subplot 1: Semi-major Axis
subplot(2,3,1);
hold on
h1 = plot(t, a, 'LineWidth', 4);
grid minor;
xlabel('Time (s)');
ylabel('a [km]');
title('Semi-major Axis');

% Subplot 2: Inclination
subplot(2,3,2);
hold on
plot(t, incl, 'LineWidth', 4);
grid minor;
xlabel('Time (s)');
ylabel('Inclination (°)');
title('Inclination');

% Subplot 3: RAAN
subplot(2,3,3);
hold on
plot(t, raan, 'LineWidth', 4);
grid minor;
xlabel('Time (s)');
ylabel('RAAN (°)');
title('RAAN');

% Subplot 4: Eccentricity x-component
subplot(2,3,4);
hold on
plot(t, ecc_x, 'LineWidth', 4);
grid minor;
xlabel('Time (s)');
ylabel('e_x');
title('Eccentricity x-component');

% Subplot 5: Eccentricity y-component
subplot(2,3,5);
hold on
plot(t, ecc_y, 'LineWidth', 4);
grid minor;
xlabel('Time (s)');
ylabel('e_y');
title('Eccentricity y-component');

% Subplot 6: Argument of Latitude
subplot(2,3,6);
hold on
plot(t, arglat, 'LineWidth', 4);
grid minor;
xlabel('Time (s)');
ylabel('ArgLat (°)');
title('Argument of Latitude');

set(findall(gcf, 'Type', 'axes'), 'FontSize', 20);


figure(3);
set(gcf, 'Color', 'white', 'Position', [500 600 2000 1000]);
hold on
grid minor;
e1 = plot(ecc_x, ecc_y, 'LineWidth', 4);
xlabel('e_x');
ylabel('e_y');
title('e_x vs e_y');


% =========================================================================
% =========================================================================

a = 6994.471731884392; % [km]
ex = 0;
ey = -4.976E-4;
incl = 49.981; % [°]
raan = 0;   % [°]
theta = 90;  % [°] Argument of latitude


[r_vec, v_vec] = keplerian2eci(a,ex,ey,incl,raan,theta,'alt');
% [a, ecc, incl, raan, ~, ~, ecc_x, ecc_y, arglat] = eci2keplerian(r_vec, rdot_vec);

% [r_vec, v_vec] = elements2cartesian(a, ex, ey, incl, raan, theta);

utc_datetime = datetime(2024, 10, 19, 17, 54, 19, 'TimeZone', 'UTC');

options = odeset('RelTol',1E-12,'AbsTol',1E-12);
[t_out,eci_out] = ode45(@(t,x) two_body_J2_ode(t,x,utc_datetime),...
                t, [r_vec v_vec], options);

for i = 1:length(t)
    [a(i), ecc(i), incl(i), raan(i),...
        ~, ~, ecc_x(i), ecc_y(i), arglat(i)] = ...
            eci2keplerian(eci_out(i,1:3),eci_out(i,4:6));

    % [a(i), ecc_x(i), ecc_y(i), incl(i), raan(i), arglat(i)] ...
    %     = cartesian2elements(eci_out(i,1:3), eci_out(i,4:6));

end

figure(2)
% Subplot 1: Semi-major Axis
subplot(2,3,1);
h2 = plot(t, a,'r', 'LineWidth', 4);
xlabel('Time (s)');
ylabel('a [km]');
title('Semi-major Axis');

% Subplot 2: Inclination
subplot(2,3,2);
plot(t, incl,'r', 'LineWidth', 4);
xlabel('Time (s)');
ylabel('Inclination (°)');
title('Inclination');

% Subplot 3: RAAN
subplot(2,3,3);
plot(t, raan,'r', 'LineWidth', 4);
xlabel('Time (s)');
ylabel('RAAN (°)');
title('RAAN');

% Subplot 4: Eccentricity x-component
subplot(2,3,4);
plot(t, ecc_x,'r', 'LineWidth', 4);
xlabel('Time (s)');
ylabel('e_x');
title('Eccentricity x-component');

% Subplot 5: Eccentricity y-component
subplot(2,3,5);
plot(t, ecc_y,'r', 'LineWidth', 4);
xlabel('Time (s)');
ylabel('e_y');
title('Eccentricity y-component');

% Subplot 6: Argument of Latitude
subplot(2,3,6);
plot(t, arglat,'r', 'LineWidth', 4);
xlabel('Time (s)');
ylabel('ArgLat (°)');
title('Argument of Latitude');

set(findall(gcf, 'Type', 'axes'), 'FontSize', 20);
legend([h1, h2], {'\theta = 0°', '\theta = 90°'}, 'FontSize', 30, 'Location', 'east');


figure(3)
e2 = plot(ecc_x, ecc_y, 'r', 'LineWidth', 4);
xlabel('e_x');
ylabel('e_y');
title('e_x vs e_y');
set(findall(gcf, 'Type', 'axes'), 'FontSize', 20);
legend([e1, e2], {'\theta = 0°', '\theta = 90°'}, 'FontSize', 30, 'Location', 'best');

%% Problem 5
clear;clc;close all

a = 6994.471731884392; % [km]
ex = 0;
ey = -4.976E-4;
incl = 49.981; % [°]
raan = 0;   % [°]
theta = 0;  % [°] Argument of latitude

T = 100 * 1.617608582630545; % [hours]
t = linspace(1,T*3600,10000);

[r_vec, v_vec] = keplerian2eci(a,ex,ey,incl,raan,theta,'alt');

orbit_plot = OrbitPlotSetup();
orbit_plot.Position = [0 0 1000 1000];

utc_datetime = datetime(2024, 10, 19, 17, 54, 19, 'TimeZone', 'UTC');

options = odeset('RelTol',1E-12,'AbsTol',1E-12);
[~,eci_out] = ode45(@(t,x) two_body_J2_ode(t,x,utc_datetime),...
                t, [r_vec v_vec], options);

for i = 1:length(t)
    [a(i), ecc(i), incl(i), raan(i),...
        argp(i), ~, ecc_x(i), ecc_y(i), arglat(i)] = ...
            eci2keplerian(eci_out(i,1:3),eci_out(i,4:6));

end

plot3(eci_out(:,1),eci_out(:,2),eci_out(:,3),'r','LineWidth',4)

figure(1);
set(gcf, 'Color', 'white', 'Position', [0 300 2000 1000]); 
sgtitle('Orbital Elements vs Time', 'FontSize', 20, 'FontWeight', 'bold');


% Subplot 1: Semi-major Axis
subplot(2,3,1);
hold on
h1 = plot(t, a, 'LineWidth', 4);
grid minor;
xlabel('Time (s)');
ylabel('a [km]');
title('Semi-major Axis');

% Subplot 2: Inclination
subplot(2,3,2);
hold on
plot(t, incl, 'LineWidth', 4);
grid minor;
xlabel('Time (s)');
ylabel('Inclination (°)');
title('Inclination');

% Subplot 3: RAAN
subplot(2,3,3);
hold on
plot(t, raan, 'LineWidth', 4);
grid minor;
xlabel('Time (s)');
ylabel('RAAN (°)');
title('RAAN');

% Subplot 4: Eccentricity x-component
subplot(2,3,4);
hold on
plot(t, ecc_x, 'LineWidth', 4);
grid minor;
xlabel('Time (s)');
ylabel('e_x');
title('Eccentricity x-component');

% Subplot 5: Eccentricity y-component
subplot(2,3,5);
hold on
plot(t, ecc_y, 'LineWidth', 4);
grid minor;
xlabel('Time (s)');
ylabel('e_y');
title('Eccentricity y-component');

% Subplot 6: Argument of Latitude
subplot(2,3,6);
hold on
plot(t, arglat, 'LineWidth', 4);
grid minor;
xlabel('Time (s)');
ylabel('ArgLat (°)');
title('Argument of Latitude');

set(findall(gcf, 'Type', 'axes'), 'FontSize', 20);


figure(2);
set(gcf, 'Color', 'white', 'Position', [500 600 2000 1000]);
hold on
grid minor;
e1 = plot(ecc_x, ecc_y, 'LineWidth', 4);
xlabel('e_x');
ylabel('e_y');
title('e_x vs e_y');

figure(3);
set(gcf, 'Color', 'white', 'Position', [500 1000 2000 1000]);
hold on
grid minor;
d1 = plot(t, argp, 'LineWidth', 4);
xlabel('Time (s)');
ylabel('Argument of Perigee (°)');
title('Argument of Perigee vs. Time');


% =========================================================================
% =========================================================================

a = 6994.471731884392; % [km]
ex = 0;
ey = -4.976E-4;
incl = 49.981; % [°]
raan = 0;   % [°]
theta = 90;  % [°] Argument of latitude

[r_vec, v_vec] = keplerian2eci(a,ex,ey,incl,raan,theta,'alt');
% [a, ecc, incl, raan, ~, ~, ecc_x, ecc_y, arglat] = eci2keplerian(r_vec, rdot_vec);

% [r_vec, v_vec] = elements2cartesian(a, ex, ey, incl, raan, theta);

utc_datetime = datetime(2024, 10, 19, 17, 54, 19, 'TimeZone', 'UTC');

options = odeset('RelTol',1E-12,'AbsTol',1E-12);
[t_out,eci_out] = ode45(@(t,x) two_body_J2_ode(t,x,utc_datetime),...
                t, [r_vec v_vec], options);

for i = 1:length(t)
    [a(i), ecc(i), incl(i), raan(i),...
        argp(i), ~, ecc_x(i), ecc_y(i), arglat(i)] = ...
            eci2keplerian(eci_out(i,1:3),eci_out(i,4:6));

    % [a(i), ecc_x(i), ecc_y(i), incl(i), raan(i), arglat(i)] ...
    %     = cartesian2elements(eci_out(i,1:3), eci_out(i,4:6));

end

figure(1)
% Subplot 1: Semi-major Axis
subplot(2,3,1);
h2 = plot(t, a,'r', 'LineWidth', 4);
xlabel('Time (s)');
ylabel('a [km]');
title('Semi-major Axis');

% Subplot 2: Inclination
subplot(2,3,2);
plot(t, incl,'r', 'LineWidth', 4);
xlabel('Time (s)');
ylabel('Inclination (°)');
title('Inclination');

% Subplot 3: RAAN
subplot(2,3,3);
plot(t, raan,'r', 'LineWidth', 4);
xlabel('Time (s)');
ylabel('RAAN (°)');
title('RAAN');

% Subplot 4: Eccentricity x-component
subplot(2,3,4);
plot(t, ecc_x,'r', 'LineWidth', 4);
xlabel('Time (s)');
ylabel('e_x');
title('Eccentricity x-component');

% Subplot 5: Eccentricity y-component
subplot(2,3,5);
plot(t, ecc_y,'r', 'LineWidth', 4);
xlabel('Time (s)');
ylabel('e_y');
title('Eccentricity y-component');

% Subplot 6: Argument of Latitude
subplot(2,3,6);
plot(t, arglat,'r', 'LineWidth', 4);
xlabel('Time (s)');
ylabel('ArgLat (°)');
title('Argument of Latitude');

set(findall(gcf, 'Type', 'axes'), 'FontSize', 20);
legend([h1, h2], {'\theta = 0°', '\theta = 90°'}, 'FontSize', 30, 'Location', 'east');


figure(2)
e2 = plot(ecc_x, ecc_y, 'r', 'LineWidth', 4);
xlabel('e_x');
ylabel('e_y');
title('e_x vs e_y');
set(findall(gcf, 'Type', 'axes'), 'FontSize', 20);
legend([e1, e2], {'\theta = 0°', '\theta = 90°'}, 'FontSize', 30, 'Location', 'best');


figure(3);
d2 = plot(t, argp, '--','LineWidth', 4);
set(findall(gcf, 'Type', 'axes'), 'FontSize', 20);
legend([d1, d2], {'\theta = 0°', '\theta = 90°'}, 'FontSize', 30, 'Location', 'best');






