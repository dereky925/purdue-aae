% =========================================================================
% 
% Filename:       HW7_P1.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE533 - Space Traffic Management
% Professor:      Dr. Carolin Frueh
% Contact:        cfrueh@purdue.edu
% Assignment:     HW 7
% Semester:       Fall 2024
% 
% Description:
% J2 perturbations
% Jacobian of 3rd body and J2 effects
%
% =========================================================================

clear;clc;close all

mu_Earth = 3.9860044e14; %m^3/s^2

utc_datetime = datetime(2024, 10, 16, 12, 0, 0, 'TimeZone', 'UTC');
JD = juliandate(utc_datetime);

% ISS pos vel vector 10 17 2024 12:00:00.00
r_ijk = [656.0758  6596.8028 1473.0872]' * 1000; % [m]
v_ijk = [-4.9631 -0.8127 5.7866]' * 1000; % [m/s]

[a,ecc,incl,RAAN,argp,nu,truelon,arglat,lonper] = ...
                                            ijk2keplerian(r_ijk, v_ijk)

[r_ijk, v_ijk] = keplerian2ijk(a,ecc,incl,RAAN,345,nu);

num_orbits = 100;
P_coast = num_orbits * (2*pi)/(sqrt(mu_Earth/(a^3))); %Coast time
t = 1:1:P_coast;

% COMPUTE 4 BODY MOTION AND PERTURBATIONS
options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
[Tout, Z] = ode45(@(t,x) four_body_ode(t,x, utc_datetime),...
                t, [r_ijk v_ijk], options);

% Compute change in RAAN per day
for i = 1:86400:length(Z)
    
    [~,~,~,curr_RAAN,~,~,~,~,~] = ijk2keplerian(Z(i,1:3), Z(i,4:6));


    if i ~= 1
        dRAAN((i-1)/86400) = curr_RAAN - prev_RAAN;
    end

    prev_RAAN = curr_RAAN;
    
end
disp('RAAN rate:')
disp(mean(dRAAN))
% Z(end,1:3)

raan_plot = zeros(1,length(Z));
argp_plot = zeros(1,length(Z)); 
for i = 1:1:length(Z)
    [~,~,~,raan_plot(i),argp_plot(i),~,~,~,~] = ijk2keplerian(Z(i,1:3), Z(i,4:6));
end

figure('Color','white','Position',[0 0 1400 1000])
hold on;grid minor
plot(Tout/86400,argp_plot)
title('\omega wrapping')
xlabel('Time [days]')
ylabel('\omega [°]')

ax=gca;
ax.FontSize = 20;
% plot(raan_plot)
%%

% Orbit Plot  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
plotA = figure();
hold on
imData = imread('2_no_clouds_4k.jpg');  % Load the texture image

[xS, yS, zS] = sphere(50);  % Generate a sphere for the globe
earth_radius = 6378.1370 * 1000;  % Earth radius in meters

xSE = earth_radius * xS;
ySE = earth_radius * yS;
zSE = earth_radius * zS;

% Create the rotation matrices
lon = 1;
Ry = R3(-lon+.6);

Rx = R1(0);

% Combine the rotations (Ry first for longitude, then Rx for latitude)
rotationMatrix = Rx * Ry;

% Apply the rotation to the coordinates
rotated_coords = rotationMatrix * [xSE(:)'; ySE(:)'; zSE(:)'];

% Reshape the rotated coordinates back to match the surface grid dimensions
xSE_rot = reshape(rotated_coords(1, :), size(xSE));
ySE_rot = reshape(rotated_coords(2, :), size(ySE));
zSE_rot = reshape(rotated_coords(3, :), size(zSE));

% Plot the Earth with texture mapping, using the rotated coordinates
surface(xSE_rot, ySE_rot, zSE_rot, 'FaceColor', 'texturemap', 'CData', ...
    flipud(imData), 'EdgeColor', 'none');

% Adjust axes
axis equal;
view(3);  % 3D view
grid on;
xlabel('Inertial x (m)');
ylabel('Inertial y (m)');
zlabel('Inertial z (m)');

plot3(Z(:,1), Z(:,2), Z(:,3), 'r', 'LineWidth', 4);

plotA.Position = [1600 100 1000 1000];


%%
clear;clc;close all

syms x y z u v w mu_Earth mu_Sun mu_Moon ...
    x_Sun y_Sun z_Sun x_Moon y_Moon z_Moon J2 RE real

r = sqrt(x^2 + y^2 + z^2);
r_sun = sqrt(x_Sun^2 + y_Sun^2 + z_Sun^2);
r_moon = sqrt(x_Moon^2 + y_Moon^2 + z_Moon^2);

s = x/r;
t = y/r;
u = z/r;

u_vec = [s; t; u];

u_J2 = 3/2 * [5*s*u^2 - s; 5*t*u^2 - t; 5*u^3 - 3*u];

ax = - mu_Sun*((x - x_Sun)/norm(x-x_Sun)^3 + x_Sun/r_sun^3)...
     - mu_Moon*((x - x_Moon)/norm(x-x_Moon)^3 + x_Moon/r_moon^3)...
     - mu_Earth/r^2*u_vec(1) + mu_Earth*J2/r^2*(RE/r)*u_J2(1);

ay = - mu_Sun*((y - y_Sun)/norm(y-y_Sun)^3 + y_Sun/r_sun^3)...
     - mu_Moon*((y - y_Moon)/norm(y-y_Moon)^3 + y_Moon/r_moon^3)...
     - mu_Earth/r^2*u_vec(2) + mu_Earth*J2/r^2*(RE/r)*u_J2(2);

az = - mu_Sun*((z - z_Sun)/norm(z-z_Sun)^3 + z_Sun/r_sun^3)...
     - mu_Moon*((z - z_Moon)/norm(z-z_Moon)^3 + z_Moon/r_moon^3)...
     - mu_Earth/r^2*u_vec(3) + mu_Earth*J2/r^2*(RE/r)*u_J2(3);


% X

dax_dx = diff(ax, x);
dax_dy = diff(ax, y);
dax_dz = diff(ax, z);

pretty(dax_dx)
pretty(dax_dy)
pretty(dax_dz)

% Y

day_dx = diff(ay, x);
day_dy = diff(ay, y);
day_dz = diff(ay, z);

pretty(day_dx)
pretty(day_dy)
pretty(day_dz)

% Z

daz_dx = diff(az, x);
daz_dy = diff(az, y);
daz_dz = diff(az, z);

pretty(daz_dx)
pretty(daz_dy)
pretty(daz_dz)











