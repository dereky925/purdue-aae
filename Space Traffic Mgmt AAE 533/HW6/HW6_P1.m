clear;clc;close all

mu = 3.9860044e14; %m^3/s^2

utc_datetime = datetime(2024, 10, 16, 12, 0, 0, 'TimeZone', 'UTC');
JD = juliandate(utc_datetime);

% ISS pos vel vector 10 17 2024 12:00:00.00
r_ijk = [656.0758  6596.8028 1473.0872]' * 1000; % [m]
v_ijk = [-4.9631 -0.8127 5.7866]' * 1000; % [m/s]

[a,ecc,incl,RAAN,argp,nu,truelon,arglat,lonper] = ...
                                            ijk2keplerian(r_ijk, v_ijk);

num_orbits = 100;
P_coast = num_orbits * (2*pi)/(sqrt(mu/(a^3))); %Coast time
t = 1:1:P_coast;

% COMPUTE 4 BODY MOTION AND PERTURBATIONS
options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
[Tout, Z] = ode45(@(t,x) four_body_ode(t,x, JD), t, [r_ijk v_ijk], options);


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


% 1.b. Compute life of ISS

ISS_C_D = 2.07;
% https://www.esa.int/Education/Space_In_Bytes/ATV_a_very_special_delivery_-_Lesson_notes

air_density_ISS = 3.8e-12 ; % kg/m³ air density at 400 km 
% https://www.esa.int/Education/Space_In_Bytes/ATV_a_very_special_delivery_-_Lesson_notes

ISS_DragArea = 1500; % m^2 
% https://www.esa.int/Education/Space_In_Bytes/ATV_a_very_special_delivery_-_Lesson_notes

ISS_Mass = 450000; % kg
% https://en.wikipedia.org/wiki/International_Space_Station

C_d = 0.4; % ISS is mostly made of alluminum
C = (1/4 + 1/9*C_d);

% Single step approximation
da = -2*pi*ISS_C_D*ISS_DragArea/ISS_Mass*air_density_ISS*a^2;

H = 58.2e3; % [m] Scale height from notes Table, page 80

L = -H/da;
time = L * 1.5/24; % [days]

% Analytical computation
height = 400 * 1000; % [m]
Earth_Radius = 6.371e6; % [m]
a_ISS = height+Earth_Radius; % [m]

% Air density table
ref_heights   = [0, 100, 150, 200, 250, 300, 350, 400]*1000; % m
air_densities = [1.2, 4.79e-7, 1.81e-9, 2.53e-10,...
                 6.24e-11, 1.95e-11, 6.98e-12, 2.72e-12]; % kg/m^3


%%
for i = 1:100000

    % Find the closest reference height to the input height
    [~, idx] = min(abs(ref_heights - height));
    closest_ref_height = ref_heights(idx) * 1000; % Convert back to meters
    
    % Get the corresponding air density
    air_density = air_densities(idx);
  
    % Compute air density at ISS altitude
    air_density_ISS_h = (air_density)*exp(-height/closest_ref_height);

    % Compute change in semi-major axis in one revolution - circular orbit
    da = -2*pi * ISS_C_D*ISS_DragArea/ISS_Mass * air_density_ISS_h * a_ISS;
    height = a_ISS + da - Earth_Radius % [m]

end
