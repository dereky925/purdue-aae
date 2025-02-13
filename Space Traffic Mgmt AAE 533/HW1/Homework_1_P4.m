clear;clc;close all

% LLA of Armstrong Hall
lat = 40.431;
lon = -86.915;
alt = 0;

% Put LLA into vector
lla = [lat,lon,alt];

% Compute ECI XYZ from LLA and UTC time, Aug 27, 8 PM 2024
pos_eci = lla2eci(lla, [2024 8 27 20 00 00])

% Compute geocentric latitude
geoc_lat = atan2d(pos_eci(2),pos_eci(1));

% Add pi to geocentric latitude to get geocentric RA, eqn 23 in notes
geoc_RA_deg = (pi + geoc_lat) * 180/pi;


mu = 3.9860044e14; %m^3/s^2

% GEO Parameters =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

a = 42378137; % Semi-major axis [m] GEO
e = 0.01; % Eccentricity
i = 0.01; % Inclination [deg]
RAAN = 0; % Right Ascension of Ascending Node [deg]
w = 0; % Argument of Perigee [deg]
v = 0; % True Anomaly [deg]

num_orbits = 2;
P_coast = num_orbits * (2*pi)/(sqrt(mu/(a^3))); %Coast time
t = 1:1:P_coast;

[r_ijk,v_ijk] = keplerian2ijk(a,e,i,RAAN,w,v);

options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
[Tout, Z] = ode45(@two_body_ode,t,[r_ijk v_ijk],options);


% Orbit Plot =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
imData = imread('2_no_clouds_4k.jpg');
[xS,yS,zS] = sphere(50);
earth_radius = 6378.1370*1000;  % meters
xSE = earth_radius*xS;
ySE = earth_radius*yS;
zSE = earth_radius*zS;
surface(xSE,ySE,zSE,'FaceColor','texturemap','cdata',flipud(imData),'edgecolor','none');
axis equal
a = view(3);
grid on
xlabel('Inertial x (m)')
ylabel('Inertial y (m)')
zlabel('Inertial z (m)')
ch = get(gca,'children');
hold on 
% Plot observer location in ECI
scatter3(pos_eci(1), pos_eci(2), pos_eci(3), 500,'r.');

% Plot GEO satellite orbit
plot3(Z(:,1), Z(:,2), Z(:,3), 'b','Linewidth', 4);
plot3([0 pos_eci(1)*6], [0 pos_eci(2)*6], [0 pos_eci(3)*6], 'k')
plot3([0 -pos_eci(1)*6], [0 -pos_eci(2)*6], [0 -pos_eci(3)*6], 'k')
plot3([-pos_eci(1)*9 pos_eci(1)*9], [0 0], [0 0], 'm')
legend('','Observer','GEO','FontSize',15)

% Create the tangent plane at the observer's ECI location
normal_vector = pos_eci / norm(pos_eci); % Normalized normal vector to the plane

plane_size = earth_radius * 8; % Set size for plane
[u, v] = meshgrid(-plane_size:plane_size/10:plane_size);

% Create arbitrary vectors perpendicular to the normal vector
tangent1 = null(normal_vector(:).'); 

plane_points = pos_eci' + tangent1(:,1) * u(:)' + tangent1(:,2) * v(:)';

plane_points = reshape(plane_points', [size(u,1), size(u,2), 3]);

% Plot the plane
surf(plane_points(:,:,1), plane_points(:,:,2), plane_points(:,:,3), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', 'cyan');

hold off

%% ECI to RA

point1 = [-15140200, 39732800, 693538];
point2 = [-2148310, -42334500, -738951];

% Compute geocentric latitude
geoc_lat_1 = atan2(point1(2),point1(1))
% Add pi to geocentric latitude to get geocentric RA, eqn 23 in notes
geoc_RA_deg_1 = (pi + geoc_lat_1) * 180/pi

% Compute geocentric latitude
geoc_lat_2 = atan2(point2(2),point2(1))
% Add pi to geocentric latitude to get geocentric RA, eqn 23 in notes
geoc_RA_deg_2 = (pi + geoc_lat_2) * 180/pi

