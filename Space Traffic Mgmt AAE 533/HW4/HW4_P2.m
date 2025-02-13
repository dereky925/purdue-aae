clear;clc;close all

mu = 3.9860044e14; %m^3/s^2

% Speed of light [m/s]
c = 3E8;

lat = 33.275825;   %[°]
lon = -111.740806; %[°]
alt = 377;         %[m]

lla = [lat, lon, alt];

offset = 500;

%2024-09-25T12:00:00.000 
utc_datetime = datetime(2024, 9, 25, 0, 00, 000, 'TimeZone', 'UTC');

pos_ISS_ECI = [5117.808565310730 -2877.911425706390 3417.677000888100]' .* 1000; %[m] 
vel_ISS_ECI = [5.04434922610630 3.47434703163623 -4.60296327886324]' .* 1000; %[m/s] 

[a,ecc,incl,RAAN,argp,nu,truelon,arglat,lonper] = ijk2keplerian(pos_ISS_ECI, vel_ISS_ECI);

num_orbits = 1.96;
% num_orbits = 2;
P_coast = num_orbits * (2*pi)/(sqrt(mu/(a^3))); %Coast time
t = 1:1:P_coast;

options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
[Tout, Z] = ode45(@two_body_ode,t,[pos_ISS_ECI vel_ISS_ECI],options);

earth_rot = 7.272E-5; % Earth rotation [rad/s]
RE = 6371*1000; % Earth Radius [m]
max_rho_guess = 4; % [RE]

% OBSERVATION 1 =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
obs1time = utc_datetime + seconds(Tout(end-offset));
JD1 = juliandate(obs1time);
pos_eci1 = lla2eci(lla, datevec(obs1time)); %[m]

% Nutation Transform
N = nutation(datevec(obs1time));

% Precession Transform
P = precession(JD1);

% Rotate to J2000
station_vec1 = P' * N' * pos_eci1';

[az1,h1,ra1,dec1] = getAzElRaDec(JD1,station_vec1', Z(end-offset,1:3));

% Compute finite difference
[~,~,ra1_FD_1,dec1_FD_1] = getAzElRaDec(JD1,station_vec1', Z(end-offset-1,1:3));
[~,~,ra1_FD_2,dec1_FD_2] = getAzElRaDec(JD1,station_vec1', Z(end-offset+1,1:3));

% Right Ascension finite difference
ra1_dot_1 = (ra1 - ra1_FD_1);
ra1_dot_2 = (ra1_FD_2 - ra1);
ra1_dot = mean([ra1_dot_1,ra1_dot_2]) * pi/180; % Convert to rad/s

% Declination finite difference
dec1_dot_1 = (dec1 - dec1_FD_1);
dec1_dot_2 = (dec1_FD_2 - dec1);
dec1_dot = mean([dec1_dot_1,dec1_dot_2]) * pi/180; % Convert to rad/s

% ra1_dot = 3.1893e-04;
% dec1_dot =  -7.9012e-05;
% ra1 = -0.6112;
% dec1 = 0.4954;

% Compute velocity of station in ECI
R_dot_ECI = cross([0 0 earth_rot],station_vec1);

num_orbits = 1;
P_coast = num_orbits * (2*pi)/(sqrt(mu/(a^3))); %Coast time
t = 1:1:P_coast;

e = 0.4;

[ecc_1,ecc_2,ecc_3, ecc_4, a_upper, a_lower, zero_upper, zero_lower, rho, orbit_rand, x_random, y_random] = get_a_ecc_constraints(station_vec1,R_dot_ECI,ra1*pi/180,ra1_dot,dec1*pi/180,dec1_dot,max_rho_guess,t,e);

% [rho_dot_test rho_test] = eccentricityConstraint ( ra1, dec1, 0, 0, station_vec1, R_dot_ECI);

% OBSERVATION 2 =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
obs2time = utc_datetime + seconds(Tout(end));
JD2 = juliandate(obs2time);

% Compute station ECI at end of orbit propogation
pos_eci2 = lla2eci(lla, datevec(obs2time));

% Nutation Transform
N = nutation(datevec(obs2time));

% Precession Transform
P = precession(JD2);

% Rotate to J2000
station_vec2 = P' * N' * pos_eci2';

[az2,h2,ra2,dec2] = getAzElRaDec(JD2,station_vec2',Z(end,1:3));

% Compute finite difference
[~,~,ra2_FD_1,dec2_FD_1] = getAzElRaDec(JD2,station_vec2', Z(end-1,1:3));

% Right Ascension finite difference
ra2_dot_1 = (ra2 - ra2_FD_1);
ra2_dot = mean(ra2_dot_1) * pi/180; % Convert to rad/s

% Declination finite difference
dec2_dot_1 = (dec2 - dec2_FD_1);
dec2_dot = mean(dec2_dot_1) * pi/180; % Convert to rad/s

[ecc_1_obs2,ecc_2_obs2, ecc_3_obs2, ecc_4_obs2, a_upper_obs2, a_lower_obs2, zero_upper_obs2, zero_lower_obs2, rho_obs2, orbit_rand_obs2, x_random_obs2, y_random_obs2] = get_a_ecc_constraints(station_vec2,R_dot_ECI,ra2*pi/180,ra2_dot,dec1*pi/180,dec2_dot,max_rho_guess,t,e);


%% Observation 1 Admissible Region  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
close all
a = figure(1);
hold on
grid minor
plot(NaN,NaN,'r','LineWidth',6);
plot(NaN,NaN,'b','LineWidth',6);
plot(NaN,NaN,'k','LineWidth',6);
scatter(NaN,NaN,250,'p','filled','MarkerEdgeColor','k','MarkerFaceColor','y');
plot(rho/RE,zero_upper./1000,'r')
plot(rho/RE,zero_lower./1000,'r')

plot(rho/RE,a_upper./1000,'b')
plot(rho/RE,a_lower./1000,'b')

plot(rho/RE,ecc_1./1000,'k')
plot(rho/RE,ecc_2./1000,'k')
plot(rho/RE,ecc_3./1000,'k')
plot(rho/RE,ecc_4./1000,'k')

% plot(rho_dot_test(1,:)./1000,'k')
% plot(rho_dot_test(2,:)./1000,'k')
% plot(rho_dot_test(3,:)./1000,'k')
% plot(rho_dot_test(4,:)./1000,'k')

rho_truth = norm(Z(end-offset,1:3)' - station_vec1);
rho_dot_truth = norm(Z(end-offset,4:6) - R_dot_ECI);
scatter(rho_truth/RE,rho_dot_truth/1000,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor','y')

scatter(x_random/RE,y_random./1000)

% plot(drho)
xlabel('Range [# of Earth Radii]','FontSize',25)
ylabel('Range rate [km/s]','FontSize',25)
legend('Zero-energy Constraint','Semi-Major Axis Constraint','Eccentricity Constraint','Truth Orbit','FontSize',25)

a.Position = [100 100 1000 1000];

%% Orbit Plot  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
b = figure(3);
hold on
imData = imread('2_no_clouds_4k.jpg');  % Load the texture image

[xS, yS, zS] = sphere(50);  % Generate a sphere for the globe
earth_radius = 6378.1370 * 1000;  % Earth radius in meters

xSE = earth_radius * xS;
ySE = earth_radius * yS;
zSE = earth_radius * zS;

% Create the rotation matrices
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
surface(xSE_rot, ySE_rot, zSE_rot, 'FaceColor', 'texturemap', 'CData', flipud(imData), 'EdgeColor', 'none');

% Adjust axes
axis equal;
view(3);  % 3D view
grid on;
xlabel('Inertial x (m)');
ylabel('Inertial y (m)');
zlabel('Inertial z (m)');

plot3(Z(:,1), Z(:,2), Z(:,3), 'r', 'LineWidth', 4);

j = 1;
for i=1:length(x_random)
    scatter3(orbit_rand(1,j), orbit_rand(1,j+1), orbit_rand(1,j+2) , 1000, 'k.');
    plot3(orbit_rand(:,j), orbit_rand(:,j+1), orbit_rand(:,j+2), 'LineWidth', 0.5);
    scatter3(orbit_rand_obs2(1,j), orbit_rand_obs2(1,j+1), orbit_rand_obs2(1,j+2) , 1000, 'r.');
    j = j+6;
end

% Plot observer location in ECI
scatter3(pos_eci1(1), pos_eci1(2), pos_eci1(3), 1000, 'c.');
scatter3(pos_eci2(1), pos_eci2(2), pos_eci2(3), 1000, 'm.');

% Mark the end point of the orbit
plot3(Z(end-offset,1), Z(end-offset,2), Z(end-offset,3), 'pm', 'MarkerSize', 20,'MarkerFaceColor','r','MarkerEdgeColor','k');
plot3(Z(end,1), Z(end,2), Z(end,3), 'pm', 'MarkerSize', 20,'MarkerFaceColor','r','MarkerEdgeColor','k');

% Get axis children (just as in the original code)
ch = get(gca, 'children');

b.Position = [1600 100 1000 1000];



















