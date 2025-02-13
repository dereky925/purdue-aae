clear;clc;close all

mu = 3.9860044e14; %m^3/s^2

% Speed of light [m/s]
c = 3E8;

% LLA of Armstrong Hall
lat = 40.431;
lon = -86.915;
alt = 10;

lla = [lat, lon, alt];

% Starlink satellite pos vel vector
% space-track.org
% SpaceX_Ephemeris_552_SpaceX_2024-09-16UTC05_21_07_01
r_ijk = [-3888.1214732056 2828.0607278016 -4851.8696685028]' * 1000; % [m]
v_ijk = [-1.5642842726 -6.9319980689 -2.7897628635]' * 1000; % [m/s]

% 2024-09-19 02:17:42
utc_datetime = datetime(2024, 9, 16, 2, 17, 42, 'TimeZone', 'UTC');

[a,ecc,incl,RAAN,argp,nu,truelon,arglat,lonper] = ijk2keplerian(r_ijk, v_ijk);

num_orbits = 1.55;
P_coast = num_orbits * (2*pi)/(sqrt(mu/(a^3))); %Coast time
t = 1:1:P_coast;

options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
[Tout, Z] = ode45(@two_body_ode,t,[r_ijk v_ijk],options);


% OBSERVATION 1 =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
offset = 300;
obs1time = utc_datetime + seconds(Tout(end-offset));
JD1 = juliandate(obs1time);

% Compute station ECI 5 minutes before end of orbit propogation
pos_eci1 = lla2eci(lla, datevec(obs1time));

% Nutation Transform
N = nutation(datevec(obs1time));

% Precession Transform
P = precession(JD1);

% Rotate to J2000
station_vec1 = P' * N' * pos_eci1';

[az1,h1,ra1,dec1] = getAzElRaDec(JD1,station_vec1', Z(end-offset,1:3))

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

[az2,h2,ra2,dec2] = getAzElRaDec(JD2,station_vec2',Z(end,1:3))

% Initial orbit determination - Restricted IOD =-=-=-=-=-=-=-=-=-=-=-=-=-

L_t1 = [ cosd(ra1)*cosd(dec1);
      sind(ra1)*cosd(dec1);
            sind(dec1)];

L_t2 = [ cosd(ra2)*cosd(dec2);
      sind(ra2)*cosd(dec2);
            sind(dec2)];


% Assume e = 0 and w = 0

% Guess values of a (semi-major axis)
R_Earth = 6.3781E6;
step = 100;
a_guess = R_Earth : step : R_Earth*1.2;

rho_t1 = -L_t1'*station_vec1 + sqrt( (L_t1'*station_vec1)^2 - ( norm(station_vec1)^2 - a_guess.^2 ) );
rho_t2 = -L_t2'*station_vec2 + sqrt( (L_t2'*station_vec2)^2 - ( norm(station_vec2)^2 - a_guess.^2 ) );

r_t1 = rho_t1'*L_t1' + station_vec1';
r_t2 = rho_t2'*L_t2' + station_vec2';

phi_g = acos( dot(r_t1',r_t2') ./ ( vecnorm(r_t1,2,2).*vecnorm(r_t2,2,2))' );

phi_CM = sqrt(mu/a^3) * (offset - rho_t2/c + rho_t1/c);

F = phi_CM - phi_g;

for i=2:length(a_guess)
    if abs(F(i))<abs(F(i-1))
        F_ans = F(i);
        index = i;
    end
end

a_IOD = R_Earth + index*step;

figA = figure(1);
hold on
grid minor
scatter(NaN,NaN, 500, 'pm','MarkerFaceColor','y','MarkerEdgeColor','k')
plot(F)
scatter(index,F_ans,500, 'pm','MarkerFaceColor','y','MarkerEdgeColor','k')
figA.Position = [100, 100, 1000, 1000];
xlabel('a (semi-major axis) [m]','FontSize',15)
ylabel('F Cost Function','FontSize',15)
title('F Cost Function vs Semi-major axis','FontSize',15)
legend('a = 6814800 [m]','FontSize',25);

% Compute angular momentum
rho_t1_IOD = -L_t1'*station_vec1 + sqrt( (L_t1'*station_vec1)^2 - ( norm(station_vec1)^2 - a_IOD.^2 ) );
rho_t2_IOD = -L_t2'*station_vec2 + sqrt( (L_t2'*station_vec2)^2 - ( norm(station_vec2)^2 - a_IOD.^2 ) );

r_t1_IOD = rho_t1_IOD'*L_t1' + station_vec1';
r_t2_IOD = rho_t2_IOD'*L_t2' + station_vec2';
angular_momentum = cross(r_t1_IOD,r_t2_IOD);

Omega_IOD = atan2(angular_momentum(1),-angular_momentum(2));
i_IOD = acos(angular_momentum(3)/norm(angular_momentum));

r_Omega_t1 = R1(i_IOD)*R3(Omega_IOD)*r_t1_IOD';
r_Omega_t2 = R1(i_IOD)*R3(Omega_IOD)*r_t2_IOD';

v_IOD = atan2d(norm(r_Omega_t2),norm(r_Omega_t1));

w_IOD = 0;
e_IOD = 0.001;
i_IOD = i_IOD*180/pi;
Omega_IOD = 360 + Omega_IOD*180/pi;

% [a,ecc,incl,RAAN,argp,nu]
RIOD_15min = [(a_IOD-a)/a*100 (e_IOD-ecc)/ecc*100 (i_IOD-incl)/incl*100 (Omega_IOD-RAAN)/RAAN*100 (w_IOD-argp)/argp*100 (v_IOD-nu)/nu*100];

[r_ijk_IOD, v_ijk_IOD] = keplerian2ijk(a_IOD, e_IOD, i_IOD, Omega_IOD, w_IOD, v_IOD);

[Tout_IOD, Z_IOD] = ode45(@two_body_ode,t,[r_ijk_IOD v_ijk_IOD],options);

% Gibbs' method =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% Take 3 position measurements from propogation

r_1_Gibbs = [Z(end-offset*2,1), Z(end-offset*2,2), Z(end-offset*2,3)];
r_2_Gibbs = [Z(end-offset,1), Z(end-offset,2), Z(end-offset,3)];
r_3_Gibbs = [Z(end,1), Z(end,2), Z(end,3)];

plane_verification = (r_1_Gibbs/norm(r_1_Gibbs)) * (cross(r_2_Gibbs,r_3_Gibbs)/norm(cross(r_2_Gibbs,r_3_Gibbs)))';

n = norm(r_1_Gibbs) * cross(r_2_Gibbs,r_3_Gibbs) + norm(r_2_Gibbs)*cross(r_3_Gibbs,r_1_Gibbs) + norm(r_3_Gibbs)*cross(r_1_Gibbs,r_2_Gibbs);
d = cross(r_1_Gibbs,r_2_Gibbs) + cross(r_2_Gibbs,r_3_Gibbs) + cross(r_3_Gibbs,r_1_Gibbs);
s = r_1_Gibbs*(norm(r_2_Gibbs)-norm(r_3_Gibbs)) + r_2_Gibbs*(norm(r_3_Gibbs)-norm(r_1_Gibbs)) + r_3_Gibbs*(norm(r_1_Gibbs)-norm(r_2_Gibbs));

v2 = sqrt( mu/(norm(n)*norm(d)) ) * ( (cross(d,r_2_Gibbs)/norm(r_2_Gibbs)) + s);


[a_Gibbs,ecc_Gibbs,incl_Gibbs,RAAN_Gibbs,argp_Gibbs,nu_Gibbs,truelon,arglat,lonper] = ijk2keplerian(r_2_Gibbs, v2)
Gibbs_5min = [a_Gibbs-a ecc_Gibbs-ecc incl_Gibbs-incl RAAN_Gibbs-RAAN argp_Gibbs-argp nu_Gibbs-nu];

[Tout_Gibbs, Z_Gibbs] = ode45(@two_body_ode,t,[r_2_Gibbs v2],options);

% Orbit Plot - Restricted IOD =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
b = figure(2);
imData = imread('2_no_clouds_4k.jpg');  % Load the texture image

[xS, yS, zS] = sphere(50);  % Generate a sphere for the globe
earth_radius = 6378.1370 * 1000;  % Earth radius in meters

xSE = earth_radius * xS;
ySE = earth_radius * yS;
zSE = earth_radius * zS;

% Create the rotation matrices
Ry = R3(-lon-0.15);

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


% Plot the orbit on top of the globe
hold on;
plot3(Z(:,1), Z(:,2), Z(:,3), 'r', 'LineWidth', 4);
plot3(Z_IOD(:,1), Z_IOD(:,2), Z_IOD(:,3), ':g', 'LineWidth', 4);

% Plot observer location in ECI
scatter3(pos_eci1(1), pos_eci1(2), pos_eci1(3), 1000, 'c.');
scatter3(pos_eci2(1), pos_eci2(2), pos_eci2(3), 1000, 'm.');

% Mark the end point of the orbit
plot3(Z(end-offset,1), Z(end-offset,2), Z(end-offset,3), 'pm', 'MarkerSize', 20,'MarkerFaceColor','r','MarkerEdgeColor','k');
plot3(Z(end,1), Z(end,2), Z(end,3), 'pm', 'MarkerSize', 20,'MarkerFaceColor','r','MarkerEdgeColor','k');

% Mark the estimated true anomaly of IOD
% plot3(Z_IOD(1,1), Z_IOD(1,2), Z_IOD(1,3), 'pm', 'MarkerSize', 20,'MarkerFaceColor','g','MarkerEdgeColor','k');

% Get axis children (just as in the original code)
get(gca, 'children');

b.Position = [1500 100 1000 1000];

% Orbit Plot - Gibbs =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
b = figure(3);
imData = imread('2_no_clouds_4k.jpg');  % Load the texture image

[xS, yS, zS] = sphere(50);  % Generate a sphere for the globe
earth_radius = 6378.1370 * 1000;  % Earth radius in meters

xSE = earth_radius * xS;
ySE = earth_radius * yS;
zSE = earth_radius * zS;


% Create the rotation matrices
Ry = R3(-lon-0.15);

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

% Plot the orbit on top of the globe
hold on;
plot3(Z(:,1), Z(:,2), Z(:,3), 'r', 'LineWidth', 4);
plot3(Z_Gibbs(:,1), Z_Gibbs(:,2), Z_Gibbs(:,3), ':g', 'LineWidth', 4);

% Plot observer location in ECI
scatter3(pos_eci1(1), pos_eci1(2), pos_eci1(3), 1000, 'c.');
scatter3(pos_eci2(1), pos_eci2(2), pos_eci2(3), 1000, 'm.');

% Mark the end point of the orbit
plot3(Z(end-offset*2,1), Z(end-offset*2,2), Z(end-offset*2,3), 'pm', 'MarkerSize', 20,'MarkerFaceColor','r','MarkerEdgeColor','k');
plot3(Z(end-offset,1), Z(end-offset,2), Z(end-offset,3), 'pm', 'MarkerSize', 20,'MarkerFaceColor','r','MarkerEdgeColor','k');
plot3(Z(end,1), Z(end,2), Z(end,3), 'pm', 'MarkerSize', 20,'MarkerFaceColor','r','MarkerEdgeColor','k');

% Get axis children (just as in the original code)
ch = get(gca, 'children');

b.Position = [1600 100 1000 1000];




