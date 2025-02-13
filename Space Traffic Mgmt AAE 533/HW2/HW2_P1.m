clear;clc;close all

% Speed of light [m/s]
c = 3E8;

% LLA of Armstrong Hall
lat = 40.431;
lon = -86.915;
alt = 0;

% Put LLA into vector
lla = [lat,lon,alt];

% Compute ECI XYZ from LLA and UTC time, Sep 6, 0:30 AM 2024 [m]
pos_eci = lla2eci(lla, [2024 9 6 4 30 00]);

% Compute geocentric latitude
% geoc_lat = atan2d(pos_eci(2),pos_eci(1));
geoc_lat = 90 - acosd (pos_eci(3) ./( sqrt(pos_eci(1) .^2 + pos_eci(2) .^2 + pos_eci(3) .^2)));

% Add pi to geocentric latitude to get geocentric RA, eqn 23 in notes
geoc_RA_deg = (pi + geoc_lat) * 180/pi;

mu = 3.9860044e14; %m^3/s^2

% Object Parameters =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

a = 6782.5E3; % Semi-major axis [m] GEO
e = 0.001; % Eccentricity
i = 51.5; % Inclination [deg]
RAAN = 10; % Right Ascension of Ascending Node [deg]
w = 20; % Argument of Perigee [deg]
v = 5; % True Anomaly [deg]

% noon UTC = 8 AM UTC
P_coast = 16.5 * 3600; % Propogate orbit for 16.5 hours to Sep 6 0:30 AM
t = 1:1:P_coast;

[r_ijk,v_ijk] = keplerian2ijk(a,e,i,RAAN,w,v);

options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
[Tout, Z] = ode45(@two_body_ode,t,[r_ijk v_ijk],options);

satFinalPosECI = Z(end,1:3);

station_x_eci = pos_eci(1); % [m]
station_y_eci = pos_eci(2); % [m]
station_z_eci = pos_eci(3); % [m]

% Compute range magnitude to station in ECI
R_eci = sqrt(station_x_eci.^2 + station_y_eci.^2 + station_z_eci.^2);

r_xy = sqrt(satFinalPosECI(1).^2 + satFinalPosECI(2).^2);

% Satellite Declination, geocentric 
sat_decl_geoc = atan2d(satFinalPosECI(3), r_xy);

% Satellite Right Ascension, geocentric
sat_r_ascen_geoc = atan2d(satFinalPosECI(2), satFinalPosECI(1));

% ====== Compute Julian Date ==========================================

% Define the date and time
year = 2024;
month = 9;
day = 6;
hour = 4;
minute = 30;

% If the month is January or February, adjust year and month
if month <= 2
    year = year - 1;
    month = month + 12;
end

% Calculate A and B
A = floor(year / 100);
B = 2 - A + floor(A / 4);

% Calculate the Julian Day
JD = floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + day + B - 1524.5;

% Adjust for the time in hours and minutes (convert to fractional day)
fractional_day = (hour + minute / 60) / 24;
JD = JD + fractional_day;


% Compute Julian Centuries, T1
T1 = (JD - 2451545) ./ 36525; 

% Compute JD0, keep only integer portion
JD0 = floor(JD); 

% Compute T0
T0 = (JD0 - 2451545) ./ 36525; 

% Mean sidereal time of Greenwich 0h Earth time (UT)
UT = 24110.54841 + 8640184.812866*T0 + 0.093104*T1.^2 - 0.0000062*T1.^3 + (1.0027279093*(mod(JD,1)-.5).*24).*3600;
UT = UT ./ 3600;

% Compute sidereal time
sidereal_time = mod( UT + lon./15 , 24);

% % Compute sidereal angle [deg]
sidereal_angle = sidereal_time * 15;

% Compute range magnitude to satellite
r = sqrt( (satFinalPosECI(1)).^2 ...
        + (satFinalPosECI(2)).^2 ...
        + (satFinalPosECI(3)).^2 ) ;

% Compute range magnitude from station to satellite
rho = sqrt( (satFinalPosECI(1) - station_x_eci)^2 + (satFinalPosECI(2) - station_y_eci)^2 + (satFinalPosECI(3) - station_z_eci)^2 );

station_vec = [ station_x_eci station_y_eci station_z_eci];

rho2 = vecnorm(satFinalPosECI - station_vec,2,2);


% Nutation Transform
N = nutation([2024,9,6,4,30]);

% Precession Transform
P = precession(JD);

% Rotate to J2000
station_vec = P' * N' * station_vec';

% satFinalPosECI = P' * N' * satFinalPosECI'

% Compute topocentric declination [deg]
decl_topo = asind( (r.*sind(sat_decl_geoc) ...
            - R_eci.*sind(geoc_lat)) ...
            ./ rho )

% Compute topocentric right ascension [deg]
r_ascen_topo_y = ( (r.*cosd(sat_decl_geoc).*sind(sat_r_ascen_geoc) ...
                       - R_eci.*cosd(geoc_lat).*sind(sidereal_angle)) ...
                       ./ (rho.*cosd(decl_topo)) );

r_ascen_topo_x = ( (r.*cosd(sat_decl_geoc).*cosd(sat_r_ascen_geoc) ...
                        -R_eci.*cosd(geoc_lat).*cosd(sidereal_angle)) ...
                        ./ (rho.*cosd(decl_topo)));

r_ascen_topo = atan2d(r_ascen_topo_y,r_ascen_topo_x)

% Compute hour angle
tau = sidereal_angle - r_ascen_topo;

% Compute elevation
h = asind(sind(lat).*sind(decl_topo) ...
          + cosd(lat).*cosd(decl_topo).*cosd(tau))

% Compute azimuth
az_x = ( (sind(lat).*cosd(decl_topo).*cosd(tau) ...
    - cosd(lat).*sind(decl_topo)) ...
    ./ cosd(h) );

az_y = (cosd(decl_topo).*sind(tau)./cosd(h));

az = atan2d(az_y, az_x)

% [EL,AZ] = getEleAzi(lat,decl_topo,tau)

[az,h,ra,dec] = getAzElRaDec(JD,station_vec',satFinalPosECI)

%%
% 1.b.=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

% Fixed point iteration of tau
dt = 0.000001;
t_small_steps = 0:dt:0.1;

% Propogate orbit backwards from final satellite position
[Tout_small, orbit_small_t] = ode45(@two_body_ode,t_small_steps,[satFinalPosECI' -v_ijk],options);

% Fixed-point iteration to solve for tau, signal travel time [s]
travel_time = 0.001;
for i=1:10
    % Compute station location as function of time
    station_vec = lla2eci(lla, [2024 9 6 4 29 (60-travel_time)]);

    % Position of satellite, signal-travel-time-seconds-ago
    r_tau = orbit_small_t( round(round(travel_time,5)/dt + 1) ,1:3);

    % Compute signal travel time
    travel_time = 1/c*( vecnorm( r_tau -  station_vec, 2, 2)  );
end


satFinalPosECI = r_tau;

% Compute range magnitude from station to satellite
rho = sqrt( (satFinalPosECI(1) - station_x_eci)^2 + (satFinalPosECI(2) - station_y_eci)^2 + (satFinalPosECI(3) - station_z_eci)^2 );

station_vec = [ station_x_eci station_y_eci station_z_eci];

rho2 = vecnorm(satFinalPosECI - station_vec,2,2);


% Nutation Transform
N = nutation([2024,9,6,4,30]);

% Precession Transform
P = precession(JD);

% Rotate to J2000
station_vec = P' * N' * station_vec';

% satFinalPosECI = P' * N' * satFinalPosECI'

% Compute topocentric declination [deg]
decl_topo = asind( (r.*sind(sat_decl_geoc) ...
            - R_eci.*sind(geoc_lat)) ...
            ./ rho );

% Compute topocentric right ascension [deg]
r_ascen_topo_y = ( (r.*cosd(sat_decl_geoc).*sind(sat_r_ascen_geoc) ...
                       - R_eci.*cosd(geoc_lat).*sind(sidereal_angle)) ...
                       ./ (rho.*cosd(decl_topo)) );

r_ascen_topo_x = ( (r.*cosd(sat_decl_geoc).*cosd(sat_r_ascen_geoc) ...
                        -R_eci.*cosd(geoc_lat).*cosd(sidereal_angle)) ...
                        ./ (rho.*cosd(decl_topo)));

r_ascen_topo = atan2d(r_ascen_topo_y,r_ascen_topo_x);

% Compute hour angle
tau = sidereal_angle - r_ascen_topo;

% Compute elevation
h = asind(sind(lat).*sind(decl_topo) ...
          + cosd(lat).*cosd(decl_topo).*cosd(tau))

% Compute azimuth
az_x = ( (sind(lat).*cosd(decl_topo).*cosd(tau) ...
    - cosd(lat).*sind(decl_topo)) ...
    ./ cosd(h) );

az_y = (cosd(decl_topo).*sind(tau)./cosd(h));

az = atan2d(az_y, az_x)



% Orbit Plot =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
imData = imread('2_no_clouds_4k.jpg');  % Load the texture image

[xS, yS, zS] = sphere(50);  % Generate a sphere for the globe
earth_radius = 6378.1370 * 1000;  % Earth radius in meters

xSE = earth_radius * xS;
ySE = earth_radius * yS;
zSE = earth_radius * zS;


% Create the rotation matrices
Ry = R3(-lon+0.1);

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
plot3(Z(:,1), Z(:,2), Z(:,3), 'b', 'LineWidth', 4);

% Plot observer location in ECI
scatter3(pos_eci(1), pos_eci(2), pos_eci(3), 500, 'r.');

% Mark the end point of the orbit
plot3(Z(end,1), Z(end,2), Z(end,3), '*r', 'MarkerSize', 10);

% Get axis children (just as in the original code)
ch = get(gca, 'children');






