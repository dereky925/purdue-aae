function [az,h,ra,dec] = getAzElRaDec(JD,station_position_eci,satellite_position_eci)

% Outputs
% az - Azimuth [degrees]
% h - Elevation [degrees]
% ra - Right Ascension [degrees]
% dec - Declination [degrees]

station_x_eci = station_position_eci(1);
station_y_eci = station_position_eci(2);
station_z_eci = station_position_eci(3);

sat_x_eci = satellite_position_eci(1);
sat_y_eci = satellite_position_eci(2);
sat_z_eci = satellite_position_eci(3);

% Compute geocentric latitude
geoc_lat = 90 - acosd (station_z_eci ./( sqrt(station_x_eci .^2 + station_y_eci .^2 + station_z_eci .^2)));

utc_datetime = datevec(datetime(JD, 'convertfrom', 'juliandate', 'TimeZone', 'UTC'));
lla = eci2lla(station_position_eci,utc_datetime);

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
sidereal_time = mod( UT + lla(2)./15 , 24);

% % Compute sidereal angle [deg]
sidereal_angle = sidereal_time * 15;

% Compute range magnitude to satellite
r = sqrt( (sat_x_eci).^2 ...
        + (sat_y_eci).^2 ...
        + (sat_z_eci).^2 ) ;

r_xy = sqrt(sat_x_eci.^2 + sat_y_eci.^2);

% Satellite Declination, geocentric 
sat_decl_geoc = atan2d(sat_z_eci, r_xy);

% Satellite Right Ascension, geocentric
sat_r_ascen_geoc = atan2d(sat_y_eci, sat_x_eci);

% Compute range magnitude to station in ECI
R_eci = sqrt(station_x_eci.^2 + station_y_eci.^2 + station_z_eci.^2);

% Compute range magnitude from station to satellite
rho = sqrt( (sat_x_eci - station_x_eci)^2 + (sat_y_eci - station_y_eci)^2 + (sat_z_eci - station_z_eci)^2 );

% Compute topocentric declination [deg]
dec = asind( (r.*sind(sat_decl_geoc) ...
            - R_eci.*sind(geoc_lat)) ...
            ./ rho );

% Compute topocentric right ascension [deg]
r_ascen_topo_y = ( (r.*cosd(sat_decl_geoc).*sind(sat_r_ascen_geoc) ...
                       - R_eci.*cosd(geoc_lat).*sind(sidereal_angle)) ...
                       ./ (rho.*cosd(dec)) );

r_ascen_topo_x = ( (r.*cosd(sat_decl_geoc).*cosd(sat_r_ascen_geoc) ...
                        -R_eci.*cosd(geoc_lat).*cosd(sidereal_angle)) ...
                        ./ (rho.*cosd(dec)));

% Topocentric right ascension
ra = atan2d(r_ascen_topo_y,r_ascen_topo_x);


% Compute hour angle
tau = sidereal_angle - ra;

% Compute elevation
h = asind(sind(lla(1)).*sind(dec) ...
          + cosd(lla(1)).*cosd(dec).*cosd(tau));

% Compute azimuth
az_x = ( (sind(lla(1)).*cosd(dec).*cosd(tau) ...
    - cosd(lla(1)).*sind(dec)) ...
    ./ cosd(h) );

az_y = (cosd(dec).*sind(tau)./cosd(h));

az = atan2d(az_y, az_x);

% ra = sat_r_ascen_geoc;
% 
% dec = sat_decl_geoc;

end
