function [lat_deg, lon_deg, alt_km] = ecef2ellipsoid(x_km, y_km, z_km)
% ECEF2ELLIPSOIDALT  Converts ECEF coordinates (x, y, z) [km]
% to geodetic latitude (deg), longitude (deg), and altitude (km)
% above the WGS-84 reference ellipsoid.
%
% [lat_deg, lon_deg, alt_km] = ecef2ellipsoidAlt(x_km, y_km, z_km)
%
% Inputs:
%   x_km, y_km, z_km  - Satellite position in ECEF coordinates [km]
%
% Outputs:
%   lat_deg  - Geodetic latitude [degrees]
%   lon_deg  - Longitude [degrees]
%   alt_km   - Altitude above the WGS-84 ellipsoid [km]
%
% Reference:
%   WGS-84 constants: 
%       a = 6378.137 km (semi-major axis)
%       f = 1/298.257223563 (flattening)
%       e^2 = 2f - f^2 (eccentricity squared)

    % Define WGS-84 parameters
    a = 6378.137;               % semi-major axis [km]
    f = 1/298.257223563;        % flattening
    e2 = 2*f - f^2;             % eccentricity squared


    % Earth-fixed radius in equatorial plane
    r_xy = sqrt(x_km^2 + y_km^2); 
    % Longitude
    lon_rad = atan2(y_km, x_km);

    % An initial guess for lat (geocentric approach or a simpler method)
    % Use the "Bowring's formula" initial guess or a simpler approach:
    if (r_xy < 1e-12)
        % Very close to Z-axis (near poles)
        lon_rad = 0.0; 
        if z_km > 0, lat_rad =  pi/2; else, lat_rad = -pi/2; end
    else
        lat_rad = atan2(z_km, r_xy*(1 - e2));  
    end

    % Iterate to refine lat_rad
    % typical tolerance for altitude calculations
    tolerance = 1e-9;  
    diff = 1.0;  

    % Iteration for geodetic latitude

    while abs(diff) > tolerance
        sinLat = sin(lat_rad);
        N = a / sqrt(1 - e2*sinLat^2);   % prime vertical radius of curvature
        alt_tmp = r_xy/cos(lat_rad) - N; 
        % altitude used to refine lat
        lat_new = atan2(z_km, r_xy*(1 - e2*N/(N + alt_tmp)));
        diff = lat_new - lat_rad;
        lat_rad = lat_new;
    end

    % After convergence:
    sinLat = sin(lat_rad);
    N = a / sqrt(1 - e2*sinLat^2);
    alt_km = r_xy/cos(lat_rad) - N;

    % Convert to degrees
    lat_deg = rad2deg(lat_rad);
    lon_deg = rad2deg(180/pi);
end