function [lat, lon, h] = ecef2geodetic(x, y, z)

    % INPUTS:
    %   x, y, z - ECEF coordinates in meters (can be vectors)
    
    % OUTPUTS:
    %   lat - Geodetic latitude in degrees
    %   lon - Longitude in degrees
    %   h   - Altitude above the WGS-84 ellipsoid in meters

    % WGS-84 ellipsoid parameters
    a = 6378137.0;                 % Semi-major axis [m]
    f = 1 / 298.257223563;          % Flattening
    e2 = 2*f - f^2;                 % Square of first eccentricity

    % Compute longitude
    lon = atan2d(y, x);  % Longitude in degrees

    % Initialize iteration for geodetic latitude
    p = sqrt(x.^2 + y.^2);  % Projection onto XY plane
    lat = atan2d(z, p * (1 - e2));  % Initial estimate of latitude
    h = 0;  % Initial altitude

    % Iterative solution for latitude and altitude
    max_iter = 10;
    tol = 1e-6;  % Convergence tolerance (about 1 mm)

    for iter = 1:max_iter
        sinLat = sind(lat);
        N = a / sqrt(1 - e2 * sinLat.^2);  % Prime vertical radius
        h_old = h;
        h = p ./ cosd(lat) - N;
        lat = atan2d(z + e2 * N * sinLat, p);
        
        % Check for convergence
        if max(abs(h - h_old)) < tol
            break;
        end
    end

end