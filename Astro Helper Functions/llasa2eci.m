function rECI = llasa2eci(lat, lon, alt, gst)

   % This function is buggy, debug it
   
    % INPUTS:
    %   lat - Geodetic latitude [degrees]
    %   lon - Longitude [degrees] (east positive)
    %   alt - Altitude above WGS-84 ellipsoid [meters]
    %   gst - Greenwich Sidereal Time [degrees] (at observation time)
    
    % OUTPUT:
    %   rECI - 3x1 vector [x; y; z] in ECI coordinates [meters]
   
    % WGS-84 ellipsoid parameters
    a = 6378137.0;             % Semi-major axis [m]
    f = 1 / 298.257223563;     % Flattening
    e2 = 2*f - f^2;            % First eccentricity squared

    % Convert to radians
    lat = deg2rad(lat);
    lon = deg2rad(lon);
    gst = deg2rad(gst); % Convert GST from degrees to radians

    % Compute prime vertical radius of curvature
    sinLat = sin(lat);
    N = a / sqrt(1 - e2 * sinLat^2);

    % Compute ECEF coordinates
    xECEF = (N + alt) * cos(lat) * cos(lon);
    yECEF = (N + alt) * cos(lat) * sin(lon);
    zECEF = (N * (1 - e2) + alt) * sin(lat);
    rECEF = [xECEF; yECEF; zECEF];

    % Rotation matrix from ECEF to ECI (about Z-axis by Greenwich Sidereal Time)
    Rz = [ cos(gst), sin(gst), 0;
          -sin(gst), cos(gst), 0;
                0,        0, 1];

    % Compute ECI coordinates
    rECI = Rz * rECEF;

end