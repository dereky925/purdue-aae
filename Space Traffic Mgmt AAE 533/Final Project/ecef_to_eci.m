function r_eci = ecef_to_eci(r_ecef, utc_time)
    % Inputs:
    % r_ecef: [3x1] ECEF position vector (meters)
    % utc_time: datetime object representing the time in UTC
    %
    % Output:
    % r_eci: [3x1] ECI position vector (meters)

    % Constants
    omega_earth = 7.2921150e-5; % Earth's rotation rate (rad/s)
    
    % Step 1: Compute Julian Date (JD)
    jd = juliandate(utc_time); % MATLAB's juliandate function
    
    % Step 2: Calculate Greenwich Sidereal Time (GST) in radians
    T_ut1 = (jd - 2451545.0) / 36525; % Centuries since J2000.0
    theta_gst = mod(280.46061837 + 360.98564736629 * (jd - 2451545.0) ...
                    + 0.000387933 * T_ut1^2 - T_ut1^3 / 38710000, 360);
    theta_gst_rad = deg2rad(theta_gst); % Convert degrees to radians

    % Step 3: Rotation matrix about z-axis
    R_z = [cos(theta_gst_rad), sin(theta_gst_rad), 0;
          -sin(theta_gst_rad), cos(theta_gst_rad), 0;
           0,                  0,                 1];

    % Step 4: Transform ECEF to ECI
    r_eci = R_z * r_ecef;
end