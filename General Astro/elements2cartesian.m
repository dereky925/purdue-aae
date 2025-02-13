function [r_vec, v_vec] = elements2cartesian(a, ex, ey, incl, raan, theta)
% ELEMENTS2CARTESIAN Convert orbital elements to ECI state vectors.
%
%   [r_vec, v_vec] = elements2cartesian(a, ex, ey, incl, raan, theta)
%
%   Inputs:
%       a     = Semi-major axis [km]
%       ex    = x-component of eccentricity (unitless)
%       ey    = y-component of eccentricity (unitless)
%       incl  = Inclination [°]
%       raan  = Right ascension of ascending node [°]
%       theta = Argument of latitude [°]  (argp + true anomaly, mod 360)
%
%   Outputs:
%       r_vec    = ECI position vector [km] (3x1)
%       v_vec    = ECI velocity vector [km/s] (3x1)
%
%   This function reconstructs the scalar eccentricity and argument of perigee,
%   computes the true anomaly, and then converts the Keplerian elements to
%   Cartesian coordinates using standard formulas.
%
%   Author: Derek Yu (adapted)
%   Date: Jan 21, 2025

    % Gravitational parameter [km^3/s^2]
    mu = 398600.4418;  

    % Compute scalar eccentricity and argument of perigee (in degrees)
    ecc = sqrt(ex^2 + ey^2);
    argp = rad2deg(atan2(ey, ex));  % This returns values in [-180, 180]
    % Adjust to [0, 360)
    if argp < 0
        argp = argp + 360;
    end

    % Compute true anomaly from the argument of latitude
    % (theta = argp + nu modulo 360) so that nu = mod(theta - argp, 360)
    nu = mod(theta - argp, 360);

    % Convert angles from degrees to radians for computation.
    incl_rad = deg2rad(incl);
    raan_rad = deg2rad(raan);
    argp_rad = deg2rad(argp);
    nu_rad   = deg2rad(nu);

    % Semi-latus rectum
    p = a * (1 - ecc^2);
    
    % Orbital radius at true anomaly
    r = p / (1 + ecc*cos(nu_rad));
    
    % Specific angular momentum
    h = sqrt(mu * p);
    
    % Position in perifocal (PQW) frame
    r_PQW = [ r*cos(nu_rad);
              r*sin(nu_rad);
              0 ];
    
    % Velocity in perifocal frame
    v_PQW = [ - (mu/h)*sin(nu_rad);
               (mu/h)*(ecc + cos(nu_rad));
               0 ];
    
    % Rotation matrices
    R3 = @(angle)[cos(angle), -sin(angle), 0; sin(angle), cos(angle), 0; 0, 0, 1];
    R1 = @(angle)[1, 0, 0; 0, cos(angle), -sin(angle); 0, sin(angle), cos(angle)];
    
    % Combined rotation matrix: from PQW to ECI
    % r_ECI = R_z(raan) * R_x(incl) * R_z(argp) * r_PQW
    R_ECI = R3(raan_rad) * R1(incl_rad) * R3(argp_rad);
    
    % Compute ECI position and velocity
    r_vec = R_ECI * r_PQW;
    v_vec = R_ECI * v_PQW;
end