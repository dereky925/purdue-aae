function [a, ex, ey, incl, raan, theta] = cartesian2elements(r_vec, v_vec)
% CARTESIAN2ELEMENTS Convert ECI state vectors to orbital elements.
%
%   [a, ex, ey, incl, raan, theta] = cartesian2elements(r_vec, v_vec)
%
%   Inputs:
%       r_vec  = ECI position vector [km] (3x1)
%       v_vec  = ECI velocity vector [km/s] (3x1)
%
%   Outputs:
%       a     = Semi-major axis [km]
%       ex    = x-component of eccentricity (unitless)
%       ey    = y-component of eccentricity (unitless)
%       incl  = Inclination [°]
%       raan  = Right ascension of ascending node [°]
%       theta = Argument of latitude [°]  (i.e. argp + true anomaly mod 360)
%
%   This function converts the ECI state (r_vec, v_vec) into classical Keplerian
%   elements, then extracts the x- and y-components of the eccentricity vector and
%   computes the argument of latitude.
%
%   Author: Derek Yu (adapted)
%   Date: Jan 21, 2025

    mu = 398600.4418;  % Earth's gravitational parameter [km^3/s^2]
    
    % Compute norms and angular momentum vector
    r = r_vec;
    v = v_vec;
    r_norm = norm(r);
    v_norm = norm(v);
    h_vec = cross(r, v);
    h = norm(h_vec);
    
    % Eccentricity vector and its magnitude
    ecc_vec = cross(v, h_vec)/mu - r/r_norm;
    ecc = norm(ecc_vec);
    
    % Inclination (in degrees)
    incl = acosd(h_vec(3)/h);
    
    % Compute the node vector (n = k x h)
    k_vec = [0; 0; 1];
    n_vec = cross(k_vec, h_vec);
    n_norm = norm(n_vec);
    
    % Right ascension of ascending node (RAAN)
    if n_norm < 1e-6
        raan = 0;  % Undefined for equatorial orbits; set to 0
    else
        if n_vec(2) >= 0
            raan = acosd(n_vec(1)/n_norm);
        else
            raan = 360 - acosd(n_vec(1)/n_norm);
        end
    end
    
    % True anomaly (nu)
    if ecc > 1e-6
        % Use the dot product to determine quadrant.
        if dot(r, v) >= 0
            nu = acosd(dot(ecc_vec, r)/(ecc*r_norm));
        else
            nu = 360 - acosd(dot(ecc_vec, r)/(ecc*r_norm));
        end
    else
        nu = 0;  % For circular orbits, true anomaly is undefined; set to 0
    end
    
    % Argument of perigee (argp)
    if n_norm < 1e-6 || ecc < 1e-6
        argp = 0; % For circular or equatorial orbits, undefined; set to 0
    else
        if ecc_vec(3) >= 0
            argp = acosd(dot(n_vec, ecc_vec)/(n_norm*ecc));
        else
            argp = 360 - acosd(dot(n_vec, ecc_vec)/(n_norm*ecc));
        end
    end
    
    % Semi-major axis (a)
    a = 1/(2/r_norm - v_norm^2/mu);
    
    % Compute the argument of latitude (theta)
    theta = mod(argp + nu, 360);
    
    % Extract eccentricity vector components in ECI x and y directions
    ex = ecc_vec(1);
    ey = ecc_vec(2);
end