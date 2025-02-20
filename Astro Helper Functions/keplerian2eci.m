function [r_vec, rdot_vec] = keplerian2eci(varargin)

% Usage 1 (traditional):
%   [r_vec, rdot_vec] = keplerian2eci(a, ecc, incl, raan, argp, nu)
%
%   Inputs:
%       a     = Semi-major axis [km]
%       ecc   = Eccentricity (scalar)
%       incl  = Inclination [°]
%       raan  = Right ascension of ascending node [°]
%       argp  = Argument of perigee [°]
%       nu    = True anomaly [°]
%
% Usage 2 (alternative with separate eccentricity components):
%   [r_vec, rdot_vec] = keplerian2eci(a, ecc_x, ecc_y, incl, raan, arglat, 'alt')
%
%   Inputs:
%       a       = Semi-major axis [km]
%       ecc_x   = x-component of eccentricity (scalar)
%       ecc_y   = y-component of eccentricity (scalar)
%       incl    = Inclination [°]
%       raan    = Right ascension of ascending node [°]
%       arglat  = Argument of latitude [°] (i.e. argp + nu mod 360)
%       'alt'   = A flag string indicating alternative calling convention.
%
%   Outputs (both cases):
%       r_vec    = ECI position vector [km] (3x1)
%       rdot_vec = ECI velocity vector [km/s] (3x1)
%
% Reference: See
% https://downloads.rene-schwarz.com/download/M002-Cartesian_State_Vectors_to_Keplerian_Orbit_Elements.pdf

    MU = 398600.4418;  % [km^3/s^2]

    % Check input count
    if nargin == 6
        % Traditional call: a, ecc, incl, raan, argp, nu
        a     = varargin{1};
        ecc   = varargin{2};
        incl  = varargin{3};
        raan  = varargin{4};
        argp  = varargin{5};
        nu    = varargin{6};
    elseif nargin == 7 && ischar(varargin{7}) && strcmp(varargin{7}, 'alt')
        % Alternative call: a, ecc_x, ecc_y, incl, raan, arglat, 'alt'
        a       = varargin{1};
        ecc_x   = varargin{2};
        ecc_y   = varargin{3};
        incl    = varargin{4};
        raan    = varargin{5};
        arglat  = varargin{6};
        % Compute scalar eccentricity and argument of perigee.
        ecc = sqrt(ecc_x^2 + ecc_y^2);
        argp = rad2deg(atan2(ecc_y, ecc_x));
        % Compute true anomaly as difference between argument of latitude and argp.
        nu = mod(arglat - argp, 360);
    else
        error('Incorrect number of input arguments. Use either 6 inputs (traditional) or 7 inputs with the flag ''alt''.');
    end

    % Convert angles to radians.
    incl_rad = deg2rad(incl);
    raan_rad = deg2rad(raan);
    argp_rad = deg2rad(argp);
    nu_rad   = deg2rad(nu);
    
    % Compute semi-latus rectum
    p = a * (1 - ecc^2);
    
    % Compute magnitude of orbital radius at the given true anomaly.
    r = p / (1 + ecc*cos(nu_rad));
    
    % Compute specific angular momentum
    h = sqrt(MU * p);
    
    % Position in the perifocal (PQW) coordinate system.
    r_PQW = [ r*cos(nu_rad);
              r*sin(nu_rad);
              0 ];
    
    % Velocity in the perifocal frame.
    v_PQW = [ - (MU/h)*sin(nu_rad);
               (MU/h)*(ecc + cos(nu_rad));
               0 ];
    
    % Combined rotation matrix from PQW to ECI:
    % r_ECI = R_z(raan) * R_x(incl) * R_z(argp) * r_PQW
    ECI_ROT = R3(-raan_rad) * R1(-incl_rad) * R3(-argp_rad);
    
    % Compute ECI position and velocity vectors.
    r_vec = ECI_ROT * r_PQW;
    rdot_vec = ECI_ROT * v_PQW;
    
end