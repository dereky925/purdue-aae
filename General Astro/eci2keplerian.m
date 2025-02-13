function [a, ecc, incl, raan, argp, nu, ecc_x, ecc_y, arglat] = eci2keplerian(r_vec, rdot_vec)
    % ECI2KEPLERIAN Convert ECI state vectors to Keplerian orbit elements.
    %
    %   [a, ecc, incl, raan, argp, nu] = eci2keplerian(r_vec, rdot_vec)
    %   returns the standard orbital elements:
    %       a     = Semi-major axis [km]
    %       ecc   = Eccentricity
    %       incl  = Inclination [°]
    %       raan  = Right ascension of ascending node [°]
    %       argp  = Argument of perigee [°]
    %       nu    = True anomaly [°]
    %
    %   [a, ecc, incl, raan, argp, nu, ecc_x, ecc_y, arglat] = eci2keplerian(r_vec, rdot_vec)
    %   additionally returns:
    %       ecc_x = The x-component of the eccentricity vector (ECI) [unitless]
    %       ecc_y = The y-component of the eccentricity vector (ECI) [unitless]
    %       arglat = Argument of latitude, defined as (argp + nu) mod 360 [°]
    %
    %   Inputs:
    %       r_vec    = ECI position vector [km]
    %       rdot_vec = ECI velocity vector [km/s]
    %
    %   Reference: https://downloads.rene-schwarz.com/download/M002-Cartesian_State_Vectors_to_Keplerian_Orbit_Elements.pdf
      
    
    % Compute angular momentum vector
    h_vec = cross(r_vec, rdot_vec);
    
    % Compute eccentricity vector
    ecc_vec = cross(rdot_vec, h_vec) / MU - r_vec / norm(r_vec);
    
    % Compute n vector pointing to ascending node
    n_vec = cross([0 0 1]', h_vec);
    
    % Compute the dot product value
    val = dot(ecc_vec, r_vec) / (norm(ecc_vec) * norm(r_vec));
    % Clamp it to the valid domain of acosd:
    val = min(max(val, -1), 1);
    if dot(r_vec, rdot_vec) >= 0
        nu = acosd(val);
    else
        nu = 360 - acosd(val);
    end
    
    % Compute inclination (degrees)
    incl = acosd(h_vec(3) / norm(h_vec));
    
    % Compute eccentricity magnitude
    ecc = norm(ecc_vec);
    
    % Compute right ascension of ascending node (RAAN)
    if n_vec(2) >= 0
        raan = acosd(n_vec(1) / norm(n_vec));
    else
        raan = 360 - acosd(n_vec(1) / norm(n_vec));
    end
    
    % Compute argument of perigee
    if ecc_vec(3) >= 0
        argp = acosd(dot(n_vec, ecc_vec) / (norm(n_vec) * norm(ecc_vec)));
    else
        argp = 360 - acosd(dot(n_vec, ecc_vec) / (norm(n_vec) * norm(ecc_vec)));
    end
    
    % Compute semi-major axis
    a = 1 / (2 / norm(r_vec) - norm(rdot_vec)^2 / MU);
    
    % Optionally output x and y components of the eccentricity vector
    if nargout > 6
        ecc_x = ecc_vec(1);
        ecc_y = ecc_vec(2);
    end
    
    % Optionally output the argument of latitude, u = (argp + nu) mod 360
    if nargout > 8
        arglat = mod(argp + nu, 360);
    end

end