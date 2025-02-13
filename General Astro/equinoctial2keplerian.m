function [a, e, i, raan, argp, nu] = equinoctial2keplerian(p,f,g,h,k,L) 

    % Inputs
    % p = Semi-parameter [km]
    % f = Equinoctial component of eccentricity
    % g = Equinoctial component of eccentricity
    % h = Equinoctial component of inclination
    % k = Equinoctial component of inclination
    % L = True longitude [rad]

    % Outputs
    % a = Semi-major axis [km]
    % e = Eccentricity
    % i = Inclination [°]
    % raan = RAAN [°]
    % argp = Argument of perigee [°]
    % nu = Mean anomaly [°]

    a = p/(1-f^2-g^2);  % Semi-major axis [km]
    e = sqrt(f^2 + g^2);  % Eccentricity
    i = atan2d(2*sqrt(h^2+k^2),1-h^2-k^2); % Inclination [°]
    raan = atan2d(k,h);      % RAAN [°]
    argp = atan2d(g,f) - atan2d(k,h);      %  Argument of perigee [°]
    nu = rad2deg(L - (atan2(k,h) + ...
                atan2(g*h-f*k,...
                f*h+g*k))); % True anomaly [°]


end
