function [p,f,g,h,k,L] = keplerian2equinoctial(a, e, i, raan, argp, nu)

    % Inputs
    % p = Semi-parameter [km]
    % f = Equinoctial component of eccentricity
    % g = Equinoctial component of eccentricity
    % h = Equinoctial component of inclination
    % k = Equinoctial component of inclination
    % L = True longitude [°]

    % Outputs
    % p = Semi-parameter [km]
    % f = Equinoctial component of eccentricity
    % g = Equinoctial component of eccentricity
    % h = Equinoctial component of inclination
    % k = Equinoctial component of inclination
    % L = True longitude [°]

    % Convert orbital elements into radians
    i = deg2rad(i);        % Inclination [rad]
    raan = deg2rad(raan);     % RAAN [rad]
    argp = deg2rad(argp);     % Argument of perigee [rad]
    nu = deg2rad(nu);       % True anomaly [rad]

    % Compute modified equinoctial elements
    p = a*(1-e^2);          % p = Semi-parameter [km]
    f = e*cos(argp+raan);   % f = Equinoctial component of eccentricity
    g = e*sin(argp+raan);   % g = Equinoctial component of eccentricity
    h = tan(i/2)*cos(raan); % h = Equinoctial component of inclination
    k = tan(i/2)*sin(raan); % k = Equinoctial component of inclination
    L = raan + argp + nu;   % L = True longitude [rad]

end