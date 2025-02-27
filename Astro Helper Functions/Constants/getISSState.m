function [r_eci, v_eci, a, e, i, RAAN, w, nu] = getISSState()

    % Random ISS states
    % r_eci = [656.0758  6596.8028 1473.0872]'; % [km]
    % v_eci = [-4.9631 -0.8127 5.7866]'; % [km/s]

    % [a, e, i, RAAN, w, nu] = eci2keplerian(r_eci,v_eci);

    a     = EARTH_RADIUS + 400;  % Semi-major axis [km]
    e     = 0.001;       % Eccentricity
    i     = 51.64;       % Inclination [deg]
    RAAN  = 180;           % RAAN [deg]
    w     = 0;           % Argument of perigee [deg]
    nu    = 0;

    [r_eci, v_eci] = keplerian2eci(a, e, i, RAAN, w, nu);


end