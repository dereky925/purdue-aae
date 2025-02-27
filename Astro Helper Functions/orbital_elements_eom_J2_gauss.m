function dXdt = orbital_elements_eom_J2_gauss(~, X)

    %   dXdt = orbital_elements_eom_J2(t, X)

    %   Inputs:
    %     X = [a, e, i, Omega, omega, nu]
    %         where
    %           a     - Semi-major axis [km]
    %           e     - Eccentricity (dimensionless)
    %           i     - Inclination [deg]
    %           Omega - RAAN [deg]
    %           omega - Argument of perigee [deg]
    %           nu    - True anomaly [deg]
    %
    %   Output:
    %     dXdt - Time derivative of the orbital elements:
    %            dXdt(1) = da/dt [km/s]
    %            dXdt(2) = de/dt
    %            dXdt(3) = di/dt [deg/s]
    %            dXdt(4) = dOmega/dt [deg/s]
    %            dXdt(5) = domega/dt [deg/s]
    %            dXdt(6) = dnu/dt [deg/s]

    % Extract and convert state variables
    a     = X(1);
    e     = X(2);
    i     = deg2rad(X(3));
    Omega = deg2rad(X(4));
    omega = deg2rad(X(5));
    nu    = deg2rad(X(6));

    % Derived quantities
    p = a * (1 - e^2);         % Semi-latus rectum [km]
    r = p / (1 + e*cos(nu));     % Orbital radius [km]
    h = sqrt(MU * p);          % Specific angular momentum [km^2/s]

    % Convert orbital elements to ECI position and velocity
    [r_eci, v_eci] = keplerian2eci(a, e, rad2deg(i), rad2deg(Omega), rad2deg(omega), rad2deg(nu));
    
    % Compute the J2 perturbing acceleration in ECI
    r_norm = norm(r_eci);
    a_J2_ECI = -3*MU*J2*EARTH_RADIUS^2/(2*r_norm^5) * ...
        [ (1 - 5*(r_eci(3)/r_norm)^2)*r_eci(1);
          (1 - 5*(r_eci(3)/r_norm)^2)*r_eci(2);
          (3 - 5*(r_eci(3)/r_norm)^2)*r_eci(3) ];
    
    % Convert J2 acceleration to the RTN frame
    a_J2_RTN = eci2rtn(r_eci, v_eci, a_J2_ECI);
    a_R = a_J2_RTN(1);
    a_T = a_J2_RTN(2);
    a_N = a_J2_RTN(3);

    % Gauss Planetary Equations
    % Lecture 5 or 10
    da_dt = (2*a^2/h) * ( e*sin(nu)*a_R + (p/r)*a_T );
    de_dt = (1/h) * ( p*sin(nu)*a_R + ((p + r)*cos(nu) + r*e)*a_T );
    di_dt = (r*cos(nu+omega)/h) * a_N;
    dOmega_dt = (r*sin(nu+omega)/(h*sin(i))) * a_N;
    domega_dt = (1/(h*e)) * ( -p*cos(nu)*a_R + ((p + r)*sin(nu))*a_T ) ...
                - (r*sin(nu+omega)/(h*tan(i)))*a_N;
    dnu_dt = h/r^2 + (1/(h*e)) * ( p*cos(nu)*a_R - ((p + r)*sin(nu))*a_T );

    % Assemble output vector; convert angular rates to degrees per second.
    dXdt = zeros(6,1);
    dXdt(1) = da_dt;                  % [km/s]
    dXdt(2) = de_dt;                  
    dXdt(3) = rad2deg(di_dt);          % [deg/s]
    dXdt(4) = rad2deg(dOmega_dt);      % [deg/s]
    dXdt(5) = rad2deg(domega_dt);      % [deg/s]
    dXdt(6) = rad2deg(dnu_dt);         % [deg/s]

end

