function dXdt = orbital_elements_eom_J2_ECI(~, x)

    % Gauss planetary equations

    % Inputs:
    % t  - Time (not used explicitly, required by ODE solvers)
    % X  - State vector of orbital elements [a, e, i, Omega, omega, nu]
    %      in degrees
    % A_M = [km^2/kg] Spacecraft ballistic coefficient
    
    % Outputs:
    % dXdt(1) = Semi-major axis [km]
    % dXdt(2) = Eccentricity
    % dXdt(3) = Inclination [°]
    % dXdt(4) = RAAN [°]
    % dXdt(5) = Argument of perigee [°]
    % dXdt(6) = True anomaly [°]

    [r_vector,v_vector] = keplerian2eci(x(1),x(2),x(3),x(4),x(5),x(6));

    r = sqrt(r_vector(1)^2 + r_vector(2)^2 + r_vector(3)^2);

    % Extract orbital elements
    a = x(1);                 % Semi-major axis [km]
    e = x(2);                 % Eccentricity
    i = deg2rad(x(3));        % Inclination [rad]
    raan = deg2rad(x(4));     % RAAN [rad]
    argp = deg2rad(x(5));    % Argument of perigee [rad]
    nu = deg2rad(x(6));       % True anomaly [rad]

    % ECI J2 perturbation implementation
    a_J2_ECI = -3*MU*J2*EARTH_RADIUS^2/(2*r^5) * ...
                   [(1-5*r_vector(3)^2/r^2)*r_vector(1) ;...
                    (1-5*r_vector(3)^2/r^2)*r_vector(2) ; ...
                    (3-5*r_vector(3)^2/r^2)*r_vector(3)];

    % Rotate J2 perturbation to RTN frame (radial, transverse, normal)
    a_J2_RTN = eci2rtn(r_vector, v_vector, a_J2_ECI);

    % Derived quantities
    n = sqrt(MU / a^3);        % Mean motion [rad/s]
    b = a * sqrt(1 - e^2);     % Semi-minor axis [km]
    p = a * (1 - e^2);         % Semi-latus rectum [km]
    
    % Compute orbital radius
    r = p / (1 + e * cos(nu));

    % Transformation matrix B(x)
    B = (1 / sqrt(MU * p)) * ...
        [ 2*a^2*e*sin(nu),    2*a^2*p/r,      0;
          p*sin(nu), (p + r)*cos(nu) + r*e,   0;
          0, 0, r*cos(nu + argp);
          0, 0, r*sin(nu + argp)/sin(i);
         -p*cos(nu)/e, (p + r)*sin(nu)/e,   -r*sin(nu+argp)/tan(i);
          b*p*cos(nu)/(a*e)-2*b*r/a, -(b*(p + r)*sin(nu))/(a*e), 0];
    
    % Nominal Keplerian dynamics
    f0 = [0; 0; 0; 0; 0; n];

    % Perturbed dynamics due to acceleration ad in ECI
    ad = a_J2_RTN;
    dXdt = f0 + B * ad;

    dXdt = dXdt(:);
    dXdt(3) = rad2deg(dXdt(3));
    dXdt(4) = rad2deg(dXdt(4));
    dXdt(5) = rad2deg(dXdt(5));

    % Convert M to 
    conversion_factor = ((1 + e*cos(nu))^2) / ((1 - e^2)^(3/2));
    dXdt(6) = dXdt(6) * conversion_factor;
    dXdt(6) = rad2deg(dXdt(6));

    % Outputs:
    % dXdt(1) = Semi-major axis [km]
    % dXdt(2) = Eccentricity
    % dXdt(3) = Inclination [°]
    % dXdt(4) = RAAN [°]
    % dXdt(5) = Argument of perigee [°]
    % dXdt(6) = True anomaly [°]

end
