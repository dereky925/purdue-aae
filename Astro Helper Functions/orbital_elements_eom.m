function dXdt = orbital_elements_eom(t, x, A_M, utcTime)

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
    % dXdt(6) = Mean anomaly [°]

    mu = 398600.4418; % [km^3/s^2]

    [r_vector,v_vector] = keplerian2eci(x(1),x(2),x(3),x(4),x(5),x(6));

    % Extract orbital elements
    a = x(1);                 % Semi-major axis [km]
    e = x(2);                 % Eccentricity
    i = deg2rad(x(3));        % Inclination [rad]
    raan = deg2rad(x(4));     % RAAN [rad]
    argp = deg2rad(x(5));    % Argument of perigee [rad]
    nu = deg2rad(x(6));       % True anomaly [rad]

    % Increment utcTime
    current_utcTime = utcTime + seconds(t);
    JD = juliandate(current_utcTime);
    
    % Compute J2 Effects ==================================================
    %https://www.vcalc.com/wiki/eng/J2
    J2 = 0.00108262668;
    
    RE = 6371; % Earth Radius [km]
    
    % Compute time in seconds since reference epoch
    time_since_epoch = JD * 24*60*60;
    
    % Calculate the Earth rotation angle since epoch
    omega = 7.2921159e-5; % Earth's angular velocity in rad/s, sidereal day
    theta = omega * time_since_epoch;
    
    % Rotate ECI position vector to ECEF
    r_vec_ECEF = R3(theta) * r_vector;
    
    % Compute ECEF position norm
    r_vec_ECEF_norm = norm(r_vec_ECEF);
    
    s = r_vec_ECEF(1)/r_vec_ECEF_norm;
    t = r_vec_ECEF(2)/r_vec_ECEF_norm;
    u = r_vec_ECEF(3)/r_vec_ECEF_norm;
    
    % Compute J2 effects at ECEF position
    u_J2 = 3/2 * [5*s*u^2 - s; 5*t*u^2 - t; 5*u^3 - 3*u];
    a_J2_ECEF = mu*J2/r_vec_ECEF_norm^2*(RE/r_vec_ECEF_norm)*u_J2; % J2
    
    % Rotate back to ECI
    a_J2_ECI = R3(-theta) * a_J2_ECEF;

    % Rotate to RTN frame (radial, transverse, normal)
    a_J2_RTN = eci2rtn(r_vector, v_vector, a_J2_ECI);

    % Derived quantities
    n = sqrt(mu / a^3);        % Mean motion [rad/s]
    b = a * sqrt(1 - e^2);     % Semi-minor axis [km]
    p = a * (1 - e^2);         % Semi-latus rectum [km]
    
    % Compute orbital radius
    r = p / (1 + e * cos(nu));

    % Transformation matrix B(x)
    B = (1 / sqrt(mu * p)) * ...
        [ 2*a^2*e*sin(nu),    2*a^2*p/r,      0;
          p*sin(nu), (p + r)*cos(nu) + r*e,   0;
          0, 0, r*cos(nu + argp);
          0, 0, r*sin(nu + argp)/sin(i);
         -p*cos(nu)/e, (p + r)*sin(nu)/e,   -r*sin(nu+argp)/tan(i);
          b*p*cos(nu)/(a*e)-2*b*r/a, -(b*(p + r)*sin(nu))/(a*e), 0];
    
    % Nominal Keplerian dynamics
    f0 = [0; 0; 0; 0; 0; n];

    % Solar Radiation Pressure ============================================
    G_0 = 1.02E14; % [kg*km/s^2] solar flux constant
    au = 1.496e11; % Astronomical Unit
    r_sun = sun(JD)'*au; % Distance to sun
    sun_unit = (x(1:3)-r_sun)/norm(x(1:3)-r_sun); % Sun unit vector
    a_SRP_ECI = -A_M * G_0 / (norm(x(1:3)-r_sun)^2) * sun_unit; % Cannon-ball model

    % Rotate SRP to RTN frame (radial, transverse, normal)
    a_SRP_RTN = eci2rtn(r_vector, v_vector, a_SRP_ECI);

    % Perturbed dynamics due to acceleration ad in ECI
    ad = a_J2_RTN + a_SRP_RTN;
    dXdt = f0 + B * ad;

    dXdt = dXdt(:);
    dXdt(3) = rad2deg(dXdt(3));
    dXdt(4) = rad2deg(dXdt(4));
    dXdt(5) = rad2deg(dXdt(5));
    dXdt(6) = rad2deg(dXdt(6));

end