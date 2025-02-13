function dXdt = milankovitch_elements_eom(t, x, A_M, utcTime)

    % Inputs
    % x - Vector of Milankovitch elements
    % ijk - Vector of cartesian position and velocity

    % Outputs
    % dXdt - Propagated Milankovitch elements

    mu = 398600.4418; % [km^3/s^2]

    % Unpack Milankovitch elements
    h_vec = [x(1) x(2) x(3)]';
    e_vec = [x(4) x(5) x(6)]';
    L = x(7);

    [r_vector, v_vector] = milankovitch2eci(h_vec, e_vec, L);
    
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

    r_cross = cross_matrix(r_vector);
    v_cross = cross_matrix(v_vector);
    h_cross = cross_matrix(h_vec);

    B1 = r_cross;
    B2 = 1/mu * (v_cross*r_cross - h_cross);

    a = r_vector(3) * h_vec';
    b = norm(h_vec) * (norm(h_vec) + h_vec(3));
    B3 = a / b;

    % Transformation matrix B(x)
    B = [B1; B2; B3];
    

    f0 = [0; 0; 0; 0; 0; 0; norm(h_vec)/norm(r_vector)^2];

    % Solar Radiation Pressure ============================================
    G_0 = 1.02E14; % [kg*km/s^2] solar flux constant
    au = 1.496e11; % Astronomical Unit
    r_sun = sun(JD)'*au; % Distance to sun
    sun_unit = (x(1:3)-r_sun)/norm(x(1:3)-r_sun); % Sun unit vector
    a_SRP_ECI = -A_M * G_0 / (norm(x(1:3)-r_sun)^2) * sun_unit; % Cannon-ball model

    % Rotate SRP to RTN frame (radial, transverse, normal)
    a_SRP_RTN = eci2rtn(r_vector, v_vector, a_SRP_ECI);

    ad = a_J2_RTN + a_SRP_RTN;

    % Perturbed dynamics due to acceleration ad in ECI
    dXdt = f0 + B * ad;

end