function xdot = two_body_J2_SRP_ode(t,x,A_M,utcTime)

    % Input 
    % A_M = [km^2/kg] Spacecraft ballistic coefficient
    
    mu = 398600.4418; % [km^3/s^2]
    r_vector = x(1:3,1); % 3x1 vector
    vel_vector = x(4:6,1); % 3X1 vector
    r = sqrt(r_vector(1)^2 + r_vector(2)^2 + r_vector(3)^2);
    
    % Increment utcTime
    current_utcTime = utcTime + seconds(t);
    JD = juliandate(current_utcTime);
    
    % Compute J2 Effects ======================================================
    %https://www.vcalc.com/wiki/eng/J2
    J2 = 0.00108262668;
    
    RE = EARTH_RADIUS; % Earth Radius [km]
    
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
    
    % Solar Radiation Pressure ============================================
    G_0 = 1.02E14; % [kg*km/s^2] solar flux constant
    au = 1.496e11; % Astronomical Unit
    r_sun = sun(JD)'*au; % Distance to sun
    sun_unit = (x(1:3)-r_sun)/norm(x(1:3)-r_sun); % Sun unit vector
    a_SRP = -A_M * G_0 / (norm(x(1:3)-r_sun)^2) * sun_unit; % Cannon-ball model
    
    % Compute two-body acceleration =======================================
    two_body_accel = -mu/r^3 * r_vector; % 3x1 vector
    
    rddot =  two_body_accel + a_J2_ECI + a_SRP;
    xdot = [vel_vector;rddot];

end