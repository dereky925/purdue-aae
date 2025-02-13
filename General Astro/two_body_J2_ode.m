function xdot = two_body_J2_ode(t,x,utcTime)

    r_vector = x(1:3,1); % 3x1 vector
    vel_vector = x(4:6,1); % 3X1 vector
    r = sqrt(r_vector(1)^2 + r_vector(2)^2 + r_vector(3)^2);
    
    % % Increment utcTime
    % current_utcTime = utcTime + seconds(t);
    % JD = juliandate(current_utcTime);
    % 
    % % Compute time in seconds since reference epoch
    % time_since_epoch = JD * 24*60*60;
    % 
    % % Calculate the Earth rotation angle since epoch
    % omega = 7.2921159e-5; % Earth's angular velocity in rad/s, sidereal day
    % theta = omega * time_since_epoch;
    % 
    % % Rotate ECI position vector to ECEF
    % r_vec_ECEF = R3(theta) * r_vector;
    % 
    % % Compute ECEF position norm
    % r_vec_ECEF_norm = norm(r_vec_ECEF);
    % 
    % s_comp = r_vec_ECEF(1)/r_vec_ECEF_norm;
    % t_comp = r_vec_ECEF(2)/r_vec_ECEF_norm;
    % u_comp = r_vec_ECEF(3)/r_vec_ECEF_norm;
    % 
    % % Compute J2 effects at ECEF position
    % u_J2 = 3/2 * [5*s_comp*u_comp^2 - s_comp; 5*t_comp*u_comp^2 - t_comp; 5*u_comp^3 - 3*u_comp];
    % a_J2_ECEF = MU*J2/r_vec_ECEF_norm^2*(EARTH_RADIUS/r_vec_ECEF_norm)*u_J2; % J2
    % 
    % 
    % % Rotate back to ECI
    % a_J2_ECI = R3(-theta) * a_J2_ECEF;

    % ECI implementation
    a_J2_ECI = -3*MU*J2*EARTH_RADIUS^2/(2*r^5) * ...
                   [(1-5*x(3)^2/r^2)*x(1) ;...
                   (1-5*x(3)^2/r^2)*x(2) ; ...
                   (3-5*x(3)^2/r^2)*x(3)];
    
    % Compute two-body acceleration ===========================================
    two_body_accel = -MU/r^3 * r_vector; % 3x1 vector
    
    rddot =  two_body_accel + a_J2_ECI;
    xdot = [vel_vector;rddot];

end