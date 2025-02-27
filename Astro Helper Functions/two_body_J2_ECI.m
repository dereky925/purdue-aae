function xdot = two_body_J2_ECI(~,x)

    r_vector = x(1:3,1); % 3x1 vector
    vel_vector = x(4:6,1); % 3X1 vector
    r = sqrt(r_vector(1)^2 + r_vector(2)^2 + r_vector(3)^2);

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