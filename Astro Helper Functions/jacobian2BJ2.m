function A_2B_J2 = jacobian2BJ2(states)

    x = states(1);
    y = states(2);
    z = states(3);

    mu = MU;
    J2 = 0.00108262668;
    r_0 = EARTH_RADIUS;

    r = sqrt(x^2 + y^2 + z^2);
    r_vec = [x y z]';
    
    % 2 Body accel jacobian
    jacob_a2B = -mu*(eye(3)/r^3 + r_vec*(-3/r^5*r_vec'));
    
    % jacobian aJ2 x y term
    A_xy = 1/r^5 - 5*z^2/r^7;
    jacob_aJ2_XY = r_vec *[ (35*z^2*x/r^9-5*x/r^7) ...
                               (35*z^2*y/r^9-5*y/r^7)...
                               (35*z^3/r^9 - 15*z/r^7)  ] + A_xy*eye(3);
    
    % jacobian aJ2 z term
    A_z = 3/r^5 - 5*z^2/r^7;
    jacob_aJ2_Z = r_vec *[ (35*z^2*x/r^9-15*x/r^7) ...
                               (35*z^2*y/r^9-15*y/r^7)...
                               (35*z^3/r^9 - 25*z/r^7)  ] + A_z*eye(3);
    
    % Create combined a_J2 jacobian
    jacobian_aJ2 = -3*mu*J2*r_0^2/2 * ([jacob_aJ2_XY(1:2,:)
                                       + jacob_aJ2_Z(3,:)]);
    
    % Form computed A matrix
    A_2B_J2 = [zeros(3), eye(3); (jacob_a2B + jacobian_aJ2) , zeros(3)];


end