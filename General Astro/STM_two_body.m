function xdot = STM_two_body(t,states)

    % Unpack state
    x = states(1);
    y = states(2);
    z = states(3);

    mu = MU;
    r_0 = EARTH_RADIUS;
    
    % Define r (magnitude of position vector)
    r = sqrt(x^2 + y^2 + z^2);
    r_vec = [x y z]';
    
    % 2 Body accel jacobian
    jacob_a2B = -mu*(eye(3)/r^3 + r_vec*(-3/r^5*r_vec'));
    
    % Form matlab computed A matrix
    A = [zeros(3) eye(3); jacob_a2B zeros(3)];

    % Extract STM from state vector
    Phi = reshape(X(7:end), 6, 6);

    % Compute STM derivative
    Phi_dot = A * Phi;

    % Construct the derivative of the full state vector
    dX_dt = [vx; vy; vz; accel; Phi_dot(:)];


end