function dXdt = stm_2B_J2(~, X)
    % Unpack state and STM from the augmented state vector X.
    pos = X(1:3);           % position vector [km]
    vel = X(4:6);           % velocity vector [km/s]
    Phi = reshape(X(7:end), 6, 6);  % state transition matrix
    
    % Compute the norm and its square
    r2 = pos' * pos;
    r  = sqrt(r2);

    RE = EARTH_RADIUS;
    
    % Two-body acceleration
    a_grav = -MU * pos / r^3;

    % J2 perturbation acceleration
    C = 1.5 * J2 * MU * RE^2;
    fact = C / r^5;
    x = pos(1); y = pos(2); z = pos(3);
    a_J2 = fact * [ x*(5*(z^2/r2) - 1);
                    y*(5*(z^2/r2) - 1);
                    z*(5*(z^2/r2) - 3) ];

    % Total acceleration
    acc = a_grav + a_J2;

    A  = jacobian2BJ2([pos vel]);
    
    % Propagate the STM
    dPhi = A * Phi;
    
    dXdt = [vel;    % derivative of position is velocity
            acc;    % derivative of velocity is acceleration
            dPhi(:)];% STM derivative (flattened into a 36x1 column vector)

    
end