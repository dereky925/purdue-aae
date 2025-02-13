function dXdt = cr3bp_EOM(~, X, mu)

    %  The state X = [x; y; z; vx; vy; vz] is in nondimensional units, where:
    %  - Lengths are scaled by the Earth-Moon distance.
    %  - Time is scaled so that one nondimensional unit is one "Earth-Moon orbital period / 2π".
    %  - mu is the mass ratio: mu = M_moon / (M_earth + M_moon)

    % Unpack the state
    x  = X(1);
    y  = X(2);
    z  = X(3);
    vx = X(4);
    vy = X(5);
    vz = X(6);

    % Reference: Orbital Mechanics for engineering students pg 89-90
    % Note: these are non-dimensionalized

    % Distances from spacecraft to Earth
    r1 = sqrt( (x + mu)^2 + y^2 + z^2 ); 

    % Distances from spacecraft to Moon
    r2 = sqrt( (x - (1 - mu))^2 + y^2 + z^2 );

    % CR3BP accelerations
    ax =  2*vy + x - (1 - mu)/r1^3*(x + mu) - mu/r2^3*(x - (1 - mu));
    ay = -2*vx + y - (1 - mu)/r1^3*y - mu/r2^3*y;
    az = -(1 - mu)*z/r1^3 - mu*z/r2^3;

    dXdt = [vx; vy; vz; ax; ay; az];
    
end