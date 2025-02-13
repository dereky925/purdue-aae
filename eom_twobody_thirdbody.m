function dXdt = eom_twobody_thirdbody(t, Y, mu_E, mu_M, dEM, n)

    % Inputs:
    %   t    = time, sec
    %   Y    = [x, y, z, vx, vy, vz] w.r.t. Earth at origin
    %   muE  = Earth GM [km^3/s^2]
    %   muM  = Moon GM [km^3/s^2]
    %   dEM  = Earth-Moon distance [km] (circular orbit assumed)
    %   n    = Earth-Moon angular velocity [rad/s]

    % Output:
    %   dXdt = time derivative of state = [vx, vy, vz, ax, ay, az]
    
    % Unpack
    x = Y(1);  y = Y(2);  z = Y(3);
    vx = Y(4); vy = Y(5); vz = Y(6);
    
    rE  = [x; y; z]; % spacecraft pos w.r.t. Earth
    rE_mag = norm(rE); % magnitude
    
    % Two-body acceleration due to Earth
    aEarth = -mu_E * rE / (rE_mag^3);
    
    % Moon position (assuming circular orbit in xy-plane, centered at Earth)
    xM = dEM*cos(n*t);
    yM = dEM*sin(n*t);
    zM = 0;
    rM = [xM; yM; zM]; % Moon position w.r.t. Earth
    rM_mag = norm(rM);
    
    % Vector from spacecraft to Moon
    rmDiff = rE - rM;
    rmD = norm(rmDiff);
    
    % Third-body acceleration due to Moon
    aMoon = -mu_M * ( rmDiff/(rmD^3) + rM/(rM_mag^3) );
    
    % Total acceleration
    aTotal = aEarth + aMoon;
    
    dXdt = [vx; vy; vz; aTotal(1); aTotal(2); aTotal(3)];
end