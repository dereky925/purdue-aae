function [nu,E] = m2nu(M, ecc)

    % Derek Yu - Jan 25, 2025
    % Reference: 
    % D. A. Vallado, Fundamentals of Astrodynamics and Applications, 4th ed.

    % Inputs
    % M: Mean anomaly [Radians]
    % ecc: Eccentricity

    % Outputs
    % nu: True anomaly [rad]
    % E: Eccentric anomaly [rad]

    M = mod(M, 2*pi);
    if ecc < 0.1
        E = M; % near-circular guess
    else
        E = pi; %for higher eccentricity
    end

    maxIter = 50;
    tol     = 1e-12;

    for k = 1:maxIter
        
        f  = M - E + ecc*sin(E);
        fp = 1 - ecc*cos(E);
        dE = f/fp;
        E  = E + dE;

        if abs(dE) < tol
            break;
        end

    end

    % Ensure within 360
    E = mod(E, 2*pi);

    % Compute true anomaly
    x = cos(E) - ecc;
    y = sin(E)*sqrt(1-ecc^2);
    nu = (atan2(y,x));

end
