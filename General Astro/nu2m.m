function [M,E] = nu2m(nu, ecc)

    % Derek Yu - Jan 25, 2025
    % Reference: 
    % D. A. Vallado, Fundamentals of Astrodynamics and Applications, 4th ed.

    % Convert true anomaly to mean anomaly
    % Inputs:
    % nu: True anomaly [rad]
    % ecc: Eccentricity

    % Outputs:
    % M: Mean anomaly [rad]
    % E: Eccentric Anomaly [rad]

    E = atan2( sin(nu)*sqrt(1-ecc^2) , (ecc+cos(nu)) );

    M = E - ecc*sin(E);

    % Ensure within 360
    M = mod(M, 2*pi);

    % Ensure within 360
    E = mod(E, 2*pi);
    
end