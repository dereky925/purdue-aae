function v2 = herrickGibbs(r1, r2, r3, t1, t2, t3)
    % Herrick-Gibbs Method for Velocity Calculation at t2 (with mu in km^3/s^2)
    % Inputs:
    % r1, r2, r3 - Position vectors at times t1, t2, t3 (in kilometers)
    % t1, t2, t3 - Times corresponding to the position vectors (in seconds)
    
    % Gravitational constant of Earth in km^3/s^2
    mu = 3.986004418e5;
    
    % Compute time differences
    dt21 = t2 - t1;
    dt32 = t3 - t2;
    dt31 = t3 - t1;
    
    %Get the velocity with influenced taylor series equation
    v2 = (-dt32 .* (1./(dt21 .* dt31) + (mu ./ (12 .* norm(r1)^3))) .* r1)...
         + ((dt32 - dt21) .* (1/(dt21 .* dt32) + (mu ./ (12 * norm(r2)^3))) .* r2) ...
        + (dt21 .* (1./(dt32 .* dt31) + (mu ./ (12 .* norm(r3)^3))) .* r3);

end
