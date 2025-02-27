function a_RTN = eci2rtn(r_eci, v_eci, a_eci)
    
    % Inputs:
    %   r_eci - Position vector in ECI [km]
    %   v_eci - Velocity vector in ECI [km/s]
    %   a_eci - Acceleration vector in ECI [3x1]
    
    % Output:
    %   a_RTN - Vector in RTN [3x1]

    % Compute radial unit vector (R)
    R_hat = r_eci / norm(r_eci);

    % Compute angular momentum vector (h = r x v)
    h_eci = cross(r_eci, v_eci);
    N_hat = h_eci / norm(h_eci);  % Normal unit vector (N)

    % Compute transverse unit vector (T)
    T_hat = cross(N_hat, R_hat);

    % Assemble rotation matrix from ECI to RTN
    R_ECI_to_RTN = [R_hat'; T_hat'; N_hat'];

    % Convert acceleration to RTN frame
    a_RTN = R_ECI_to_RTN * a_eci;

end