function dXdt = modified_equinoctial_elements_eom2(t, x, utcTime)

    % Propagate modified equinoctial elements. These elements have less singularties
    % Reference: https://spsweb.fltops.jpl.nasa.gov/portaldataops...
    % /mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf
    
    % Inputs:
    % t  - Time (not used explicitly, required by ODE solvers)
    % X  - State vector of orbital elements [a, e, i, Omega, omega, nu]
    %      in degrees
    
    % Outputs:
    % dXdt(1) = Semi-major axis [km]
    % dXdt(2) = Eccentricity
    % dXdt(3) = Inclination [°]
    % dXdt(4) = RAAN [°]
    % dXdt(5) = Argument of perigee [°]
    % dXdt(6) = Mean anomaly [°]


    mu = 398600.4418; % [km^3/s^2]

    [r_vector,v_vector] = keplerian2eci(x(1),x(2),x(3),x(4),x(5),x(6));

    % Extract orbital elements
    a = x(1);                 % Semi-major axis [km]
    e = x(2);                 % Eccentricity
    i = deg2rad(x(3));        % Inclination [rad]
    raan = deg2rad(x(4));     % RAAN [rad]
    argp = deg2rad(x(5));     % Argument of perigee [rad]
    nu = deg2rad(x(6));       % True anomaly [rad]

    % Compute modified equinoctial elements
    p = a*(1-e^2);          % p = Semi-parameter [km]
    f = e*cos(argp+raan);   % f = Equinoctial component of eccentricity
    g = e*sin(argp+raan);   % g = Equinoctial component of eccentricity
    h = tan(i/2)*cos(raan); % h = Equinoctial component of inclination
    k = tan(i/2)*sin(raan); % k = Equinoctial component of inclination
    L = raan + argp + nu;   % L = True longitude [rad]
    % L = mod(L,2*pi);

    % Increment utcTime
    current_utcTime = utcTime + seconds(t);
    JD = juliandate(current_utcTime);
    
    % Compute J2 Effects ==================================================
    %https://www.vcalc.com/wiki/eng/J2
    J2 = 0.00108262668;
    
    RE = 6371; % Earth Radius [km]
    
    % Compute time in seconds since reference epoch
    time_since_epoch = JD * 24*60*60;
    
    % Calculate the Earth rotation angle since epoch
    omega = 7.2921159e-5; % Earth's angular velocity in rad/s, sidereal day
    theta = omega * time_since_epoch;
    
    % Rotate ECI position vector to ECEF
    r_vec_ECEF = R3(theta) * r_vector;
    
    % Compute ECEF position norm
    r_vec_ECEF_norm = norm(r_vec_ECEF);
    
    s = r_vec_ECEF(1)/r_vec_ECEF_norm;
    t = r_vec_ECEF(2)/r_vec_ECEF_norm;
    u = r_vec_ECEF(3)/r_vec_ECEF_norm;
    
    % Compute J2 effects at ECEF position
    u_J2 = 3/2 * [5*s*u^2 - s; 5*t*u^2 - t; 5*u^3 - 3*u];
    a_J2_ECEF = mu*J2/r_vec_ECEF_norm^2*(RE/r_vec_ECEF_norm)*u_J2; % J2
    
    % Rotate back to ECI
    a_J2_ECI = R3(-theta) * a_J2_ECEF;

    % Rotate to RTN frame (radial, transverse, normal)
    a_RTN = eci2rtn(r_vector, v_vector, a_J2_ECI);

    ad = a_RTN;

    % Derived quantities
    q = 1 + f*cos(L) + g*sin(L);
    s = sqrt(1+h^2+k^2);
   
    % Transformation matrix B(x)
    B = sqrt(p/mu) * ...
        [0 2*p/q 0;
         sin(L)   ((q+1)*cos(L)+f)/q   -(g*(h*sin(L)-k*cos(L)))/q;
         -cos(L)  ((q+1)*sin(L)+g)/q      f*(h*sin(L)-k*cos(L))/q;
         0 0 (s^2)/(2*q)*cos(L);
         0 0 (s^2)/(2*q)*sin(L);
         0 0 (h*sin(L)-k*cos(L))/q];
    
    % Nominal Keplerian dynamics
    f0 = [0; 0; 0; 0; 0; sqrt(mu*p)*(q/p)^2];

    % Perturbed dynamics due to acceleration ad in ECI
    % mod_eq_el = f0 + B * ad;
    mod_eq_el = f0;
    p_out = mod_eq_el(1);
    f_out = mod_eq_el(2);
    g_out = mod_eq_el(3);
    h_out = mod_eq_el(4);
    k_out = mod_eq_el(5);
    L_out = mod_eq_el(6);

    % Convert back into orbital elements for output
    dXdt(1) = p_out/(1-f_out^2-g_out^2);  % Semi-major axis [km]
    dXdt(2) = sqrt(f_out^2 + g_out^2);  % Eccentricity
    dXdt(3) = atan2d(2*sqrt(h_out^2+k_out^2),1-h_out^2-k_out^2); % Inclination [°]
    dXdt(4) = atan2d(k_out,h_out);      % RAAN [°]
    dXdt(5) = atan2d(g_out,f_out) - atan2d(k_out,h_out);      %  Argument of perigee [°]
    dXdt(6) = rad2deg(L_out - (atan2(k_out,h_out) + ...
                atan2(g_out*h_out-f_out*k_out,...
                f_out*h_out+g_out*k_out))); % True anomaly [°]

    dXdt = dXdt';
    
end
