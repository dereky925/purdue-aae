function dXdt = modified_equinoctial_elements_eom(t, x, A_M, utcTime)

    % Propagate modified equinoctial elements. These elements have less singularties
    % Reference: https://spsweb.fltops.jpl.nasa.gov/portaldataops...
    % /mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf
    
    % Inputs:
    % t  - Time (not used explicitly, required by ODE solvers)
    % x  - State vector of modified equinoctial elements
    % x(1) = p = Semi-parameter [km]
    % x(2) = f = Equinoctial component of eccentricity
    % x(3) = g = Equinoctial component of eccentricity
    % x(4) = h = Equinoctial component of inclination
    % x(5) = k = Equinoctial component of inclination
    % x(6) = L = True longitude [rad]
    % A_M = [km^2/kg] Spacecraft ballistic coefficient
    
    % Outputs:
    % dXdt(1) = Semi-major axis [km]
    % dXdt(2) = Eccentricity
    % dXdt(3) = Inclination [°]
    % dXdt(4) = RAAN [°]
    % dXdt(5) = Argument of perigee [°]
    % dXdt(6) = Mean anomaly [°]

    mu = 398600.4418; % [km^3/s^2]

    % Store modified equinoctial elements
    p = x(1);         
    f = x(2);   
    g = x(3);    
    h = x(4);  
    k = x(5);   
    L = x(6);   

    [a, e, i, raan, argp, nu] = equinoctial2keplerian(p,f,g,h,k,L);
    [r_vector,v_vector] = keplerian2eci(a, e, i, raan, argp, nu);

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
    a_J2_RTN = eci2rtn(r_vector, v_vector, a_J2_ECI);

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

    % Solar Radiation Pressure ============================================
    G_0 = 1.02E14; % [kg*km/s^2] solar flux constant
    au = 1.496e11; % Astronomical Unit
    r_sun = sun(JD)'*au; % Distance to sun
    sun_unit = (x(1:3)-r_sun)/norm(x(1:3)-r_sun); % Sun unit vector
    a_SRP_ECI = -A_M * G_0 / (norm(x(1:3)-r_sun)^2) * sun_unit; % Cannon-ball model

    % Rotate SRP to RTN frame (radial, transverse, normal)
    a_SRP_RTN = eci2rtn(r_vector, v_vector, a_SRP_ECI);

    ad = a_J2_RTN + a_SRP_RTN;

    % Perturbed dynamics due to acceleration ad in ECI
    dXdt = f0 + B * ad;
    
end
