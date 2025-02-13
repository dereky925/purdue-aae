function dXdt = modified_equinoctial_elements_eom2(t, x, utcTime)
% MODIFIED_EQUINOCTIAL_ELEMENTS_EOM2
%
%  This function computes the time derivative of the modified equinoctial
%  elements including the J2 perturbation.
%
%  INPUTS:
%    t       - current time (s)
%    x       - state vector. On the very first call, it is expected to be
%              the classical orbital elements:
%                  x = [a; e; i; RAAN; omega; nu]
%              where:
%                  a     = semi-major axis (m)
%                  e     = eccentricity
%                  i     = inclination (deg)
%                  RAAN  = right ascension of ascending node (deg)
%                  omega = argument of perigee (deg)
%                  nu    = true anomaly (deg)
%
%              After the first call, x is assumed to be in the modified 
%              equinoctial element set:
%                  x = [p; f; g; h; k; L]
%              where L (the true longitude) is in radians.
%
%    utcTime - current UTC time (not used in this example, but provided 
%              for extensibility)
%
%  OUTPUT:
%    dXdt - time derivative of the state vector (modified equinoctial element
%           rates)
%
%  The function first converts classical orbital elements (in degrees) to
%  modified equinoctial elements using
%
%      p = a*(1-e^2)
%      f = e*cos(RAAN+omega)
%      g = e*sin(RAAN+omega)
%      h = tan(i/2)*cos(RAAN)
%      k = tan(i/2)*sin(RAAN)
%      L = RAAN + omega + nu
%
%  and then computes the equations of motion under the sole influence of 
%  the Earth's J2 perturbation.
%
%  If you wish to have the output in classical elements, you must convert 
%  the propagated modified equinoctial elements back to classical form.
%
%  NOTE: Because a persistent flag is used to perform the conversion only on
%  the first call, if you restart the integration within the same MATLAB
%  session, clear the function (e.g. "clear functions") to force a re–conversion.
%

    %% Constants (Earth parameters)
    mu = 3.986004418e14;    % Earth's gravitational parameter (m^3/s^2)
    Re = 6378137;           % Earth's equatorial radius (m)
    J2 = 1.08262668e-3;     % Earth's J2 coefficient

    %% Conversion from classical elements to modified equinoctial elements
    % The conversion is performed only on the first call.
    persistent conversionDone
    if isempty(conversionDone)
        % Assume input state x = [a; e; i; RAAN; omega; nu] with angles in degrees.
        deg2rad = pi/180;
        a    = x(1);
        e    = x(2);
        i    = x(3)   * deg2rad;
        RAAN = x(4)   * deg2rad;
        omega= x(5)   * deg2rad;
        nu   = x(6)   * deg2rad;
        
        % Convert classical to modified equinoctial elements:
        p = a*(1 - e^2);
        f = e * cos(RAAN + omega);
        g = e * sin(RAAN + omega);
        h = tan(i/2) * cos(RAAN);
        k = tan(i/2) * sin(RAAN);
        L = RAAN + omega + nu;  % true longitude (radians)
        
        % Replace the input state with the converted modified equinoctial state.
        x = [p; f; g; h; k; L];
        
        conversionDone = true;
    else
        % Assume state vector x is already in modified equinoctial elements:
        p = x(1);
        f = x(2);
        g = x(3);
        h = x(4);
        k = x(5);
        L = x(6);
    end

    %% Pre-compute trigonometric functions and auxiliary parameter w
    cosL = cos(L);
    sinL = sin(L);
    w = 1 + f*cosL + g*sinL;
    
    % Radial distance from the central body
    r_val = p / w;
    
    %% Transform from modified equinoctial elements to inertial coordinates
    % Define the transformation matrix Q.
    denom = 1 + h^2 + k^2;
    Q = 1/denom * [ 1 - h^2 + k^2,    2*h*k,       -2*h;
                    2*h*k,           1 + h^2 - k^2,  2*k;
                    2*h,             -2*k,          1 - h^2 - k^2 ];
                
    % Inertial position vector (m)
    r_ijk = Q * [r_val*cosL; r_val*sinL; 0];
    
    %% Compute the orbital (RTN) frame unit vectors
    % Radial unit vector
    r_hat = r_ijk / norm(r_ijk);
    % Transverse unit vector (direction of increasing L)
    t_hat = Q * [-sinL; cosL; 0];
    % Normal unit vector (completing the right–handed triad)
    n_hat = Q * [0; 0; 1];
    
    %% Compute the J2 perturbing acceleration in the inertial frame
    r_norm = norm(r_ijk);
    z = r_ijk(3);
    
    factor = (3/2) * J2 * mu * Re^2 / (r_norm^5);
    
    ax_J2 = factor * (r_ijk(1)/r_norm) * (5*(z/r_norm)^2 - 1);
    ay_J2 = factor * (r_ijk(2)/r_norm) * (5*(z/r_norm)^2 - 1);
    az_J2 = factor * (r_ijk(3)/r_norm) * (5*(z/r_norm)^2 - 3);
    
    a_J2 = [ax_J2; ay_J2; az_J2];
    
    % Project the inertial acceleration into the orbital (RTN) frame.
    A_r = dot(a_J2, r_hat);
    A_t = dot(a_J2, t_hat);
    A_n = dot(a_J2, n_hat);
    
    %% Equations of motion for modified equinoctial elements
    % The following equations (e.g., from Montenbruck & Gill or Vallado)
    % give the rates of change of the modified equinoctial elements.
    sqrt_p_over_mu = sqrt(p/mu);
    
    pdot = 2 * sqrt_p_over_mu * p * A_t / w;
    
    fdot = sqrt_p_over_mu * ( A_r*sinL + (((w+1)*cosL + f)/w)*A_t ...
                              - (g*(h*sinL - k*cosL)/w)*A_n );
    
    gdot = sqrt_p_over_mu * (-A_r*cosL + (((w+1)*sinL + g)/w)*A_t ...
                              + (f*(h*sinL - k*cosL)/w)*A_n );
    
    hdot = sqrt_p_over_mu * (w * A_n * cosL) / (2*p);
    
    kdot = sqrt_p_over_mu * (w * A_n * sinL) / (2*p);
    
    Ldot = sqrt(mu/p)*w^2 + sqrt_p_over_mu * ((h*sinL - k*cosL)/w)*A_n;
    
    %% Assemble the derivatives into the output vector
    dXdt = [pdot; fdot; gdot; hdot; kdot; Ldot];
    
end