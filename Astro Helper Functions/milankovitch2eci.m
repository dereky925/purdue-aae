function [r_eci, v_eci] = milankovitch2eci(h_vec, e_vec, L)

    % Inputs:
    % h_vec  : 3x1 specific angular momentum vector [km^2/s]
    % e_vec  : 3x1 eccentricity vector (dimensionless)
    % L : true longitude, i.e. Omega + omega + f (radians)
    
    % Outputs:
    % r_eci  : 3x1 position in ECI frame [km]
    % v_eci  : 3x1 velocity in ECI frame [km/s]

    % Gravitational parameter for Earth [km^3/s^2]
    mu = 398600.4418;

    % Magnitudes
    h = norm(h_vec);
    e = norm(e_vec);
    
    % Semilatus rectum
    p = h^2 / mu;

    % Unit normal to orbital plane
    h_hat = h_vec / h;
    
    % Reference direction along z-axis
    k_hat = [0; 0; 1];

    n_vec = cross(k_hat, h_hat);
    n_norm = norm(n_vec);
    
    n_hat = n_vec / n_norm;
    Omega = atan2(n_hat(2), n_hat(1));

    e_plane = e_vec - dot(e_vec, h_hat)*h_hat;
    e_plane = e_plane / norm(e_plane);  % e_hat

    omega = atan2( dot(h_hat, cross(n_hat, e_plane)), dot(n_hat, e_plane) );
 
    % True anomaly
    f = L - Omega - omega;

    % semilatus rectum
    r_mag = p / (1 + e*cos(f));

    e_hat = e_plane;

    p_hat = cross(h_hat, e_hat);  % automatically in-plane, orthonormal if e_hat is unit

    r_eci = r_mag * ( cos(f)*e_hat + sin(f)*p_hat );
    v_eci = (mu/h) * ( -sin(f)*e_hat + (e + cos(f))*p_hat );


end