%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% HW1 Q1
% Author: Travis Hastreiter 
% Created On: 25 January, 2024
% Description: Simulate orbits with J2 and SRP perturbations
% Most Recent Change: 25 January, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assumed dynamical parameter values
R_E = 6378.1; % [km] Earth radius
mu_E = 398600; % [km3 / s2] Earth gravitational parameter

% Initial conditions for s/c in Earth orbit (in ECI frame)
a_0 = 8000; % [km] semi-major axis
e_0 = 0.15; % [] eccentricity
i_0 = deg2rad(80); % [rad] inclination
Omega_0 = deg2rad(30); % [rad] right ascension of ascending node
omega_0 = deg2rad(20); % [rad] argument of periapsis
nu_0 = deg2rad(0); % [rad] true anomaly at epoch

M_0 = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu_0, e_0), e_0);
x0_keplerian = [a_0; e_0; i_0; Omega_0; omega_0; M_0];

tspan = linspace(0, 15 * period(a_0, mu_E), 1e5);
t_orbits = linspace(0, 15, numel(tspan));

default_tolerance = 1e-12;

%% a) Simulate orbit under Earth J2 perturbation using different representations 
%a_d_RTN = @(x) [0; 0; 0.001]; 
J_2_val = 1.0826e-3; % [] Earth J2

a_d_J2 = @(x) J2_perturbation(x(1:3), mu_E, R_E, J_2_val);

% Simulate 
[x_cartesians, x_elements, ~] = propagate_orbit_cart_kepl_modeq_mila(x0_keplerian, nu_0, mu_E, a_d_J2, tspan, default_tolerance);

%% a) Plot 

earthy(R_E, "Earth", 0.25); hold on;
axis equal

plot_cartesian_orbit(x_cartesians.cartesian, 'r', 0.3, 1); hold on;
plot_cartesian_orbit(x_cartesians.keplerian, 'g', 0.3, 1); hold on;
plot_cartesian_orbit(x_cartesians.modified_equinoctial, 'b', 0.3, 1); hold on;
plot_cartesian_orbit(x_cartesians.milankovitch, 'c', 0.3, 1); hold off;

title("Orbit Propagated with J2 for 15 Orbits Using Different Orbital Elements")
subtitle("Travis Hastreiter")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("", "", "Cartesian", "", "Keplerian", "", "Modified Equinoctial", "", "Milankovitch")

%% b) Evaluate Milankovitch Constraint
% Constraint is hvec is orthogonal to evec
hvec = x_elements.milankovitch(:, 1:3);
evec = x_elements.milankovitch(:, 4:6);
milankovitch_constraint_violation = dot(hvec, evec, 2);
milankovitch_constraint_violation_percent = milankovitch_constraint_violation ./ (vecnorm(hvec, 2, 2) .* vecnorm(evec, 2, 2)) * 100;

plot(t_orbits, milankovitch_constraint_violation)
title("Milankovitch Constraint Violation")
subtitle("Travis Hastreiter")
xlabel("Time [orbits]")
ylabel("Constraint Violation $\vec{h}\cdot\vec{e}$", Interpreter="latex")

plot(t_orbits, milankovitch_constraint_violation_percent)
title("Milankovitch Constraint Violation Percent")
subtitle("Travis Hastreiter")
xlabel("Time [orbits]")
ylabel("Constraint Violation $\frac{\vec{h}\cdot\vec{e}}{he}$", Interpreter="latex")

%% c) Compare with 1e-9, 1e-6, 1e-3 Error Tolerances
[x_cartesians_eneg9, ~, ~] = propagate_orbit_cart_kepl_modeq_mila(x0_keplerian, nu_0, mu_E, a_d_J2, tspan, 1e-9);

[x_cartesians_eneg6, ~, ~] = propagate_orbit_cart_kepl_modeq_mila(x0_keplerian, nu_0, mu_E, a_d_J2, tspan, 1e-6);

[x_cartesians_eneg3, ~, ~] = propagate_orbit_cart_kepl_modeq_mila(x0_keplerian, nu_0, mu_E, a_d_J2, tspan, 1e-3);

% Do some error comparison and graphs - compare to 1e-12

earthy(R_E, "Earth", 0.25); hold on;
axis equal

plot_cartesian_orbit(x_cartesians_eneg6.cartesian, 'r', 0.3, 1); hold on;
plot_cartesian_orbit(x_cartesians_eneg6.keplerian, 'g', 0.3, 1); hold on;
plot_cartesian_orbit(x_cartesians_eneg6.modified_equinoctial, 'b', 0.3, 1); hold on;
plot_cartesian_orbit(x_cartesians_eneg6.milankovitch, 'c', 0.3, 1); hold off;

title("Orbit Propagated with J2 for 15 Orbits Using Different Orbital Elements")
subtitle("Travis Hastreiter")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("", "", "Cartesian", "", "Keplerian", "", "Modified Equinoctial", "", "Milankovitch")


%% d) Simulate orbit under Cannonball SRP and Earth J2 perturbations using different representations 
A_over_m = 5.4e-6; % [km2 / kg] spacecraft area to mass ratio
G_0 = 1.02e14; % [kg km / s2] solar flux constant
AU = 149597898; % [km] astronautical unit, Earth-Sun difference

a_d_J2_SRP = @(x) J2_perturbation(rvec, mu_E, R_E, J_2_val) ... 
                + SRP_perturbation(A_over_m, G_0, AU, svec); % CALCULATE Svec

% Simulate 
[x_cartesians_SRP, ~, ~] = propagate_orbit_cart_kepl_modeq_mila(x0_keplerian, nu_0, mu_E, a_d_J2_SRP, tspan, default_tolerance);

%% d) Plot 

earthy(R_E, "Earth", 0.25); hold on;
axis equal

plot_cartesian_orbit(x_cartesians_SRP.cartesian, 'r', 0.3, 1); hold on;
plot_cartesian_orbit(x_cartesians_SRP.keplerian, 'g', 0.3, 1); hold on;
plot_cartesian_orbit(x_cartesians_SRP.modified_equinoctial, 'b', 0.3, 1); hold on;
plot_cartesian_orbit(x_cartesians_SRP.milankovitch, 'c', 0.3, 1); hold off;

title("Orbit Propagated with J2 and SRP for 15 Orbits Using Different Orbital Elements")
subtitle("Travis Hastreiter")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("", "", "Cartesian", "", "Keplerian", "", "Modified Equinoctial", "", "Milankovitch")

%% e) Simulate orbit under Cannonball SRP, Earth J2, and drag perturbations using different representations 
C_D = 2.1; % [] drag coefficient

a_d_J2_SRP_drag = @(x) J2_perturbation(rvec, mu_E, R_E, J_2_val) ... 
                     + SRP_perturbation(A_over_m, G_0, AU, svec) ... % CALCULATE Svec
                     + drag_perturbation(rvec, vvec, R_E, C_D, A_over_m);

% Simulate 
[x_cartesians_SRP_drag, ~, ~] = propagate_orbit_cart_kepl_modeq_mila(x0_keplerian, nu_0, mu_E, a_d_J2_SRP_drag, tspan, default_tolerance);

%% e) Plot 

earthy(R_E, "Earth", 0.25); hold on;
axis equal

plot_cartesian_orbit(x_cartesians_SRP_drag.cartesian, 'r', 0.3, 1); hold on;
plot_cartesian_orbit(x_cartesians_SRP_drag.keplerian, 'g', 0.3, 1); hold on;
plot_cartesian_orbit(x_cartesians_SRP_drag.modified_equinoctial, 'b', 0.3, 1); hold on;
plot_cartesian_orbit(x_cartesians_SRP_drag.milankovitch, 'c', 0.3, 1); hold off;

title("Orbit Propagated with J2 and SRP for 15 Orbits Using Different Orbital Elements")
subtitle("Travis Hastreiter")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("", "", "Cartesian", "", "Keplerian", "", "Modified Equinoctial", "", "Milankovitch")


%% Helper Functions

function [f_0] = f0_cartesian(x, mu)
    rvec = x(1:3);
    vvec = x(4:6);
    r = norm(rvec);

    a_cartesian = -mu ./ r .^ 3 * rvec;

    f_0 = [vvec; a_cartesian];
end

function [f_0] = f0_keplerian(x, mu)
    a = x(1);

    n = sqrt(mu / a ^ 3);

    f_0 = [zeros([5, 1]); ...
           n];
end

function [f_0] = f0_modified_equinoctial(x, mu)
    p = x(1, :);
    f = x(2, :);
    g = x(2, :);
    L = x(6, :);

    q = 1 + f .* cos(L) + g .* sin(L);

    f_0 = [zeros(5, 1); ...
           sqrt(mu * p) .* (q ./ p) .^ 2];
end

function [f_0] = f0_milankovitch(x, mu)
    hvec = x(1:3);
    Khat = [0; 0; 1];
    Ihat = [1; 0; 0];

    [x_keplerian, nu] = milankovitch_to_keplerian(x, Khat, Ihat, mu);
    [rvec, r] = rvec_from_keplerian(x_keplerian, nu);

    h = norm(hvec);

    f_0 = [zeros(6, 1); ...
           h ./ r^2];
end

function [B] = B_cartesian(x, mu)
    B = [zeros(3); eye(3)];
end

function [B] = B_keplerian(x, mu)
    a = x(1);
    e = x(2);
    i = x(3);
    Omega = x(4);
    omega = x(5);
    M = x(6);

    p = a * (1 - e ^ 2);
    b = a * sqrt(1 - e ^ 2);
    E = mean_to_eccentric_anomaly(M, e);
    nu = eccentric_to_true_anomaly(E, e);
    r = a * (1 - e * cos(E));
    h = sqrt(mu * p);

    B = 1 / h * [2 * a ^ 2 * e * sin(nu), 2 * a ^ 2 * p / r, 0;
         p * sin(nu), (p + r) * cos(nu) + r * e, 0;
         0, 0, r * cos(nu + omega);
         0, 0, r * sin(nu + omega) / sin(i);
         -p * cos(nu) / e, (p + r) * sin(nu) / e, -r * sin(nu + omega) / tan(i);
         b * p * cos(nu) / (a * e) - 2 * b * r / a, -b * (p + r) * sin(nu) / (a * e), 0];
    
    % Make B convert disturbance into RTN frame from cartesian frame
    B = B * cartesian_to_RTN_DCM(i, Omega, omega, nu)';
end

function [B] = B_modified_equinoctial(x, mu)
    p = x(1);
    f = x(2);
    g = x(3);
    h = x(4);
    k = x(5);
    L = x(6);

    q = 1 + f * cos(L) + g * sin(L);
    s_sqr = 1 + h ^ 2 + k ^ 2;

    cons_1 = (h * sin(L) - k * cos(L)) / q;
    B = sqrt(p / mu) * [0, 2 * p / q, 0;
        sin(L), ((q + 1) * cos(L) + f) / q, -g * cons_1;
        -cos(L), ((q + 1) * sin(L) + g) / q, f * cons_1;
        0, 0, s_sqr / (2 * q) * cos(L);
        0, 0, s_sqr / (2 * q) * sin(L);
        0, 0, cons_1];

    % Make B convert disturbance into RTN frame from cartesian frame
    i = atan2(2 * sqrt(h ^ 2 + k ^ 2), 1 - h ^ 2 - k ^ 2);
    Omega = atan2(k, h);
    omega = atan2(g * h - f * k, f * h + g * k);
    nu = L - (Omega + omega);
    B = B * cartesian_to_RTN_DCM(i, Omega, omega, nu)';
end

function [B] = B_milankovitch(x, mu)
    hvec = x(1:3);
    evec = x(4:6);
    L = x(7);

    Khat = [0; 0; 1];
    Ihat = [1; 0; 0];

    [x_keplerian, nu] = milankovitch_to_keplerian(x, Khat, Ihat, mu);
    x_cartesian = keplerian_to_cartesian(x_keplerian, nu, mu);
    rvec = x_cartesian(1:3);
    vvec = x_cartesian(4:6);
    
    h = norm(hvec);

    rcross = cross_operator(rvec);
    vcross = cross_operator(vvec);
    hcross = cross_operator(hvec);

    B = [rcross; 
        1 / mu * (vcross * rcross - hcross); 
        rvec(3) / (h * (h + hvec(3))) * hvec'];
end

function [x_dot] = gauss_planetary_eqn(f_0, B, a_d)
    x_dot = f_0 + B * a_d;
end

function [a_J2] = J2_perturbation(rvec, mu, r_o, J_2)
    r = norm(rvec);
    x = rvec(1);
    y = rvec(2);
    z = rvec(3);

    cons = 1 - 5 * z .^ 2 ./ r .^ 2;

    a_J2 = -3 * mu * J_2 * r_o ^ 2 ./ (2 * r .^ 5) ...
        .* ([cons .* x; cons .* y; (2 + cons) .* z]);
end

function [a_SRP] = SRP_perturbation(A_over_m, G_0, d, svec)
    a_SRP = A_over_m * G_0 / d ^ 2 * svec;
end

function [a_drag] = drag_perturbation(rvec, vvec, r_o, C_D, A_over_m)
    v = norm(vvec);
    r = norm(rvec);
    alt = r - r_o;

    rho = 1.02e7 * alt ^ -7.172 * 1e9; % Formula (up to 1000 km ish) from 

    a_drag = -1 / 2 * rho * v * vvec * A_over_m * C_D;
end

function [x_cartesians, x_elements, t] = propagate_orbit_cart_kepl_modeq_mila(x0_keplerian, nu_0, mu, a_d, tspan, error_tolerance)
    tolerances = odeset(RelTol=error_tolerance, AbsTol=error_tolerance);

    %% i) Cartesian
    x0_cartesian = keplerian_to_cartesian(x0_keplerian, nu_0, mu);

    [t_cartesian,x_cartesian] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, mu), B_cartesian(x, mu), a_d(x)), tspan, x0_cartesian, tolerances);
    
    %% ii) Keplerian Orbital Elements
    
    [t_keplerian,x_keplerian] = ode45(@(t,x) gauss_planetary_eqn(f0_keplerian(x, mu), B_keplerian(x, mu), a_d(keplerian_to_cartesian(x, [], mu))), tspan, x0_keplerian, tolerances);
    x_keplerian_cartesian = keplerian_to_cartesian_array(x_keplerian, [], mu);
    
    %plot(linspace(0, 15, numel(t_keplerian)), x_keplerian(:, 1))
    
    %% iii) Modified Equinoctial Orbital Elements
    x0_modified_equinoctial = keplerian_to_modified_equinoctial(x0_keplerian, nu_0);
    
    [t_modified_equinoctial,x_modified_equinoctial] = ode45(@(t,x) gauss_planetary_eqn(f0_modified_equinoctial(x, mu), B_modified_equinoctial(x, mu), a_d(modified_equinoctial_to_cartesian(x, mu))), tspan, x0_modified_equinoctial, tolerances);
    x_modified_equinoctial_cartesian = modified_equinoctial_to_cartesian_array(x_modified_equinoctial, mu);
    
    %plot(linspace(0, 15, numel(t_modified_equinoctial)), x_modified_equinoctial(:, 1))
    
    %% iv) Milankovitch Orbital Elements
    Khat = [0; 0; 1];
    Ihat = [1; 0; 0];
    
    x0_milankovitch = cartesian_to_milankovitch(x0_cartesian,Khat,Ihat,mu)
    
    [t_milankovitch,x_milankovitch] = ode45(@(t,x) gauss_planetary_eqn(f0_milankovitch(x, mu), B_milankovitch(x, mu), a_d(milankovitch_to_cartesian(x, Khat, Ihat, mu))), tspan, x0_milankovitch, tolerances);
    x_milankovitch_cartesian = milankovitch_to_cartesian_array(x_milankovitch, Khat, Ihat, mu);

    %% Package Output
    x_cartesians.cartesian = x_cartesian;
    x_cartesians.keplerian = x_keplerian_cartesian;
    x_cartesians.modified_equinoctial = x_modified_equinoctial_cartesian;
    x_cartesians.milankovitch = x_milankovitch_cartesian;

    x_elements.cartesian = x_cartesian;
    x_elements.keplerian = x_keplerian;
    x_elements.modified_equinoctial = x_modified_equinoctial;
    x_elements.milankovitch = x_milankovitch;

    t = t_cartesian;
end