function fullOrbitTransfer
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FULL ORBIT TRANSFER SIMULATION
    %
    % This script propagates the full six-dimensional orbital state,
    % including the fast true anomaly, using Gauss's variational equations.
    % A control input is computed to drive the slow orbital elements
    % [a, e, i, Omega, omega] to their target values. Two simulations
    % are performed: one without and one with a control-magnitude constraint.
    %
    % The control law is:
    %   u = - pinv(B_slow) * ( K*(x_slow - x_target) )
    % and the full dynamics are:
    %   dX/dt = f0(X) + B_full(X) * u,
    % where f0 = [0;0;0;0;0; h/r^2] with h = sqrt(mu*a(1-e^2)).
    % The 6x3 matrix B_full is computed from Gauss's equations.
    %
    % NOTE: Angles in the state X are stored in degrees; within the ODE,
    % they are converted to radians.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clc; clear; close all;

    %% 1) Constants and Time Span
    mu = 398600;             % Earth's gravitational parameter [km^3/s^2]
    DAY2SEC = 86400;         % seconds per day
    tspan = [0 10*DAY2SEC];   % simulate for 10 days

    %% 2) Initial and Target Orbital Elements
    % (State vector: X = [a; e; i; Omega; omega; nu])
    % Initial orbit:
    a0 = 10000;   % [km]
    e0 = 0.3;
    i0 = 10;      % [deg]
    Omega0 = 20;  % [deg]
    omega0 = 0;   % [deg]
    nu0 = 0;      % [deg]  (initial true anomaly)
    X0 = [a0; e0; i0; Omega0; omega0; nu0];

    % Target (desired) slow orbital elements:
    a_target = 12000;   % [km]
    e_target = 0.7;
    i_target = 130;     % [deg]
    Omega_target = 180; % [deg]
    omega_target = 270; % [deg]
    x_target = [a_target; e_target; i_target; Omega_target; omega_target];

    %% 3) Control Parameters
    K = eye(5);                % Gain matrix (5x5 identity)
    u_max = 0.1/1000;          % max acceleration: 0.1 m/s^2 -> km/s^2
    % Set up a parameter structure:
    param.mu = mu;
    param.x_target = x_target;
    param.K = K;
    param.u_max = u_max;  
    % flag_constrained = false for unconstrained; true to apply saturation.
    param.flag_constrained = false;

    %% 4) Integrate Full Dynamics (Unconstrained Control)
    disp('Running unconstrained simulation...');
    [t_uncon, X_uncon] = ode45(@(t,X) odeFullOrbitTransfer(t, X, param), tspan, X0);

    %% 5) Integrate Full Dynamics (Constrained Control)
    param.flag_constrained = true;
    disp('Running constrained simulation...');
    [t_con, X_con] = ode45(@(t,X) odeFullOrbitTransfer(t, X, param), tspan, X0);

    %% 6) Plot Evolution of the Slow State Errors (Unconstrained Example)
    fig1 = figure('Name','Slow State Errors (Unconstrained)','NumberTitle','off',...
                  'Color','w','Position',[100,100,1400,900]);
    fs = 20; lw = 2;
    % Compute error: (slow state minus target)
    err_uncon = X_uncon(:,1:5) - repmat(x_target.', size(X_uncon,1), 1);
    t_days_uncon = t_uncon / DAY2SEC;
    
    subplot(3,2,1); plot(t_days_uncon, err_uncon(:,1), 'LineWidth',lw);
    xlabel('Time (days)','FontSize',fs); ylabel('\delta a (km)','FontSize',fs); set(gca,'FontSize',fs);
    
    subplot(3,2,2); plot(t_days_uncon, err_uncon(:,2), 'LineWidth',lw);
    xlabel('Time (days)','FontSize',fs); ylabel('\delta e','FontSize',fs); set(gca,'FontSize',fs);
    
    subplot(3,2,3); plot(t_days_uncon, err_uncon(:,3), 'LineWidth',lw);
    xlabel('Time (days)','FontSize',fs); ylabel('\delta i (deg)','FontSize',fs); set(gca,'FontSize',fs);
    
    subplot(3,2,4); plot(t_days_uncon, err_uncon(:,4), 'LineWidth',lw);
    xlabel('Time (days)','FontSize',fs); ylabel('\delta \Omega (deg)','FontSize',fs); set(gca,'FontSize',fs);
    
    subplot(3,2,5); plot(t_days_uncon, err_uncon(:,5), 'LineWidth',lw);
    xlabel('Time (days)','FontSize',fs); ylabel('\delta \omega (deg)','FontSize',fs); set(gca,'FontSize',fs);
    
    sgtitle('Slow State Errors Over Time (Unconstrained)','FontSize',fs);

    %% 7) 3D Orbit Plots
    % For plotting the orbit shapes we generate the full ellipse from the
    % orbital elements (ignoring the fast evolution of nu).
    fig2 = figure('Name','3D Orbits','NumberTitle','off',...
                  'Color','w','Position',[200,200,1200,800]);
    nu_plot = linspace(0, 2*pi, 300);   % [rad] sweep for plotting

    % Target orbit: convert target COEs to ECI positions.
    coe_target = [x_target(1), x_target(2), x_target(3), x_target(4), x_target(5)];
    r_target = coe2eci(coe_target, nu_plot);
    
    % For the spacecraft, take the final slow states from the unconstrained run.
    X_final = X_uncon(end, :).';  % full state at final time
    coe_sc = [X_final(1), X_final(2), X_final(3), X_final(4), X_final(5)];
    r_sc = coe2eci(coe_sc, nu_plot);
    
    plot3(r_target(1,:), r_target(2,:), r_target(3,:), 'b','LineWidth',lw); hold on;
    plot3(r_sc(1,:), r_sc(2,:), r_sc(3,:), 'r--','LineWidth',lw);
    xlabel('X (km)','FontSize',fs); ylabel('Y (km)','FontSize',fs); zlabel('Z (km)','FontSize',fs);
    set(gca,'FontSize',fs); grid on; axis equal;
    legend('Target Orbit','Spacecraft Orbit (Unconstrained)', 'FontSize',fs, 'Location','best');
    title('3D Orbit Comparison (Unconstrained)','FontSize',fs);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ODE Function: Full Orbit Dynamics with Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dXdt = odeFullOrbitTransfer(~, X, param)
    % The state vector is:
    %   X = [a; e; i; Omega; omega; nu]
    % with a in [km], e unitless, and i, Omega, omega, nu in [deg].
    %
    % Within this function we convert angles to radians for the dynamics.
    
    % Extract states:
    a       = X(1);
    e       = X(2);
    i_deg   = X(3);
    Omega_deg = X(4);
    omega_deg = X(5);
    nu_deg  = X(6);
    
    % Convert angles from degrees to radians for computation:
    i     = deg2rad(i_deg);
    Omega = deg2rad(Omega_deg);
    omega = deg2rad(omega_deg);
    nu    = deg2rad(nu_deg);
    
    % Compute derived quantities:
    p = a * (1 - e^2);
    r = p / (1 + e*cos(nu));
    b = a * sqrt(1 - e^2);
    
    % Compute the natural (Keplerian) drift.
    % In a Kepler orbit, the slow elements are constant so:
    f0_slow = zeros(5,1);
    % The true anomaly evolves at the Keplerian rate:
    h = sqrt(param.mu * p);  % specific angular momentum [km^2/s]
    f0_nu = h / (r^2);        % [rad/s]
    % Assemble f0 (note: f0_nu is in rad/s)
    f0 = [f0_slow; f0_nu];
    
    % Compute the full 6x3 B matrix (derived from Gauss's equations)
    % The provided formulation:
    B_full = (1/sqrt(param.mu*p)) * [ ...
        2*a^2*e*sin(nu),           2*a^2*p/r,                  0;
        p*sin(nu),                 (p+r)*cos(nu)+r*e,            0;
        0,                         0,                          r*cos(nu+omega);
        0,                         0,                          r*sin(nu+omega)/sin(i);
       -p*cos(nu)/e,               (p+r)*sin(nu)/e,           -r*sin(nu+omega)/tan(i);
        b*p*cos(nu)/(a*e)-2*b*r/a,  -b*(p+r)*sin(nu)/(a*e),       0 ];
    
    % Separate the slow part from the full B matrix:
    B_slow = B_full(1:5, :);  % rows corresponding to [a, e, i, Omega, omega]
    
    % Compute the control input.
    % Use the error in the slow elements:
    x_slow = [a; e; i_deg; Omega_deg; omega_deg];
    error = x_slow - param.x_target;
    
    % Control law (we omit any natural drift in slow elements because f0_slow = 0)
    u_raw = - pinv(B_slow) * (param.K * error);
    
    % If control magnitude constraint is enabled, saturate the control:
    if param.flag_constrained && (norm(u_raw) > param.u_max)
        u = u_raw * (param.u_max / norm(u_raw));
    else
        u = u_raw;
    end
    
    % Compute the full state derivative:
    dXdt = f0 + B_full * u;
    
    % dX/dt for the slow angular states (i, Omega, omega) and nu have been
    % computed in rad/s; convert these rates to deg/s (since the state is in deg)
    dXdt(3:6) = rad2deg(dXdt(3:6));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function: Convert COEs to ECI Positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r_eci = coe2eci(coe, nu_vec)
    % Inputs:
    %   coe = [a, e, i, Omega, omega]
    %         where a is [km], e is unitless, and i, Omega, omega are in [deg]
    %   nu_vec = vector of true anomaly values [rad] over which to generate the orbit
    % Output:
    %   r_eci: 3 x N matrix of ECI position vectors [km]
    
    a = coe(1);
    e = coe(2);
    i = deg2rad(coe(3));
    Omega = deg2rad(coe(4));
    omega = deg2rad(coe(5));
    
    nPts = length(nu_vec);
    r_eci = zeros(3, nPts);
    for k = 1:nPts
        nu = nu_vec(k);
        r = a*(1-e^2) / (1+e*cos(nu));
        % Position in perifocal coordinates:
        r_pf = [ r*cos(nu); r*sin(nu); 0 ];
        % Rotation matrix from perifocal to ECI coordinates:
        R = [ cos(Omega)*cos(omega)-sin(Omega)*sin(omega)*cos(i), -cos(Omega)*sin(omega)-sin(Omega)*cos(omega)*cos(i),  sin(Omega)*sin(i);
              sin(Omega)*cos(omega)+cos(Omega)*sin(omega)*cos(i), -sin(Omega)*sin(omega)+cos(Omega)*cos(omega)*cos(i), -cos(Omega)*sin(i);
              sin(omega)*sin(i),                                  cos(omega)*sin(i),                                  cos(i) ];
        r_eci(:,k) = R * r_pf;
    end
end

% Call the main function if the script is run
fullOrbitTransfer