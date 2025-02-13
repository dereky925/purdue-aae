clc; clear; close all;

%% MAIN DRIVER SCRIPT

%% Define Common (Earth) Parameters
common.mu   = 3.986004418e14;   % [m^3/s^2] Earth's gravitational parameter
common.Re   = 6371e3;           % [m]     Earth's radius
common.rho0 = 1.225;            % [kg/m^3] Sea-level air density
common.H    = 8500;             % [m]     Atmospheric scale height

%% Define Stage‑Specific Parameters
% Stage 1 (Booster) Parameters:
stage1.m0     = 2e6;           % [kg] Initial mass at liftoff
stage1.T_mag  = 20e6;          % [N]  Engine thrust
stage1.Isp    = 350;           % [s]  Specific impulse
stage1.t_burn = 300;           % [s]  Burn time
stage1.g0     = 9.81;          % [m/s^2] Standard gravity
stage1.A_ref  = pi*(9/2)^2;      % [m^2] Reference area (from a 9-m diameter)
stage1.Cd     = 0.3;           % Drag coefficient

% Stage 2 (Upper Stage) Parameters:
stage2.m0     = 2.5e5;         % [kg] Mass at booster separation (upper stage mass)
stage2.T_mag  = 5.5e6;           % [N]  Engine thrust (for orbit insertion)
stage2.Isp    = 600;           % [s]  Specific impulse
stage2.t_burn = 200;           % [s]  Burn time
stage2.g0     = 9.81;          % [m/s^2] Standard gravity
stage2.A_ref  = stage1.A_ref;  % Use the same reference area for simplicity
stage2.Cd     = stage1.Cd;     % Same drag coefficient

% Coast Phase: No engine thrust.
coast.T_mag  = 0;            % [N] Zero thrust
coast.g0     = 9.81;
coast.A_ref  = stage1.A_ref;
coast.Cd     = stage1.Cd;

%% Launch-Site and Initial Conditions
% Starbase, Texas (approximate LLA)
lat = deg2rad(25.997);       % [rad] Latitude
lon = deg2rad(-97.156);      % [rad] Longitude (negative for West)
alt_launch = 0;              % [m] Sea level
r0 = (common.Re + alt_launch) * [cos(lat)*cos(lon); cos(lat)*sin(lon); sin(lat)];

v0 = [0; 0; 0];              % [m/s] Initial velocity (assumed at rest)
% Set initial attitude so that the vehicle’s body‐z (thrust axis) aligns with the local vertical
r0_unit = r0 / norm(r0);
phi0   = 0;
theta0 = acos(r0_unit(3));
psi0   = atan2(r0_unit(2), r0_unit(1));

% Assume zero angular rates.
p0 = 0;  q0 = 0;  r0_ang = 0;
% Assemble the state vector:
% [ position(3); velocity(3); Euler angles (phi, theta, psi); angular rates (p, q, r); mass ]
state0 = [r0; v0; phi0; theta0; psi0; p0; q0; r0_ang; stage1.m0];

%% PHASE 1: Stage 1 (Booster) Burn
tspan1 = [0, stage1.t_burn];
options = odeset('RelTol',1e-8, 'AbsTol',1e-9);
[T1, X1] = ode45(@(t,state) dynamics_stage1(t, state, stage1, common), tspan1, state0, options);
state_stage1_end = X1(end,:)';
t_stage1_end = T1(end);

%% STAGE SEPARATION and Initialization of Stage 2
% Jettison the booster by resetting the mass to stage 2's value.
state_stage2_init = state_stage1_end;
state_stage2_init(13) = stage2.m0;

%% PHASE 2: Stage 2 (Upper Stage) Burn
% In Stage 2 we force the thrust to be directed along the horizontal component of the velocity.
tspan2_rel = [0, stage2.t_burn];
[T2_rel, X2] = ode45(@(t,state) dynamics_stage2(t, state, stage2, common), tspan2_rel, state_stage2_init, options);
T2 = T2_rel + t_stage1_end;
t_stage2_end = T2(end);
state_stage2_end = X2(end,:)';

%% PHASE 3: Coast (No Thrust)
t_final = 10800;  % [s] Total simulation time (~3 hrs to show ~2 orbital revolutions)
t_coast = t_final - t_stage2_end;
tspan3 = [0, t_coast];
[T3_rel, X3] = ode45(@(t,state) dynamics_coast(t, state, coast, common), tspan3, X2(end,:)', options);
T3 = T3_rel + t_stage2_end;

%% Combine All Phases
T_all = [T1; T2; T3];
X_all = [X1; X2; X3];
x  = X_all(:,1);  y  = X_all(:,2);  z  = X_all(:,3);
vx = X_all(:,4);  vy = X_all(:,5);  vz = X_all(:,6);
r_norm   = sqrt(x.^2 + y.^2 + z.^2);
altitude = r_norm - common.Re;  % [m]
speed    = sqrt(vx.^2 + vy.^2 + vz.^2);  % [m/s]
rho = common.rho0 * exp(-altitude/common.H);
rho(altitude > 100e3) = 0;
q_dyn = 0.5 * rho .* speed.^2;

%% PLOTTING

% Figure 1: Altitude, Speed, and Dynamic Pressure vs. Time (placed in top left)
figure('Position', [50, 600, 1200, 800]);
subplot(3,1,1);
plot(T_all, altitude/1000, 'LineWidth', 1.5);
hold on;
xline(t_stage1_end, 'k--', 'Booster Separation', 'LineWidth', 2, 'FontSize', 14);
xline(t_stage2_end, 'r--', 'Upper Stage Cutoff', 'LineWidth', 2, 'FontSize', 14);
xlabel('Time (s)', 'FontSize', 14); ylabel('Altitude (km)', 'FontSize', 14);
title('Altitude vs. Time', 'FontSize', 16); grid on;
subplot(3,1,2);
plot(T_all, speed, 'LineWidth', 1.5);
hold on;
xline(t_stage1_end, 'k--', 'Booster Separation', 'LineWidth', 2, 'FontSize', 14);
xline(t_stage2_end, 'r--', 'Upper Stage Cutoff', 'LineWidth', 2, 'FontSize', 14);
xlabel('Time (s)', 'FontSize', 14); ylabel('Speed (m/s)', 'FontSize', 14);
title('Speed vs. Time', 'FontSize', 16); grid on;
subplot(3,1,3);
plot(T_all, q_dyn, 'LineWidth', 1.5);
hold on;
xline(t_stage1_end, 'k--', 'Booster Separation', 'LineWidth', 2, 'FontSize', 14);
xline(t_stage2_end, 'r--', 'Upper Stage Cutoff', 'LineWidth', 2, 'FontSize', 14);
xlabel('Time (s)', 'FontSize', 14); ylabel('Dynamic Pressure (Pa)', 'FontSize', 14);
title('Dynamic Pressure vs. Time', 'FontSize', 16); grid on;

% Figure 2: 3D Trajectory (placed in top right)
figure('Position', [1300, 600, 1200, 800]);
[Xe,Ye,Ze] = sphere(50);
surf(common.Re*Xe, common.Re*Ye, common.Re*Ze, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
colormap([0 0 1]);  % Blue Earth
hold on;
plot3(x, y, z, 'r', 'LineWidth', 2);
scatter3(state0(1), state0(2), state0(3), 120, 'g', 'filled');
text(state0(1), state0(2), state0(3), '  Initial Position', 'FontSize', 14, 'Color', 'k');
scatter3(state_stage2_end(1), state_stage2_end(2), state_stage2_end(3), 120, 'r', 'filled');
text(state_stage2_end(1), state_stage2_end(2), state_stage2_end(3), '  Stage 2 Cutoff', 'FontSize', 14, 'Color', 'k');
axis equal;
xlabel('X (m)', 'FontSize', 14); ylabel('Y (m)', 'FontSize', 14); zlabel('Z (m)', 'FontSize', 14);
title('3D Trajectory', 'FontSize', 16); grid on;

% Figure 3: Flight–Path (Elevation) Angles for Stage 1 and Stage 2 (placed in bottom left)
% Elevation angle: gamma = asin( (v dot (r/|r|)) / |v| )
gamma1 = zeros(length(T1),1);
for i = 1:length(T1)
    r_vec = X1(i,1:3)'; v_vec = X1(i,4:6)';
    if norm(v_vec) > 1e-3
        gamma1(i) = asin(dot(v_vec, r_vec/norm(r_vec))/norm(v_vec));
    else
        gamma1(i) = pi/2;
    end
end
gamma1_deg = gamma1 * (180/pi);

gamma2 = zeros(length(T2),1);
for i = 1:length(T2)
    r_vec = X2(i,1:3)'; v_vec = X2(i,4:6)';
    if norm(v_vec) > 1e-3
        gamma2(i) = asin(dot(v_vec, r_vec/norm(r_vec))/norm(v_vec));
    else
        gamma2(i) = pi/2;
    end
end
gamma2_deg = gamma2 * (180/pi);

figure('Position', [50, 50, 1200, 600]);
subplot(2,1,1);
plot(T1, gamma1_deg, 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 14); ylabel('Elevation Angle (deg)', 'FontSize', 14);
title('Stage 1 (Booster) Elevation Angle', 'FontSize', 16); grid on;
subplot(2,1,2);
plot(T2, gamma2_deg, 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 14); ylabel('Elevation Angle (deg)', 'FontSize', 14);
title('Stage 2 Elevation Angle', 'FontSize', 16); grid on;

%% LOCAL FUNCTIONS

function dstate = dynamics_stage1(t, state, stage, common)
    % STAGE 1 DYNAMICS: Booster Burn with Aggressive Pitch-Over
    % This version commands a very rapid pitch-over so that the desired
    % elevation angle drops to near 0 (horizontal) by about 20% of the burn time.
    
    r = state(1:3);
    v = state(4:6);
    phi = state(7);
    theta = state(8);
    psi = state(9);
    m = state(13);
    
    r_norm = norm(r);
    g_acc = -common.mu * r / r_norm^3;
    
    % Atmospheric density (exponential model)
    altitude_local = r_norm - common.Re;
    if altitude_local < 100e3
        rho_local = common.rho0 * exp(-altitude_local/common.H);
    else
        rho_local = 0;
    end
    
    v_norm = norm(v);
    if v_norm > 1e-3
        F_drag_mag = 0.5 * stage.Cd * stage.A_ref * rho_local * v_norm^2;
        F_drag = -F_drag_mag * (v/v_norm);
        % Compute the current flight–path (elevation) angle relative to local vertical:
        alpha = acos(dot(v, r)/(v_norm*r_norm));
        current_gamma = pi/2 - alpha;
    else
        F_drag = [0;0;0];
        current_gamma = pi/2;
    end
    
    % Aggressive booster pitch program:
    % Command the desired elevation angle to drop rapidly toward 0.
    if t < stage.t_burn * 0.2
        desired_gamma = (pi/2) * (1 - 20*t/stage.t_burn);
        if desired_gamma < 0
            desired_gamma = 0;
        end
    else
        desired_gamma = 0;
    end
    gamma_error = desired_gamma - current_gamma;
    % Increase the gain to force a rapid change (note the gain of 2 here)
    gimbal_pitch = max(min(2 * gamma_error, deg2rad(5)), -deg2rad(5));
    
    T = stage.T_mag;
    thrust_body = [T * sin(gimbal_pitch); 0; T * cos(gimbal_pitch)];
    R_mat = euler_to_dcm(phi, theta, psi);
    T_inertial = R_mat * thrust_body;
    
    F_total = T_inertial + F_drag;
    a = F_total/m + g_acc;
    
    dphi = 0; dtheta = 0; dpsi = 0;
    dp = 0; dq = 0; dr_ang = 0;
    m_dot = -T/(stage.Isp * stage.g0);
    
    dstate = zeros(13,1);
    dstate(1:3) = v;
    dstate(4:6) = a;
    dstate(7:9) = [dphi; dtheta; dpsi];
    dstate(10:12) = [dp; dq; dr_ang];
    dstate(13) = m_dot;
end

function dstate = dynamics_stage2(t, state, stage, common)
    % STAGE 2 DYNAMICS: Upper Stage Burn with Forced Reorientation
    % In this version, we force the thrust direction to be the projection of
    % the current velocity onto the local horizontal plane so that the vehicle
    % accelerates tangentially (i.e. with near-zero elevation).
    r = state(1:3);
    v = state(4:6);
    m = state(13);
    
    r_norm = norm(r);
    unit_r = r / r_norm;
    g_acc = -common.mu * r / r_norm^3;
    
    altitude_local = r_norm - common.Re;
    if altitude_local < 100e3
        rho_local = common.rho0 * exp(-altitude_local/common.H);
    else
        rho_local = 0;
    end
    
    v_norm = norm(v);
    if v_norm > 1e-3
        F_drag_mag = 0.5 * stage.Cd * stage.A_ref * rho_local * v_norm^2;
        F_drag = -F_drag_mag * (v/v_norm);
    else
        F_drag = [0;0;0];
    end
    
    % Compute the horizontal (tangential) component of the velocity:
    v_horiz = v - (dot(v, unit_r))*unit_r;
    if norm(v_horiz) > 1e-3
        thrust_dir = v_horiz / norm(v_horiz);
    else
        % If the velocity is nearly vertical, default to an arbitrary horizontal direction (e.g., east)
        east = cross([0;0;1], unit_r);
        if norm(east) < 1e-3
            east = [1;0;0];
        else
            east = east / norm(east);
        end
        thrust_dir = east;
    end
    
    % Force the thrust to be exactly horizontal (tangential)
    T_inertial = stage.T_mag * thrust_dir;
    
    F_total = T_inertial + F_drag;
    a = F_total/m + g_acc;
    
    dphi = 0; dtheta = 0; dpsi = 0;
    dp = 0; dq = 0; dr_ang = 0;
    m_dot = -stage.T_mag/(stage.Isp * stage.g0);
    
    dstate = zeros(13,1);
    dstate(1:3) = v;
    dstate(4:6) = a;
    dstate(7:9) = [dphi; dtheta; dpsi];
    dstate(10:12) = [dp; dq; dr_ang];
    dstate(13) = m_dot;
end

function dstate = dynamics_coast(t, state, coast, common)
    % COAST PHASE DYNAMICS: No Thrust
    r = state(1:3);
    v = state(4:6);
    m = state(13);
    
    r_norm = norm(r);
    g_acc = -common.mu * r / r_norm^3;
    
    altitude_local = r_norm - common.Re;
    if altitude_local < 100e3
        rho_local = common.rho0 * exp(-altitude_local/common.H);
    else
        rho_local = 0;
    end
    
    v_norm = norm(v);
    if v_norm > 1e-3
        F_drag_mag = 0.5 * coast.Cd * coast.A_ref * rho_local * v_norm^2;
        F_drag = -F_drag_mag * (v/v_norm);
    else
        F_drag = [0;0;0];
    end
    
    T_inertial = [0; 0; 0];
    F_total = T_inertial + F_drag;
    a = F_total/m + g_acc;
    
    dphi = 0; dtheta = 0; dpsi = 0;
    dp = 0; dq = 0; dr_ang = 0;
    m_dot = 0;
    
    dstate = zeros(13,1);
    dstate(1:3) = v;
    dstate(4:6) = a;
    dstate(7:9) = [dphi; dtheta; dpsi];
    dstate(10:12) = [dp; dq; dr_ang];
    dstate(13) = m_dot;
end

function R = euler_to_dcm(phi, theta, psi)
    % Compute the Direction Cosine Matrix (DCM) using a 3-2-1 (yaw-pitch-roll) sequence.
    R_roll = [1, 0, 0;
              0, cos(phi), -sin(phi);
              0, sin(phi), cos(phi)];
    R_pitch = [cos(theta), 0, sin(theta);
               0, 1, 0;
              -sin(theta), 0, cos(theta)];
    R_yaw = [cos(psi), -sin(psi), 0;
             sin(psi), cos(psi), 0;
             0, 0, 1];
    R = R_yaw * R_pitch * R_roll;
end