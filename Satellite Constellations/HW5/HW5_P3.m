
clear; close all; clc;
mu = 398600;              % Earth's gravitational parameter [km^3/s^2]
a = 7000;                 % semi-major axis [km]
e = 0;                    % eccentricity (circular orbit)
i_deg = 60;               % inclination in degrees
i = deg2rad(i_deg);
Omega_ref_deg = 0;        % RAAN for reference satellite in degrees
Omega_ref = deg2rad(Omega_ref_deg);
Omega_target_deg = 0.1;   % RAAN for target satellite in degrees
Omega_target = deg2rad(Omega_target_deg);
w_deg = 0;                % argument of perigee (irrelevant for circular orbit)
w = deg2rad(w_deg);
theta0_ref_deg = 0;       % initial true anomaly for reference [deg]
theta0_ref = deg2rad(theta0_ref_deg);
theta0_target_deg = 0.1;  % initial true anomaly for target [deg]
theta0_target = deg2rad(theta0_target_deg);

n = sqrt(mu/a^3);         % mean motion [rad/s]
T = 2*pi/n;               % orbital period [s]

N = 10000;
t = linspace(0, T, N);

theta_ref = theta0_ref + n*t;
theta_target = theta0_target + n*t;

% Define rotation matrix functions (PQW to inertial)
R3 = @(angle) [cos(angle) sin(angle) 0; -sin(angle) cos(angle) 0; 0 0 1];
R1 = @(angle) [1 0 0; 0 cos(angle) sin(angle); 0 -sin(angle) cos(angle)];

% Compute inertial positions and velocities using orbital elements
R_ref = R3(-Omega_ref) * R1(-i) * R3(-w);
r_ref_PQW = [a*cos(theta_ref); a*sin(theta_ref); zeros(1, length(theta_ref))];
v_ref_PQW = sqrt(mu/a) * [-sin(theta_ref); cos(theta_ref); zeros(1, length(theta_ref))];
r_ref_inertial = R_ref * r_ref_PQW;
v_ref_inertial = R_ref * v_ref_PQW;

% For the target satellite
R_target = R3(-Omega_target) * R1(-i) * R3(-w);
r_target_PQW = [a*cos(theta_target); a*sin(theta_target); zeros(1, length(theta_target))];
v_target_PQW = sqrt(mu/a) * [-sin(theta_target); cos(theta_target); zeros(1, length(theta_target))];
r_target_inertial = R_target * r_target_PQW;
v_target_inertial = R_target * v_target_PQW;

% (a) Plot the relative position in the inertial frame over one revolution
r_rel_inertial = r_target_inertial - r_ref_inertial;
figure('Position',[100,100,1000,1000],'Color','w');
plot3(r_rel_inertial(1,:), r_rel_inertial(2,:), r_rel_inertial(3,:), 'LineWidth',2);
grid on; axis equal;
xlabel('X [km]','FontSize',20);
ylabel('Y [km]','FontSize',20);
zlabel('Z [km]','FontSize',20);
title('Relative Position in Inertial Frame','FontSize',20);

% (b) Compute and plot the relative position in the LVLH frame
r_rel_LVLH = zeros(3, N);
R_LVLH = zeros(3, N);
S_LVLH = zeros(3, N);
W_LVLH = zeros(3, N);

for k = 1:N
    r_ref_k = r_ref_inertial(:,k);
    v_ref_k = v_ref_inertial(:,k);
    Rhat = r_ref_k / norm(r_ref_k);
    What = cross(r_ref_k, v_ref_k); 
    What = What / norm(What);
    Shat = cross(What, Rhat);
    
    % Save the LVLH basis (used later for transforming the CW solution)
    R_LVLH(:,k) = Rhat;
    S_LVLH(:,k) = Shat;
    W_LVLH(:,k) = What;
    
    % Express the relative inertial position vector in the LVLH frame
    r_rel = r_rel_inertial(:,k);
    r_rel_LVLH(:,k) = [dot(r_rel, Rhat); dot(r_rel, Shat); dot(r_rel, What)];
end

figure('Position',[100,100,1000,1000],'Color','w');
plot3(r_rel_LVLH(1,:), r_rel_LVLH(2,:), r_rel_LVLH(3,:), 'LineWidth',2);
grid on; axis tight;
xlabel('Radial (R) [km]','FontSize',20);
ylabel('Along-track (S) [km]','FontSize',20);
zlabel('Cross-track (W) [km]','FontSize',20);
title('Relative Position in LVLH Frame','FontSize',20);

% (c) Determine the initial conditions (relative position and velocity) in the LVLH frame
r0_inertial = r_rel_inertial(:,1);
v0_inertial = v_target_inertial(:,1) - v_ref_inertial(:,1);
R0 = R_LVLH(:,1);
S0 = S_LVLH(:,1);
W0 = W_LVLH(:,1);
% Angular velocity of the LVLH frame for a circular orbit 
omega0 = n * W0;  % [rad/s]

% Compute the relative position in LVLH (simple rotation)
r0_LVLH = [dot(r0_inertial, R0); dot(r0_inertial, S0); dot(r0_inertial, W0)];

% Apply the transport theorem to get the LVLH relative velocity
v0_LVLH = [ dot(v0_inertial - cross(omega0, r0_inertial), R0);
            dot(v0_inertial - cross(omega0, r0_inertial), S0);
            dot(v0_inertial - cross(omega0, r0_inertial), W0) ];

fprintf('Initial relative position in LVLH [km]: [%f, %f, %f]\n', r0_LVLH);
fprintf('Initial relative velocity in LVLH [km/s]: [%.10f, %.10f, %.10f]\n', v0_LVLH);

% (d) CW solution in LVLH frame and its difference with the Keplerian solution
x0 = r0_LVLH(1); 
y0 = r0_LVLH(2); 
z0 = r0_LVLH(3);
vx0 = v0_LVLH(1); 
vy0 = v0_LVLH(2); 
vz0 = v0_LVLH(3);

x_CW = (4 - 3*cos(n*t))*x0 + (sin(n*t)/n)*vx0 + (2*(1-cos(n*t))/n)*vy0;
y_CW = 6*(sin(n*t) - n*t)*x0 + y0 - (2*(1-cos(n*t))/n)*vx0 + (sin(n*t)/n)*vy0;
z_CW = z0*cos(n*t) + (sin(n*t)/n)*vz0;
r_CW_LVLH = [x_CW; y_CW; z_CW];


diff_LVLH = r_CW_LVLH - r_rel_LVLH;

figure('Position',[100,100,1000,1000],'Color','w');
subplot(3,1,1);
plot(t/3600, diff_LVLH(1,:), 'LineWidth',2);
ylabel('Radial Error [km]','FontSize',20);
grid on;
subplot(3,1,2);
plot(t/3600, diff_LVLH(2,:), 'LineWidth',2);
ylabel('Along-track Error [km]','FontSize',20);
grid on;
subplot(3,1,3);
plot(t/3600, diff_LVLH(3,:), 'LineWidth',2);
xlabel('Time [Hours]','FontSize',20);
ylabel('Cross Error [km]','FontSize',20);
grid on;
sgtitle('CW - Keplerian (LVLH) vs Time','FontSize',20);

% (e) Transform the CW solution from LVLH to the inertial frame
r_CW_inertial = zeros(3, N);
diff_inertial = zeros(3, N);
for k = 1:N
    Rhat = R_LVLH(:,k);
    Shat = S_LVLH(:,k);
    What = W_LVLH(:,k);
    % Transform the CW solution to the inertial frame
    r_CW_inertial(:,k) = x_CW(k)*Rhat + y_CW(k)*Shat + z_CW(k)*What;


    diff_inertial(:,k) = r_CW_inertial(:,k) - r_rel_inertial(:,k);
end

figure('Position',[100,100,1000,1000],'Color','w');
subplot(3,1,1);
plot(t/3600, diff_inertial(1,:), 'LineWidth',2);
ylabel('Diff X [km]','FontSize',20);
grid on;
subplot(3,1,2);
plot(t/3600, diff_inertial(2,:), 'LineWidth',2);
ylabel('Diff Y [km]','FontSize',20);
grid on;
subplot(3,1,3);
plot(t/3600, diff_inertial(3,:), 'LineWidth',2);
xlabel('Time [hours]','FontSize',20);
ylabel('Diff Z [km]','FontSize',20);
grid on;
sgtitle('CW - Keplerian (Inertial) vs. Time','FontSize',20);

% (f) Propagate over 100 orbital revolutions and analyze the error in the inertial frame
N_long = 100 * 10000;      % 100 revolutions with 10000 points per orbit
t_long = linspace(0, 100*T, N_long);

% Full Keplerian propagation for both satellites
theta_ref_long = theta0_ref + n*t_long;
theta_target_long = theta0_target + n*t_long;
r_ref_PQW_long = [a*cos(theta_ref_long); a*sin(theta_ref_long); zeros(1, length(t_long))];
r_target_PQW_long = [a*cos(theta_target_long); a*sin(theta_target_long); zeros(1, length(t_long))];
r_ref_inertial_long = R_ref * r_ref_PQW_long;
r_target_inertial_long = R_target * r_target_PQW_long;
r_rel_inertial_long = r_target_inertial_long - r_ref_inertial_long;

% CW propagation (using the same initial LVLH conditions) over 100 orbits
x_CW_long = (4 - 3*cos(n*t_long))*x0 + (sin(n*t_long)/n)*vx0 + (2*(1-cos(n*t_long))/n)*vy0;
y_CW_long = 6*(sin(n*t_long) - n*t_long)*x0 + y0 - (2*(1-cos(n*t_long))/n)*vx0 + (sin(n*t_long)/n)*vy0;
z_CW_long = z0*cos(n*t_long) + (sin(n*t_long)/n)*vz0;

% Transform the CW solution to inertial coordinates
r_CW_inertial_long = zeros(3, length(t_long));
for k = 1:length(t_long)
    % Recompute the LVLH basis for the reference at time t_long(k)
    r_ref_k = r_ref_inertial_long(:,k);
    v_ref_k = R_ref * (sqrt(mu/a) * [-sin(theta_ref_long(k)); cos(theta_ref_long(k)); 0]);
    Rhat = r_ref_k / norm(r_ref_k);
    What = cross(r_ref_k, v_ref_k); What = What / norm(What);
    Shat = cross(What, Rhat);
    
    r_CW_inertial_long(:,k) = x_CW_long(k)*Rhat + y_CW_long(k)*Shat + z_CW_long(k)*What;
end

% Compute the error in the relative distance (norm difference)
rel_distance_full = vecnorm(r_rel_inertial_long);
rel_distance_CW = vecnorm(r_CW_inertial_long);
error_distance = abs(rel_distance_full - rel_distance_CW);

figure('Position',[100,100,1000,1000],'Color','w');
plot(t_long/(3600*24), error_distance, 'LineWidth',2);
xlabel('Time [days]','FontSize',20);
ylabel('Error in Relative Distance [km]','FontSize',20);
title('Error in Relative Distance over 100 Orbits','FontSize',20);
grid on;