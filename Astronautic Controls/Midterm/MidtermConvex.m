clear;clc;close all

% Configurable Parameters
N  = 1;         % Number of intervals
t0 = 0;          % Initial time (seconds)
tF = 15;         % Final time (seconds)
T  = (tF - t0) / N;  % Sampling time

a = 6728;            % km
MU = 398600;         % km^3/s^2, gravitational parameter
n = sqrt(MU/a^3);

% Continuous-time A_cw matrix
A = [ 0    0    0    1      0    0;
         0    0    0    0      1    0;
         0    0    0    0      0    1;
         3*n^2  0    0    0    2*n   0;
         0      0    0   -2*n   0    0;
         0      0   -n^2  0     0    0 ];

% Continuous-time B_cw matrix
B = [ 0   0   0;
         0   0   0;
         0   0   0;
         1   0   0;
         0   1   0;
         0   0   1 ];

c = zeros(6,1);

% Discretize the System over One Interval (T seconds)
n = size(A,1);  % State dimension (7)
m = size(B,2);  % Input dimension (4)

sys_c = ss(A, B, eye(n), zeros(n,m));
sys_d = c2d(sys_c, T, 'zoh');
A_d = sys_d.A;
B_d = sys_d.B;

% Compute Overall "Jump" Matrices from t0 to tF (N intervals)
A_jump = A_d^N;   % Overall state transition matrix

B_jump = zeros(n, m);
c_jump = zeros(n, 1);
for i = 0:(N-1)
    B_jump = B_jump + A_d^i * B_d;
    % c_jump = c_jump + A_d^i * c_d;
end

% Display the Results
% disp('Discrete one-step state matrix A_d:');
% disp(A_d)
% disp('Discrete one-step input matrix B_d:');
% disp(B_d)
% disp('Discrete one-step constant vector c_d:');
% disp(c_d)

disp('Overall state transition matrix A_jump (from t0 to tF):');
disp(A_jump)
disp('Overall input matrix B_jump (from t0 to tF):');
disp(B_jump)
disp('Overall constant term contribution c_jump (from t0 to tF):');
disp(c_jump)

%% 
clear; clc; close all;

% Given
mu = 398600;      % [km^3/s^2] Earth's gravitational parameter
a  = 6728;        % [km] Orbit radius (approx. ISS altitude + Earth radius)
n  = sqrt(mu/a^3); % Mean motion

% Continuous-time state-space matrices for CWH dynamics:
Ac = [0    0    0    1    0    0;
      0    0    0    0    1    0;
      0    0    0    0    0    1;
      3*n^2 0    0    0    2*n  0;
      0    0    0   -2*n  0    0;
      0    0   -n^2  0    0    0];

Bc = [zeros(3); eye(3)];

% Discretization using Zero-Order Hold (ZOH)
t0 = 0;            
tf = 3600;           % [s] Total simulation time (e.g., 1 hour)
N  = 1000;           % Number of time steps
dt = (tf - t0)/N;    % Time step size [s]

% Create continuous-time state-space system and discretize:
sysC = ss(Ac, Bc, eye(6), zeros(6,3));
sysD = c2d(sysC, dt);  % Using zero-order hold
Ad = sysD.A;
Bd = sysD.B;

% Clear YALMIP environment
yalmip('clear');

% Problem dimensions
nx = 6;  % Number of states
nu = 3;  % Number of control inputs

% Define decision variables for state trajectory x and control u:
x = sdpvar(nx, N+1); % States over the horizon (including initial and final)
u = sdpvar(nu, N);   % Control inputs

% Initialize constraints and objective
Constraints = [];
Objective = 0;

% Maximum control acceleration
umax = 5e-6;   

% Initial condition
x0 = [-1.5; -0.5; 0.2; 5E-3; 0; 0];
Constraints = [Constraints, x(:,1) == x0];

% Final state target - Rendezvous to the chief
x_target = [-0.25; 0; 0; 0.2E-3; 0; 0];
Constraints = [Constraints, x(:,N+1) == x_target];
rmax = 2; % km
vmax = 0.8E-3; % km/s

% Loop through the horizon to impose dynamics, control constraints, and cost
for k = 1:N

    % Position state path constraint
    Constraints = [Constraints, norm(x(1:3,k),2) <= rmax];

    if k >= floor(N/2) + 1
        % Velocity State Path Constraint
        Constraints = [Constraints, norm(x(5,k),2) <= vmax];
    end
    
    % Dynamics: x_{k+1} = Ad*x_k + Bd*u_k
    Constraints = [Constraints, x(:,k+1) == Ad*x(:,k) + Bd*u(:,k)];
    
    % Control constraint: ||u_k||_2 <= umax
    Constraints = [Constraints, norm(u(:,k), 2) <= umax];
    
    % Cost - Sum the Euclidean norm of the control inputs (integrated over time)
    Objective = Objective + norm(u(:,k), 2)*dt;
end

% Set MOSEK options
ops = sdpsettings('solver','mosek','verbose',1);

% Solve
sol = optimize(Constraints, Objective, ops);

if sol.problem == 0
    disp('Optimization is feasible and converged.');
    x_opt = value(x);
    u_opt = value(u);
    J_opt = value(Objective);
    fprintf('Optimal cost J = %.6g\n', J_opt);
else
    error('Optimization failed: %s', sol.info);
end

time = linspace(t0, tf, N+1)/60;

% Plot relative position over time
figure('Color','white','Position',[100 100 1200 800]); 
subplot(3,1,1); plot(time, x_opt(1,:), '-o'); ylabel('x [km]'); grid on;
axis tight
subplot(3,1,2); plot(time, x_opt(2,:), '-o'); ylabel('y [km]'); grid on;
axis tight
subplot(3,1,3); plot(time, x_opt(3,:), '-o'); ylabel('z [km]'); grid on;
xlabel('Time [Minutes]');
sgtitle('Relative Position vs. Time','FontSize',30);
axis tight

% Plot relative velocity over time
figure('Color','white','Position',[150 150 1200 800]); 
subplot(3,1,1); plot(time, x_opt(4,:), '-o'); ylabel('vx [km/s]'); grid on;
axis tight
subplot(3,1,2); plot(time, x_opt(5,:), '-o'); ylabel('vy [km/s]'); grid on;
axis tight
subplot(3,1,3); plot(time, x_opt(6,:), '-o'); ylabel('vz [km/s]'); grid on;
xlabel('Time [Minutes]');
sgtitle('Relative Velocity vs. Time','FontSize',30);
axis tight

% Plot control inputs over time
timeu = linspace(t0, tf, N)/60;
figure('Color','white','Position',[200 200 1200 800]); 
plot(timeu, u_opt(1,:), '-o', 'LineWidth',1.2); hold on;
plot(timeu, u_opt(2,:), '-o', 'LineWidth',1.2);
plot(timeu, u_opt(3,:), '-o', 'LineWidth',1.2);
plot(timeu,vecnorm(u_opt,2), '-o', 'LineWidth',1.2);
grid on;
xlabel('Time [Minutes]');
ylabel('Control [km/s^2]');
legend('u_x','u_y','u_z','u_{norm}','Location','best');
title('Control Acceleration vs. Time','FontSize',30);
axis tight

% Define an interpolation for the computed control history (u_opt)
u_interp = @(t) [interp1(timeu*60, u_opt(1,:), t, 'previous', 'extrap');
                 interp1(timeu*60, u_opt(2,:), t, 'previous', 'extrap');
                 interp1(timeu*60, u_opt(3,:), t, 'previous', 'extrap')];

% Continuous CWH dynamics
dynamics_CWH = @(t, x) Ac*x + Bc*u_interp(t);

% Integrate
[t_CWH, x_CWH] = ode45(dynamics_CWH, [t0 tf], x0);

% Convert time to minutes for plotting
time_CWH = t_CWH/60;

% Chief's circular orbit (at t=0):
r_c0 = [a; 0; 0];         % position in km
v_c0 = [0; n*a; 0];        % velocity in km/s

% The initial relative state in LVLH from optimization:
x_rel0 = x0(1:3);          % relative position [km]
v_rel0 = x0(4:6);          % relative velocity [km/s]

dotQ0 = [-n*x_rel0(2); n*x_rel0(1); 0];

% Thus the deputy's absolute initial state is:
r_d0 = r_c0 + x_rel0;              % deputy position in km
v_d0 = v_c0 + dotQ0 + v_rel0;      % deputy velocity in km/s
X0 = [r_d0; v_d0];                 % combined state vector [r_d; v_d]


nonlinearDynamics = @(t, X) nonlinDynamics(t, X, u_interp, mu, a, n);

% Integrate
[t_nl, X_nl] = ode45(nonlinearDynamics, [t0 tf], X0);

% Convert time to minutes for plotting later
time_nl = t_nl/60;

% Convert the nonlinear simulation results back to LVLH (relative) coordinates
% For each time step, compute the chief's state, then obtain the relative state
x_nl = zeros(length(t_nl), 6);
for i = 1:length(t_nl)
    t_val = t_nl(i);
    r_d = X_nl(i,1:3)';  % deputy position (inertial)
    v_d = X_nl(i,4:6)';  % deputy velocity (inertial)
    
    % Chief state at time t_val (circular orbit)
    r_c = a * [cos(n*t_val); sin(n*t_val); 0];
    v_c = a * [-n*sin(n*t_val); n*cos(n*t_val); 0];
    
    % LVLH rotation matrix (from LVLH to inertial)
    Q = [cos(n*t_val), -sin(n*t_val), 0; 
         sin(n*t_val),  cos(n*t_val), 0;
         0,             0,            1];
    
    % Relative position in LVLH
    x_rel = Q'*(r_d - r_c);
    
    % Relative velocity:
    Omega = [0, -n, 0; n, 0, 0; 0, 0, 0];
    v_rel = Q'*(v_d - v_c) - Omega*x_rel;
    
    % Store the relative state
    x_nl(i,:) = [x_rel', v_rel'];
end

% Relative Position vs Time
figure('Color','white','Position',[100 100 1200 800]); 
subplot(3,1,1); 
plot(time, x_opt(1,:), 'b', 'LineWidth',1.5); hold on;
% plot(time_CWH, x_CWH(:,1), 'b', 'LineWidth',1.5);
plot(time_nl, x_nl(:,1), 'r', 'LineWidth',1.5);
% ylabel('x [km]'); grid on; legend('Optimized','CWH','Nonlinear','Location','best');
ylabel('x [km]'); grid on; legend('Optimized','Nonlinear','Location','best');
subplot(3,1,2); 
plot(time, x_opt(2,:), 'b', 'LineWidth',1.5); hold on;
% plot(time_CWH, x_CWH(:,2), 'b', 'LineWidth',1.5);
plot(time_nl, x_nl(:,2), 'r', 'LineWidth',1.5);
ylabel('y [km]'); grid on;
subplot(3,1,3); 
plot(time, x_opt(3,:), 'b', 'LineWidth',1.5); hold on;
% plot(time_CWH, x_CWH(:,3), 'b', 'LineWidth',1.5);
plot(time_nl, x_nl(:,3), 'r', 'LineWidth',1.5);
ylabel('z [km]'); grid on; xlabel('Time [Minutes]');
sgtitle('Relative Position vs. Time','FontSize',25);

% Figure: Relative Velocity vs Time (vx, vy, vz)
figure('Color','white','Position',[150 150 1200 800]); 
subplot(3,1,1); 
plot(time, x_opt(4,:), 'b', 'LineWidth',1.5); hold on;
% plot(time_CWH, x_CWH(:,4), 'b', 'LineWidth',1.5);
plot(time_nl, x_nl(:,4), 'r', 'LineWidth',1.5);
% ylabel('v_x [km/s]'); grid on; legend('Optimized','CWH','Nonlinear','Location','best');
ylabel('v_x [km/s]'); grid on; legend('Optimized','Nonlinear','Location','best');
subplot(3,1,2); 
plot(time, x_opt(5,:), 'b', 'LineWidth',1.5); hold on;
% plot(time_CWH, x_CWH(:,5), 'b', 'LineWidth',1.5);
plot(time_nl, x_nl(:,5), 'r', 'LineWidth',1.5);
ylabel('v_y [km/s]'); grid on;
subplot(3,1,3); 
plot(time, x_opt(6,:), 'b', 'LineWidth',1.5); hold on;
% plot(time_CWH, x_CWH(:,6), 'b', 'LineWidth',1.5);
plot(time_nl, x_nl(:,6), 'r', 'LineWidth',1.5);
ylabel('v_z [km/s]'); grid on; xlabel('Time [Minutes]');
sgtitle('Relative Velocity vs. Time','FontSize',25);


function dXdt = nonlinDynamics(t, X, u_interp, mu, a, n)
    % X is 6x1 vector: [r_d; v_d] in inertial coordinates
    r_d = X(1:3);
    v_d = X(4:6);
    
    % Compute the chief's absolute state at time t
    r_c = a * [cos(n*t); sin(n*t); 0];
    v_c = a * [-n*sin(n*t); n*cos(n*t); 0];
    
    % LVLH rotation matrix from LVLH to inertial
    Q = [cos(n*t), -sin(n*t), 0;
         sin(n*t),  cos(n*t), 0;
         0,         0,        1];
    
    % Compute control in inertial coordinates
    u = Q * u_interp(t);
    
    % Two-body dynamics for the deputy
    drdt = v_d;
    dvdt = -mu * r_d / norm(r_d)^3 + u;
    
    dXdt = [drdt; dvdt];
end