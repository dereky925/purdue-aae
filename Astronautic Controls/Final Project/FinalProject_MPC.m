% ==================== Initialization ====================
clear;          % Remove all variables from workspace
clc;            % Clear the command window
close all;      % Close all figure windows
tic             % Start timer to measure total run time

% ==================== Problem Parameters (Earth Model) ====================
alpha           = 0.0;           % Fuel consumption coefficient [s/km] (used in MPC linearized mass model)
test_alpha      = 0.3;           % Alternative fuel‐flow coefficient for theoretical metrics
g               = [-9.81e-3;     % Gravity vector in km/s^2 (only along r1 axis)
                    0;
                    0];
T_min           = 1;             % Minimum throttle scalar [kg·km/s^2]
T_max           = 500;           % Maximum throttle scalar [kg·km/s^2]
gamma_min_deg   = 1;             % Minimum glide‐slope angle [degrees]
m0              = 2000;          % Initial interceptor mass [kg]
z0              = log(m0);       % Log‐mass state (linearizes mass dynamics)

% ==================== Interceptor MPC Discretization ====================
N               = 30;            % Number of control intervals
r0              = [400;1000;1000];  % Initial interceptor position [km] in body‐fixed frame
v0              = [0;-4;-4];        % Initial interceptor velocity [km/s]
gamma_min_rad   = gamma_min_deg * pi/180;  % Glide‐slope constraint in radians
t_delay         = 20;            % Configurable delay before interceptor launch [s]

% ==================== ICBM Boost‐Phase Simulation ====================
t_boost         = 180;           % Total boost‐phase duration [s]
N_icbm          = 300;           % Number of time samples for boost simulation
pitch_max       = 45 * pi/180;   % Maximum pitch‐over angle (from vertical) [rad]
a_mag           = 0.02;          % Thrust acceleration magnitude [km/s^2]
dt_icbm         = t_boost / (N_icbm-1);  % Time step for ICBM simulation [s]

% Preallocate arrays for ICBM trajectory
t_icbm = linspace(0, t_boost, N_icbm);      % Time vector
r_icbm = zeros(3, N_icbm);                  % Position states
v_icbm = zeros(3, N_icbm);                  % Velocity states

% Simulate curved pitch‐program during boost
for i = 2:N_icbm
    tau         = t_icbm(i);                     % Current time [s]
    pitch       = pitch_max * (tau / t_boost);   % Linear ramp from 0 to pitch_max
    a_thrust    = [a_mag*cos(pitch);              % Thrust acceleration vector
                   a_mag*sin(pitch);
                   0];
    a_total     = a_thrust + g;                  % Total acceleration = thrust + gravity
    v_icbm(:,i) = v_icbm(:,i-1) + a_total*dt_icbm;   % Integrate velocity
    r_icbm(:,i) = r_icbm(:,i-1) + v_icbm(:,i-1)*dt_icbm; % Integrate position
end
% Reference for realistic boost‐phase pitch program: Wikipedia “Intercontinental ballistic missile”

% ==================== Interceptor Continuous‐Time Dynamics ====================
% States: [r; v; z] where z = log(mass)
A_ct = zeros(7);        % Continuous‐time state matrix
A_ct(1,4) = 1;          % dr/dt = v
A_ct(2,5) = 1;
A_ct(3,6) = 1;

B_ct = zeros(7,4);      % Continuous‐time input matrix
B_ct(4,1) = 1;          % dv/dt = a_control
B_ct(5,2) = 1;
B_ct(6,3) = 1;
B_ct(7,4) = -alpha;     % dz/dt = -alpha * throttle

c_ct = zeros(7,1);      % Continuous‐time affine term
c_ct(4:6) = g;          % Gravity enters acceleration dynamics

% ==================== Solver Setup & Initial Feasibility Search ====================
ops       = sdpsettings('solver','mosek','verbose',0);
tol       = 1e-3;        % Bisection convergence tolerance on burn time [s]
maxIter   = 20;          % Max bisection iterations
T_low     = 0;           % Lower bound on interceptor burn time
T_high    = N * 1.0;     % Initial upper bound guess (N seconds)

% Find a feasible upper bound T_high by doubling until constraints are met
while true
    dt_test = T_high / N;  
    [Ad0,Bd0,cd0] = zohOneStep(A_ct, B_ct, c_ct, dt_test);
    [~,~,cons0] = createDecisionVarsAndConstraints(Ad0, Bd0, cd0, dt_test, ...
                          r0, v0, z0, m0, alpha, T_min, T_max, gamma_min_rad, N, [0;0;0]);
    sol = optimize(cons0, [], ops);
    if sol.problem == 0
        break;
    end
    T_high = 2 * T_high;
    if T_high > 1e4
        error('Could not find a feasible upper bound on burn time.');
    end
end

% ==================== Bisection to Find Optimal Burn Time ====================
for iter = 1:maxIter
    % Midpoint candidate
    T_mid = 0.5 * (T_low + T_high);
    dt_mid = T_mid / N;
    
    % Discretize and build constraints
    [Ad1, Bd1, cd1] = zohOneStep(A_ct, B_ct, c_ct, dt_mid);
    [~,~,cons1] = createDecisionVarsAndConstraints(Ad1, Bd1, cd1, dt_mid, ...
        r0, v0, z0, m0, alpha, T_min, T_max, gamma_min_rad, N, [0;0;0]);
    
    % Test feasibility
    sol = optimize(cons1, [], ops);
    if sol.problem == 0
        % Feasible → can reduce upper bound
        T_high  = T_mid;
        bestAd  = Ad1;
        bestBd  = Bd1;
        bestcd  = cd1;
    else
        % Infeasible → increase lower bound
        T_low   = T_mid;
    end
    
    % Check convergence
    if (T_high - T_low) < tol
        break;
    end
end

% Optimal burn time and derived quantities
T_opt         = T_high;             % Optimal interceptor burn time [s]
dt            = T_opt / N;          % Time step for MPC
t_int_global  = t_delay + T_opt;    % Total time from ICBM launch to intercept [s]

% ==================== Compute Intercept Point ====================
% Interpolate ICBM trajectory to find where it will be at t_int_global
r_int_pt = interp1(t_icbm', r_icbm', t_int_global)';

% Solve final MPC targeting that intercept point
[X, U, consOpt] = createDecisionVarsAndConstraints(bestAd, bestBd, bestcd, dt, ...
    r0, v0, z0, m0, alpha, T_min, T_max, gamma_min_rad, N, r_int_pt);
sol = optimize(consOpt, [], ops);
if sol.problem ~= 0
    error('No feasible solution at T_opt = %.4f s', T_opt);
end

toc  % End timer

% ==================== Display Key Results ====================
fprintf('Interceptor flight time (since interceptor launch): %.2f s\n', T_opt);
fprintf('Global intercept time (since ICBM launch):         %.2f s\n\n', t_int_global);

% ==================== Extract Throttle History ====================
s_vals = zeros(1, N);        % Preallocate
for k = 1:N
    s_vals(k) = value(U{k}(4));   % Extract throttle scalar σ_k
end

% === Real Fuel Mass via Specific Impulse Conversion ===
Isp    = 300;                % Engine specific impulse [s]
g0     = 9.81;               % Standard gravity [m/s^2]
Tn     = s_vals * 1e3;       % Convert σ from [kg·km/s^2] → thrust [N]
mdot   = Tn ./ (Isp * g0);   % Mass‐flow rate [kg/s]
fuel_kg = sum(mdot) * dt;    % Total fuel consumed [kg]
fprintf('Realistic fuel consumed: %.2f kg  (Isp=%g s)\n\n', fuel_kg, Isp);

% ==================== Extract and Analyze Trajectory ====================
% Pull full state trajectory from solution
Xval = zeros(7, N+1);
for k = 1:N+1
    Xval(:,k) = value(X{k});
end

% Compute speeds and stagnation temperature
vVals   = Xval(4:6, :);           % Velocity history [km/s]
v_int   = sqrt(sum(vVals.^2, 1)); % Speed magnitude [km/s]
v_peak  = max(v_int);             % Peak speed
T_ambient = 288;                  % Ambient temperature [K]
cp = 1004;                        % Specific heat of air [J/(kg·K)]
v_m_s = v_peak * 1e3;             % Convert km/s → m/s
T_stag = T_ambient + v_m_s^2/(2*cp);
fprintf('Peak stagnation temperature: %.2f K (v_{max}=%.2f km/s)\n\n', T_stag, v_peak);

% ==================== Plotting ====================
% Find index on ICBM profile corresponding to intercept time
[~, idx_int] = min(abs(t_icbm - t_int_global));

% 1) 3D Trajectories: interceptor vs. ICBM
figure('Color','white','Position',[100,100,900,700]);
plot3(-Xval(2,:), Xval(3,:), Xval(1,:), 'b-o', 'LineWidth',2); hold on;
plot3(-r_icbm(2,1:idx_int), r_icbm(3,1:idx_int), r_icbm(1,1:idx_int), 'r','LineWidth',2);
plot3(-r_icbm(2,idx_int:end), r_icbm(3,idx_int:end), r_icbm(1,idx_int:end), ...
      'Color',[0.5 0.5 0.5],'LineWidth',2);
plot3(-r_int_pt(2), r_int_pt(3), r_int_pt(1), 'kp','MarkerSize',25,'MarkerFaceColor','y');
grid on; axis equal;
xlabel('r_2 [km]','FontSize',20); ylabel('r_3 [km]','FontSize',20); zlabel('r_1 [km]','FontSize',16);
title('3D Trajectories: Interceptor vs. ICBM','FontSize',20);
legend('Interceptor','ICBM (boost)','ICBM (post‐intercept)','Intercept Point','Location','best');

% 2) ICBM Position vs Time
figure('Color','white','Position',[200,200,800,500]);
plot(t_icbm, r_icbm(1,:), 'r-', t_icbm, r_icbm(2,:), 'b--', t_icbm, r_icbm(3,:), 'g-.', 'LineWidth',2);
grid on;
xlabel('Time [s]','FontSize',20); ylabel('Position [km]','FontSize',20);
legend('r_1','r_2','r_3'); title('ICBM Position vs. Time','FontSize',20);

% 3) ICBM Velocity vs Time
figure('Color','white','Position',[200,200,800,500]);
plot(t_icbm, v_icbm(1,:), 'r-', t_icbm, v_icbm(2,:), 'b--', t_icbm, v_icbm(3,:), 'g-.', 'LineWidth',2);
grid on;
xlabel('Time [s]','FontSize',20); ylabel('Velocity [km/s]','FontSize',20);
legend('v_1','v_2','v_3'); title('ICBM Velocity vs. Time','FontSize',20);

% 4) Interceptor diagnostics (lossless property, state & control histories)
plotLosslessProperty(U, N, dt);
plotResults(X, U, N, dt, t_delay, t_int_global);

% ==================== Local Function Definitions ====================

%------------------------------------------------------------------------------
function [X, U, cons] = createDecisionVarsAndConstraints(Ad, Bd, cd, dt, ...
    r0, v0, z0, m0, alpha, T_min, T_max, gamma, N, r_target)
% Builds optimization variables and constraints for the interceptor MPC
%
% Inputs:
%   Ad, Bd, cd: discrete-time dynamics (A_d, B_d, c_d)
%   dt         : time step
%   r0, v0, z0 : initial states
%   m0, alpha  : mass parameters
%   T_min, T_max: throttle limits
%   gamma      : minimum glide-slope angle (rad)
%   N          : number of intervals
%   r_target   : terminal intercept point (3×1)
%
% Outputs:
%   X : cell array of state sdpvars (7×1 each)
%   U : cell array of control sdpvars (4×1 each: [a; σ])
%   cons : YALMIP constraint array

    % Create state and input variables
    X = cell(N+1, 1);
    U = cell(N,   1);
    for k = 1:N+1
        X{k} = sdpvar(7,1,'full');   % States: [r; v; z]
    end
    for k = 1:N
        U{k} = sdpvar(4,1,'full');   % Inputs: [a_x;a_y;a_z;σ]
    end

    % Initial condition
    cons = (X{1} == [r0; v0; z0]);

    % Dynamics constraints
    for k = 1:N
        cons = [cons, X{k+1} == Ad*X{k} + Bd*U{k} + cd];
    end

    % Precompute lower bound on log-mass for linearization
    z_lb = zeros(N+1,1);
    for i = 0:N
        z_lb(i+1) = log(m0 - alpha * T_max * (i * dt));
    end

    % Throttle and mass‐consistency constraints
    for k = 1:N
        zk      = X{k}(7);          % log-mass state
        ak      = U{k}(1:3);        % acceleration command
        sk      = U{k}(4);          % throttle scalar
        e_lbk   = exp(-z_lb(k));    % linearization factor
        lin_app = e_lbk * (1 - (zk - z_lb(k)));  % first‐order approx of exp(-z)

        % Enforce T_min <= σ <= linearized T_max, and norm(a) ≤ σ
        cons = [cons,
            sk >= T_min * e_lbk,
            sk <= T_max * lin_app,
            norm(ak,2) <= sk
        ];
    end

    % Terminal constraint: must hit the intercept point exactly
    cons = [cons, X{end}(1:3) == r_target];

    % Glide‐slope constraint at every step: r1 ≥ tan(γ)*||[r2;r3]||
    for k = 1:N+1
        rk = X{k}(1:3);
        cons = [cons, rk(1) >= tan(gamma) * norm(rk(2:3),2)];
    end
end

%------------------------------------------------------------------------------
function [Ad, Bd, cd] = zohOneStep(A, B, c, dt)
% Zero-Order Hold discretization for continuous system ẋ = A x + B u + c
%
% Outputs discrete-time:
%   x[k+1] = Ad x[k] + Bd u[k] + cd

    n = size(A,1);
    M = expm(A * dt);       % Matrix exponential

    % Compute integral of exp(A τ) dτ
    if rank(A) == n
        % Analytical when A invertible
        IntA = A \ (M - eye(n));
    else
        % Numerical integration fallback
        steps = 50;
        dTau  = dt / steps;
        IntA  = zeros(n);
        for i = 1:steps
            tau = (i-1) * dTau;
            IntA = IntA + expm(A * tau) * dTau;
        end
    end

    Ad = M;
    Bd = IntA * B;
    cd = IntA * c;
end

%------------------------------------------------------------------------------
function plotLosslessProperty(U, N, dt)
% Plots σ_k - ||a_k|| to verify lossless property (throttle covers control norm)
    diffVals = zeros(1,N);
    for k = 1:N
        diffVals(k) = value(U{k}(4)) - norm(value(U{k}(1:3)),2);
    end
    t = linspace(0, dt*(N-1), N);
    figure('Color','white','Position',[0,100,1000,600]);
    plot(t, diffVals, 'k-o', 'LineWidth',2);
    grid on;
    xlabel('Time [s]','FontSize',20);
    ylabel('\sigma_k - ||a_k||_2','FontSize',20);
    title('Lossless Property vs. Time','FontSize',20);
end

%------------------------------------------------------------------------------
function plotResults(X, U, N, dt, t0, tF)
% Plots interceptor state and control histories over [t0, tF]
%
% States: r vs time, v vs time
% Controls: throttle percentage, raw accelerations

    % Extract numeric states & controls
    Xv = zeros(7, N+1);
    for k = 1:N+1
        Xv(:,k) = value(X{k});
    end
    r = Xv(1:3, :);
    v = Xv(4:6, :);

    % Compute accelerations and throttle %
    a = diff(v,1,2) / dt;           % finite difference
    s = zeros(1,N);
    for k = 1:N
        s(k) = value(U{k}(4));
    end
    throttlePct = (s / 500) * 100;  % example conversion to %

    % Time vectors
    t_state = linspace(t0, tF, N+1);
    t_ctrl  = linspace(t0, tF - dt, N);

    % Position vs time
    figure('Color','white','Position',[500,0,1000,500]);
    plot(t_state, r(1,:), 'r-o', t_state, r(2,:), 'b-o', t_state, r(3,:), 'g-o','LineWidth',2);
    grid on;
    xlabel('Time [s]','FontSize',20);
    ylabel('Position [km]','FontSize',20);
    legend('r_1','r_2','r_3');
    title('Interceptor Positions vs. Time','FontSize',20);

    % Velocity vs time
    figure('Color','white','Position',[3000,1000,1000,500]);
    plot(t_state, v(1,:), 'r-o', t_state, v(2,:), 'b-o', t_state, v(3,:), 'g-o','LineWidth',2);
    grid on;
    xlabel('Time [s]','FontSize',20);
    ylabel('Velocity [km/s]','FontSize',20);
    legend('v_1','v_2','v_3');
    title('Interceptor Velocities vs. Time','FontSize',20);

    % Throttle % vs time
    figure('Color','white','Position',[3000,0,1000,500]);
    plot(t_ctrl, throttlePct, 'm-o','LineWidth',2);
    grid on;
    xlabel('Time [s]','FontSize',20);
    ylabel('Throttle [%]','FontSize',20);
    title('Interceptor Throttle vs. Time','FontSize',20);

    % Control inputs (accelerations) vs time
    figure('Color','white','Position',[3000,1500,1000,500]);
    plot(t_ctrl, a(1,:), 'r-', t_ctrl, a(2,:), 'g-', t_ctrl, a(3,:), 'b-', ...
         t_ctrl, sqrt(sum(a.^2,1)), 'k--', 'LineWidth',2);
    grid on;
    xlabel('Time [s]','FontSize',20);
    ylabel('Control [km/s^2]','FontSize',20);
    legend('a_1','a_2','a_3','||a||','Location','best');
    title('Interceptor Control Inputs','FontSize',20);
end