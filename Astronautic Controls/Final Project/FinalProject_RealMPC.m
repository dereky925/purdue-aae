% ==================== Initialization ====================
clear; clc; close all; tic;
set(groot,'defaultAxesFontSize',20,'defaultTextFontSize',20);

% ==================== Problem Parameters (Earth Model) ====
alpha           = 0.0;            % Fuel consumption coefficient [s/km]
test_alpha      = 0.3;            % Alternative fuel-flow coefficient for theoretical metrics
g               = [-9.81e-3; 0; 0];% Gravity vector in km/s^2 (along r1 axis)
T_min           = 1;              % Minimum throttle scalar [kg·km/s^2]
T_max           = 500;            % Maximum throttle scalar [kg·km/s^2]
gamma_min_deg   = 1;              % Minimum glide-slope angle [degrees]
m0              = 2000;           % Initial interceptor mass [kg]
z0              = log(m0);        % Log-mass state (linearizes mass dynamics)

% ==================== Interceptor MPC Discretization =========
N               = 30;             % Number of control intervals
r0              = [400;1000;1000];% Initial interceptor position [km]
v0              = [0;-4;-4];      % Initial interceptor velocity [km/s]
gamma_min_rad   = gamma_min_deg * pi/180;  % Glide-slope constraint [rad]
t_delay         = 20;             % Delay before interceptor launch [s]

% ==================== ICBM Boost-Phase Simulation ===========
t_boost         = 180;            % Total boost-phase duration [s]
N_icbm          = 300;            % Number of time samples for boost simulation
dt_icbm         = t_boost / (N_icbm-1);  % ICBM time step [s]

t_icbm = linspace(0, t_boost, N_icbm)';   % ICBM time vector
r_icbm = zeros(3, N_icbm);                % Position states
v_icbm = zeros(3, N_icbm);                % Velocity states

pitch_max = 45 * pi/180;                  % Max pitch-over angle [rad]
a_mag     = 0.02;                         % Thrust acceleration magnitude [km/s^2]
for i = 2:N_icbm
    tau       = t_icbm(i);
    pitch     = pitch_max * (tau / t_boost);
    a_thrust  = [a_mag * cos(pitch); a_mag * sin(pitch); 0];
    a_total   = a_thrust + g;
    v_icbm(:,i) = v_icbm(:,i-1) + a_total * dt_icbm;
    r_icbm(:,i) = r_icbm(:,i-1) + v_icbm(:,i-1) * dt_icbm;
end
% Reference: Wikipedia “Intercontinental ballistic missile” for realistic boost profile

% ==================== Continuous-Time Dynamics =====================
A_ct = zeros(7);
A_ct(1,4) = 1; A_ct(2,5) = 1; A_ct(3,6) = 1;
B_ct = zeros(7,4);
B_ct(4,1) = 1; B_ct(5,2) = 1; B_ct(6,3) = 1; B_ct(7,4) = -alpha;
c_ct    = zeros(7,1);
c_ct(4:6)= g;

% ==================== Solver Setup & Time-Optimal Search ========
ops     = sdpsettings('solver','mosek','verbose',0);
tol     = 1e-3;        % Bisection tolerance [s]
maxIter = 20;          % Max bisection iterations
T_low   = 0;
T_high  = N * 1.0;     % Initial upper bound guess (s)

% Find feasible upper bound by doubling T_high
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

% Bisection to find T_opt
for iter = 1:maxIter
    T_mid = 0.5 * (T_low + T_high);
    dt_mid = T_mid / N;
    [Ad1,Bd1,cd1] = zohOneStep(A_ct, B_ct, c_ct, dt_mid);
    [~,~,cons1] = createDecisionVarsAndConstraints(Ad1, Bd1, cd1, dt_mid, ...
                     r0, v0, z0, m0, alpha, T_min, T_max, gamma_min_rad, N, [0;0;0]);
    sol = optimize(cons1, [], ops);
    if sol.problem == 0
        T_high = T_mid;
    else
        T_low = T_mid;
    end
    if (T_high - T_low) < tol
        break;
    end
end
T_opt        = T_high;
dt           = T_opt / N;
t_int_global = t_delay + T_opt;

% Compute intercept point on ICBM path
[~, idx_int] = min(abs(t_icbm - t_int_global));
r_int_pt     = interp1(t_icbm, r_icbm', t_int_global)';

% ==================== MPC Loop ===============================
[Ad,Bd,cd] = zohOneStep(A_ct, B_ct, c_ct, dt);
Xcur       = [r0; v0; z0];
Xtraj      = zeros(7, N+1);
Utraj      = zeros(4, N);
Xtraj(:,1) = Xcur;

for k = 1:N
    N_horiz = N - k + 1;

    if k == 15
        Xcur = Xcur + [0; 0; 0; 0; 0; 0; 0];
    end

    [Xh,Uh,cons] = createDecisionVarsAndConstraints(Ad, Bd, cd, dt, ...
                       Xcur(1:3), Xcur(4:6), Xcur(7), m0, alpha, T_min, T_max, gamma_min_rad, N_horiz, r_int_pt);
    sol = optimize(cons, [], ops);
    if sol.problem ~= 0
        error('MPC infeasible at step %d', k);
    end
    Uk = value(Uh{1});
    Utraj(:,k) = Uk;
    % Apply first control and update state

    Xcur = Ad * Xcur + Bd * Uk + cd;
    Xtraj(:,k+1) = Xcur;
end

toc  % End timer

% ==================== Plot MPC Results ==========================
% 3D Trajectory: Interceptor vs. ICBM
figure('Name','MPC Trajectory','Color','white');
plot3(-Xtraj(2,:), Xtraj(3,:), Xtraj(1,:), 'b-o','LineWidth',2); hold on;
plot3(-r_icbm(2,1:idx_int), r_icbm(3,1:idx_int), r_icbm(1,1:idx_int), 'r','LineWidth',2);
plot3(-r_icbm(2,idx_int:end), r_icbm(3,idx_int:end), r_icbm(1,idx_int:end), ...
      'Color',[0.5 0.5 0.5], 'LineWidth',2);
plot3(-r_int_pt(2), r_int_pt(3), r_int_pt(1), 'kp','MarkerSize',25,'MarkerFaceColor','y');
grid on; axis equal;
xlabel('r_2 [km]'); ylabel('r_3 [km]'); zlabel('r_1 [km]');
title('Interceptor vs. ICBM Paths');
legend('Interceptor','ICBM(boost)','ICBM(post)','Intercept Point','Location','best');

% Throttle Profile
t_ctrl = t_delay + (0:N-1)*dt;
throttlePct = (Utraj(4,:) / T_max) * 100;
figure('Name','Throttle Profile','Color','white');
plot(t_ctrl, throttlePct,'m-o','LineWidth',2);
grid on; xlabel('Time [s]'); ylabel('Throttle [%]');
title('MPC Throttle Profile');

% ==================== Local Functions ===========================
function [X, U, cons] = createDecisionVarsAndConstraints(Ad, Bd, cd, dt, ...
    r0, v0, z0, m0, alpha, T_min, T_max, gamma, N, r_target)
    % Build state & input vars and constraints for horizon
    X = cell(N+1,1); U = cell(N,1);
    for i = 1:N+1, X{i} = sdpvar(7,1,'full'); end
    for i = 1:N,   U{i} = sdpvar(4,1,'full');   end
    cons = (X{1} == [r0; v0; z0]);
    % Dynamics
    for k = 1:N
        cons = [cons, X{k+1} == Ad*X{k} + Bd*U{k} + cd];
    end
    % Mass linearization
    z_lb = zeros(N+1,1);
    for i = 0:N
        z_lb(i+1) = log(m0 - alpha * T_max * (i * dt));
    end
    % Throttle & accel constraints
    for k = 1:N
        zk     = X{k}(7);
        ak     = U{k}(1:3);
        sk     = U{k}(4);
        e_lbk  = exp(-z_lb(k));
        lin_app= e_lbk * (1 - (zk - z_lb(k)));
        cons  = [cons,
            sk >= T_min * e_lbk,
            sk <= T_max * lin_app,
            norm(ak,2) <= sk
        ];
    end
    % Terminal intercept if provided
    if ~isempty(r_target)
        cons = [cons, X{end}(1:3) == r_target];
    end
    % Glide-slope
    for k = 1:N+1
        rk = X{k}(1:3);
        cons = [cons, rk(1) >= tan(gamma) * norm(rk(2:3),2)];
    end
end

function [Ad, Bd, cd] = zohOneStep(A, B, c, dt)
    % Zero-order hold discretization
    n = size(A,1);
    M = expm(A * dt);
    if rank(A) == n
        IntA = A \ (M - eye(n));
    else
        steps = 50; dTau = dt/steps; IntA = zeros(n);
        for i = 1:steps
            IntA = IntA + expm(A * (i-1) * dTau) * dTau;
        end
    end
    Ad = M;
    Bd = IntA * B;
    cd = IntA * c;
end