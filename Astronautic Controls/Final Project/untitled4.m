% optimized_interceptor_real_time.m
% Time‐Optimal Interceptor with Real‐Time ICBM Feedback
clear; clc; close all; tic;
set(groot,'defaultAxesFontSize',20,'defaultTextFontSize',20);

% ==================== Problem Parameters (Earth Model) ====================
alpha           = 0.0;            % Fuel consumption coefficient [s/km]
test_alpha      = 0.3;            % (unused here, but kept)
g               = [-9.81e-3; 0; 0];% Gravity vector [km/s^2]
T_min           = 1;              % Minimum throttle scalar [kg·km/s^2]
T_max           = 500;            % Maximum throttle scalar [kg·km/s^2]
gamma_min_deg   = 1;              % Glide-slope angle [deg]
gamma_min_rad   = gamma_min_deg * pi/180;

m0              = 2000;           % Initial interceptor mass [kg]
z0              = log(m0);        % Log-mass state

% ==================== Discretization & Initial States ====================
N               = 30;             % # control intervals
r0              = [400;1000;1000];% Initial interceptor pos [km]
v0              = [0;-4;-4];      % Initial interceptor vel [km/s]
t_delay         = 20;             % Launch delay [s]

% ==================== ICBM Boost-Phase Simulation ====================
t_boost         = 180;            % Boost duration [s]
N_icbm          = 300;            % # samples
dt_icbm         = t_boost/(N_icbm-1);

t_icbm = linspace(0, t_boost, N_icbm)';
r_icbm = zeros(3, N_icbm);
v_icbm = zeros(3, N_icbm);
pitch_max = 45*pi/180;
a_mag     = 0.02;

for i = 2:N_icbm
    pitch     = pitch_max * (t_icbm(i)/t_boost);
    a_thrust  = [a_mag*cos(pitch); a_mag*sin(pitch); 0];
    a_total   = a_thrust + g;
    v_icbm(:,i) = v_icbm(:,i-1) + a_total*dt_icbm;
    r_icbm(:,i) = r_icbm(:,i-1) + v_icbm(:,i-1)*dt_icbm;
end

% ==================== Continuous-Time Dynamics ====================
A_ct = zeros(7);
A_ct(1,4)=1; A_ct(2,5)=1; A_ct(3,6)=1;
B_ct = zeros(7,4);
B_ct(4:6,1:3) = eye(3);
B_ct(7,4) = -alpha;
c_ct = zeros(7,1);
c_ct(4:6) = g;

% ==================== MPC with Real-Time ICBM Feedback ==============
Xcur  = [r0; v0; z0];
Xtraj = zeros(7, N+1);  Xtraj(:,1)=Xcur;
Utraj = zeros(4, N);
ops    = sdpsettings('solver','mosek','verbose',0);

for k = 1:N
    % 1) Current absolute time
    t_now = t_delay + (k-1)*(t_boost/(N-1));
    % 2) Interpolate ICBM position at t_now
    r_target = interp1(t_icbm, r_icbm', t_now, 'linear', 'extrap')';
    % 3) Discretize dynamics for this control step
    dt = t_boost/(N-1);
    [Ad,Bd,cd] = zohOneStep(A_ct, B_ct, c_ct, dt);

    % 4) Build & solve MPC over remaining horizon
    N_horiz = N - k + 1;
    [Xh, Uh, cons] = createDecisionVarsAndConstraints(...
        Ad, Bd, cd, dt, ...
        Xcur(1:3), Xcur(4:6), Xcur(7), ...
        m0, alpha, T_min, T_max, gamma_min_rad, ...
        N_horiz, r_target);
    sol = optimize(cons, [], ops);
    if sol.problem
        error('MPC infeasible at step %d', k);
    end
    Uk = value(Uh{1});
    Utraj(:,k) = Uk;

    % 5) Apply control & propagate  
    Xcur = Ad*Xcur + Bd*Uk + cd;
    Xtraj(:,k+1) = Xcur;
end

toc  % total runtime

% ==================== Plot Results ==========================
figure('Name','Interceptor vs ICBM','Color','white');
plot3(-Xtraj(2,:), Xtraj(3,:), Xtraj(1,:), 'b-o','LineWidth',2); hold on;
plot3(-r_icbm(2,:), r_icbm(3,:), r_icbm(1,:), 'r-','LineWidth',2);
grid on; axis equal;
xlabel('r_2 [km]','FontSize',20);
ylabel('r_3 [km]','FontSize',20);
zlabel('r_1 [km]','FontSize',20);
title('Interceptor vs. Real-Time ICBM Trajectories','FontSize',20);
legend('Interceptor','ICBM','Location','best','FontSize',16);

% ==================== Local Functions ============================
function [X, U, cons] = createDecisionVarsAndConstraints(...
    Ad, Bd, cd, dt, r0, v0, z0, m0, alpha, T_min, T_max, gamma, N, r_target)

    X = cell(N+1,1); U = cell(N,1);
    for i=1:N+1, X{i}=sdpvar(7,1,'full'); end
    for i=1:N,   U{i}=sdpvar(4,1,'full');   end

    cons = (X{1} == [r0; v0; z0]);
    for i=1:N
        cons = [cons, X{i+1} == Ad*X{i} + Bd*U{i} + cd];
    end

    z_lb = arrayfun(@(i) log(m0 - alpha*T_max*(i*dt)), 0:N)';
    for i=1:N
        zk = X{i}(7); ak = U{i}(1:3); sk = U{i}(4);
        e_lbk = exp(-z_lb(i));
        lin_app = e_lbk*(1 - (zk - z_lb(i)));
        cons = [cons,
            sk >= T_min*e_lbk,
            sk <= T_max*lin_app,
            norm(ak,2) <= sk
        ];
    end

    % Terminal position = real-time ICBM pos
    cons = [cons, X{end}(1:3) == r_target];

    for i=1:N+1
        rk = X{i}(1:3);
        cons = [cons,
            rk(1) >= tan(gamma)*norm(rk(2:3),2)
        ];
    end
end

function [Ad, Bd, cd] = zohOneStep(A, B, c, dt)
    n = size(A,1);
    M = expm(A*dt);
    if rank(A)==n
        IntA = A\(M-eye(n));
    else
        steps=50; dTau=dt/steps; IntA=zeros(n);
        for j=1:steps
            IntA = IntA + expm(A*(j-1)*dTau)*dTau;
        end
    end
    Ad = M;
    Bd = IntA * B;
    cd = IntA * c;
end