% =========================================================================
% 
% Filename:       HW4.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE590 - Applied Control in Astronautics
% Professor:      Dr. Kenshiro Oguri
% Assignment:     HW 4
% Semester:       Spring 2025
% 
% Description: Homework 4
%
% =========================================================================

%%

clear;clc;close all

mu = 1.0;   % Gravitational parameter (AU^3/year^2)
T_fixed = 8;  % Fixed final time

% Earth initial conditions
r0      = 1.0;
theta0  = 0;
vr0     = 0;
vtheta0 = sqrt(mu / r0);  % = 1 if r0=1 and mu=1

% Mars orbit info
rMars       = 1.524;
nMars       = sqrt(mu / (rMars^3));  % mean motion
vthetaMars  = sqrt(mu / rMars);      % circular orbit speed

% 2) Single-Shooting Unknowns & Initial Guesses
%   X = [lambda_r(0), lambda_theta(0), lambda_vr(0), lambda_vtheta(0)]
lr0_guess   = 0;
lth0_guess  = 0;
lvr0_guess  = 0;
lvth0_guess = 0;

X_guess = [lr0_guess; lth0_guess; lvr0_guess; lvth0_guess];

% 3) Solve the BVP using fsolve
options = optimoptions('fsolve','Display','iter','MaxFunctionEvaluations',1e5, ...
                       'TolFun',1e-8,'TolX',1e-8);

[X_sol, ~, exitflag, output] = fsolve(@(X) boundaryConditions_fixedTime( ...
    X, mu, T_fixed, r0, theta0, vr0, vtheta0, rMars, nMars, vthetaMars), ...
    X_guess, options);

disp('--------------------------------------------------');
disp('RESULTS FROM INDIRECT METHOD (Rendezvous with Mars, FIXED T=8)');
disp('--------------------------------------------------');
if exitflag > 0
    disp('fsolve converged to a solution!');
else
    disp('fsolve did NOT converge properly; try different guesses or tolerances.');
end
disp(output.message);

% Extract solution (just costates, since T is fixed)
lambda_r0   = X_sol(1);
lambda_th0  = X_sol(2);
lambda_vr0  = X_sol(3);
lambda_vth0 = X_sol(4);

disp(['Fixed Transfer Time T = ', num2str(T_fixed), ' years']);
disp('Initial Costates at t=0:');
disp(['  lambda_r(0)      = ', num2str(lambda_r0)]);
disp(['  lambda_theta(0)  = ', num2str(lambda_th0)]);
disp(['  lambda_vr(0)     = ', num2str(lambda_vr0)]);
disp(['  lambda_vtheta(0) = ', num2str(lambda_vth0)]);

% 4) Re-integrate with the Found Solution
X_init = [r0; theta0; vr0; vtheta0; lambda_r0; lambda_th0; lambda_vr0; lambda_vth0];
ode_opts = odeset('MaxStep', 0.01, 'Refine', 4, 'RelTol',1e-10, 'AbsTol',1e-12);
[tSol, X_all] = ode45(@(t,Z) stateCostateOde(t, Z, mu), [0 T_fixed], X_init, ode_opts);

% Extract state
r_sol     = X_all(:,1);
theta_sol = X_all(:,2);
vr_sol    = X_all(:,3);
vth_sol   = X_all(:,4);

% Extract costate
lr_sol    = X_all(:,5);
lth_sol   = X_all(:,6);
lvr_sol   = X_all(:,7);
lvth_sol  = X_all(:,8);

% 5) Convert spacecraft polar to Cartesian for plotting
x_sc = r_sol .* cos(theta_sol);
y_sc = r_sol .* sin(theta_sol);

% Mars's position over time
thetaMars_t = pi + nMars * tSol;  % Mars starts at π
x_mars = rMars * cos(thetaMars_t);
y_mars = rMars * sin(thetaMars_t);

% Earth start position in Cartesian
x_earthStart = r0*cos(0);  % = 1
y_earthStart = r0*sin(0);  % = 0

% Mars start position in Cartesian
x_marsStart = rMars*cos(pi);  % = -rMars
y_marsStart = rMars*sin(pi);  % = 0

% 6) Plot the Transfer Trajectory & Markers
figure('Color','white','Position',[1000 0 1000 1000]);
plot(x_sc, y_sc, 'b-', 'LineWidth', 2, 'DisplayName','Spacecraft'); hold on;
plot(x_mars, y_mars, 'r--', 'LineWidth', 1.5, 'DisplayName','Mars Trajectory');
% Earth orbit for reference:
plot(r0*cos(linspace(0,2*pi,200)), r0*sin(linspace(0,2*pi,200)), 'k--','DisplayName','Earth Orbit');
% Mark Earth start
plot(x_earthStart, y_earthStart, 'bo','MarkerFaceColor','k','MarkerSize',8, ...
    'DisplayName','Earth Start');
% Mark Mars start
plot(x_marsStart, y_marsStart, 'ro','MarkerFaceColor','r','MarkerSize',8, ...
    'DisplayName','Mars Start');
% Mark intercept
scatter(x_sc(end), y_sc(end), 1000, 'kp','MarkerFaceColor','y', ...
    'DisplayName','Intercept');

axis equal; grid on;
legend('Location','Best');
xlabel('x (AU)'); ylabel('y (AU)');
title('Minimum-Energy Transfer for t_f = 8 years');
ax = gca;
ax.FontSize = 20;

% 7) Plot State vs. Time
figure('Color','white','Position',[0 0 1000 700]);
subplot(3,1,1);
plot(tSol, r_sol, 'LineWidth', 2); grid on;
ylabel('r (AU)'); title('States vs. Time');
ax = gca;
ax.FontSize = 20;

subplot(3,1,2);
plot(tSol, vr_sol, 'LineWidth', 2); hold on;
plot(tSol, vth_sol, 'LineWidth', 2);
legend('v_r','v_\theta','Location','Best'); grid on;
ylabel('Velocity (AU/year)');
ax = gca;
ax.FontSize = 20;

subplot(3,1,3);
plot(tSol, theta_sol, 'LineWidth', 2); grid on;
xlabel('Time (years)'); ylabel('\theta (rad)');
ax = gca;
ax.FontSize = 20;

% 8) Compute & Plot H(t) - H(0)
Npts = length(tSol);
H_sol = zeros(Npts,1);
for i = 1:Npts
    % State at i
    rr      = r_sol(i);
    th      = theta_sol(i);
    vrr     = vr_sol(i);
    vthh    = vth_sol(i);
    % Costate at i
    lrr     = lr_sol(i);
    lthh    = lth_sol(i);
    lvrr    = lvr_sol(i);
    lvthh   = lvth_sol(i);
    % Compute Hamiltonian
    H_sol(i) = computeHamiltonian(rr, th, vrr, vthh, ...
                                  lrr, lthh, lvrr, lvthh, mu);
end
dH_sol = H_sol - H_sol(1);

figure('Color','white','Position',[1000 500 1000 700]);
plot(tSol, dH_sol, 'LineWidth',2);
grid on;
xlabel('Time (years)');
ylabel('H(t) - H(0)');
title('Change in Control Hamiltonian Over Time');
ax = gca;
ax.FontSize = 20;

disp('Done.');

% FIGURE 9: Control Time History
Npts = length(tSol);
u1_sol = zeros(Npts,1);
u2_sol = zeros(Npts,1);
for i = 1:Npts
    u1_sol(i) = -lvr_sol(i)/2;
    u2_sol(i) = -lvth_sol(i)/2;
end
figure('Color','white','Position',[0 500 1000 700]); hold on;
plot(tSol, u1_sol, 'b-', 'LineWidth', 2, 'DisplayName', 'u_1');
plot(tSol, u2_sol, 'r-', 'LineWidth', 2, 'DisplayName', 'u_2');
plot(tSol, vecnorm([u1_sol,u2_sol],2,2), 'k-', 'LineWidth', 2, 'DisplayName', 'u_{mag}');
grid on;
xlabel('Time (years)');
ylabel('Control Input');
title('Control Time History');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

% FIGURE 10: Costate Time History
figure('Color','white','Position',[1000 0 1000 700]);
subplot(2,1,1); hold on; grid on;
plot(tSol, lr_sol, 'b-', 'LineWidth', 2, 'DisplayName', '\lambda_r');
plot(tSol, lth_sol, 'r-', 'LineWidth', 2, 'DisplayName', '\lambda_{\theta}');
xlabel('Time (years)');
ylabel('Costate');
title('Costate Time History: \lambda_r and \lambda_{\theta}');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

subplot(2,1,2); hold on; grid on;
plot(tSol, lvr_sol, 'b-', 'LineWidth', 2, 'DisplayName', '\lambda_{vr}');
plot(tSol, lvth_sol, 'r-', 'LineWidth', 2, 'DisplayName', '\lambda_{v\theta}');
xlabel('Time (years)');
ylabel('Costate');
title('Costate Time History: \lambda_{vr} and \lambda_{v\theta}');
legend('Location','Best');
ax = gca; ax.FontSize = 20;

%% 20 initial state guesses
clear; clc; close all;

% 1) Problem Setup (same as your existing code)
mu = 1.0;  % gravitational parameter (AU^3/year^2)

% Earth initial conditions
r0      = 1.0;
theta0  = 0;
vr0     = 0;
vtheta0 = sqrt(mu / r0);  % = 1 if r0=1 and mu=1

% Mars orbit info
rMars       = 1.524;
nMars       = sqrt(mu / (rMars^3));  % mean motion
vthetaMars  = sqrt(mu / rMars);      % circular orbit speed

% We'll fix an initial guess for T (e.g., T=2.0). You can adjust as needed.
fixedT_guess = 8;

% 2) Random Trials
N = 20;               % number of random samples
rng('default');       % for reproducibility (optional)

% Storage for results
X_sol_all      = zeros(N,5);  % will store [T, lr(0), lth(0), lvr(0), lvth(0)]
exitflag_all   = zeros(N,1);
iterations_all = zeros(N,1);

% fsolve options
options = optimoptions('fsolve','Display','none','MaxFunctionEvaluations',1e5,...
                       'TolFun',1e-8,'TolX',1e-8);

for i = 1:N
    % Each costate is uniform(-10,10)
    costateGuess = unifrnd(-10,10,[4,1]);
    X_guess = [fixedT_guess; costateGuess];  % [T, lr(0), lth(0), lvr(0), lvth(0)]
    
    % Solve with fsolve
    [X_sol_i, ~, exitflag_i, output_i] = fsolve(@(X) boundaryConditions(...
        X, mu, r0, theta0, vr0, vtheta0, rMars, nMars, vthetaMars), X_guess, options);
    
    % Store results
    X_sol_all(i,:)      = X_sol_i(:)';
    exitflag_all(i)     = exitflag_i;
    iterations_all(i)   = output_i.iterations;
end

% 3) Summarize Convergence
num_converged = sum(exitflag_all > 0);
disp('------------------------------------------');
disp(['Random Trials: N = ', num2str(N)]);
disp(['Number of solutions that converged = ', num2str(num_converged)]);
if num_converged > 0
    disp('Iteration counts for converged samples:');
    disp(iterations_all(exitflag_all>0)');
else
    disp('No solutions converged with these random guesses.');
end

% 4) If at least one converged, plot up to the last 3 converged solutions
idx_converged = find(exitflag_all>0);
if isempty(idx_converged)
    disp('No converged samples to plot. Done.');
    return;
end

% We'll take up to 3 from the end
numToPlot = min(3, length(idx_converged));
indicesToPlot = idx_converged(end-numToPlot+1 : end);

% Preallocate cell arrays to store data for each solution we plot
all_tSol   = cell(numToPlot,1);
all_rSol   = cell(numToPlot,1);
all_thetaSol = cell(numToPlot,1);
all_vrSol  = cell(numToPlot,1);
all_vthSol = cell(numToPlot,1);
all_xSc    = cell(numToPlot,1);
all_ySc    = cell(numToPlot,1);
all_dH     = cell(numToPlot,1);
labels     = cell(numToPlot,1);

% We'll use some distinct line colors
colorList = {'b','g','m'};  % for up to 3 solutions

for k = 1:numToPlot
    conv_idx = indicesToPlot(k);
    X_sol_final = X_sol_all(conv_idx,:);
    T_opt       = X_sol_final(1);
    lambda_r0   = X_sol_final(2);
    lambda_th0  = X_sol_final(3);
    lambda_vr0  = X_sol_final(4);
    lambda_vth0 = X_sol_final(5);

    disp(' ');
    disp(['--- Converged sample #', num2str(conv_idx), ...
          ' (one of the last ', num2str(numToPlot), ' ) ---']);
    disp(['  T* = ', num2str(T_opt)]);
    disp(['  lambda_r(0)      = ', num2str(lambda_r0)]);
    disp(['  lambda_theta(0)  = ', num2str(lambda_th0)]);
    disp(['  lambda_vr(0)     = ', num2str(lambda_vr0)]);
    disp(['  lambda_vtheta(0) = ', num2str(lambda_vth0)]);

    % Re-integrate
    X_init_final = [r0; theta0; vr0; vtheta0; lambda_r0; lambda_th0; lambda_vr0; lambda_vth0];
    [tSol, X_all] = ode45(@(t,Z) stateCostateOde(t, Z, mu), [0 T_opt], X_init_final);

    r_sol     = X_all(:,1);
    theta_sol = X_all(:,2);
    vr_sol    = X_all(:,3);
    vth_sol   = X_all(:,4);

    % Convert polar to Cartesian
    x_sc = r_sol .* cos(theta_sol);
    y_sc = r_sol .* sin(theta_sol);

    % Compute Hamiltonian difference H(t)-H(0)
    H_sol = zeros(size(tSol));
    for i = 1:length(tSol)
        H_sol(i) = computeHamiltonian(X_all(i,1), X_all(i,2), X_all(i,3), X_all(i,4), ...
                                      X_all(i,5), X_all(i,6), X_all(i,7), X_all(i,8), mu);
    end
    dH = H_sol - H_sol(1);

    % Store in cell arrays
    all_tSol{k}     = tSol;
    all_rSol{k}     = r_sol;
    all_thetaSol{k} = theta_sol;
    all_vrSol{k}    = vr_sol;
    all_vthSol{k}   = vth_sol;
    all_xSc{k}      = x_sc;
    all_ySc{k}      = y_sc;
    all_dH{k}       = dH;
    labels{k}       = ['Sol#', num2str(conv_idx)];
end

% 5) Now we create 3 figures, each with up to 3 solutions

% Figure 1: Trajectories
figure('Color','white','Position',[1000 500 1000 1000]); hold on;
% Plot Earth orbit
plot(r0*cos(linspace(0,2*pi,200)), r0*sin(linspace(0,2*pi,200)), 'k--', 'DisplayName','Earth Orbit');
% Mars orbit param over time if you want it, or just the circle:
plot(rMars*cos(linspace(0,2*pi,200)), rMars*sin(linspace(0,2*pi,200)), 'r--', 'DisplayName','Mars Orbit');
% Mark Earth start
plot(r0, 0, 'ko','MarkerFaceColor','b','MarkerSize',12,'DisplayName','Earth Start');
% Mark Mars start (rMars, angle=pi => (-rMars,0))
plot(-rMars, 0, 'ro','MarkerFaceColor','r','MarkerSize',12,'DisplayName','Mars Start');

startValue = 8;
lineWidthStart = startValue;
decrement = 3;

for k = 1:numToPlot
    plot(all_xSc{k}, all_ySc{k}, [colorList{k},'-'], 'LineWidth',lineWidthStart, ...
        'DisplayName',[labels{k},' Traj']);
    % Mark final intercept
    scatter(all_xSc{k}(end), all_ySc{k}(end), lineWidthStart*600, [colorList{k},'p'], ...
        'filled','MarkerEdgeColor','k','DisplayName',[labels{k},' Intercept']);
    lineWidthStart = lineWidthStart - decrement;
end
axis equal; grid on;
legend('Location','best');
xlabel('x (AU)'); ylabel('y (AU)');
title(['Trajectories of Last ', num2str(numToPlot),' Converged Solutions']);
ax = gca;
ax.FontSize = 20;

% Figure 2: States vs. Time
figure('Color','white','Position',[1000 0 1000 1000])
subplot(3,1,1); hold on; grid on;
title(['States vs. Time for Last ', num2str(numToPlot),' Converged Solutions']);
ylabel('r (AU)');

lineWidthStart = startValue;
for k = 1:numToPlot
    plot(all_tSol{k}, all_rSol{k}, [colorList{k},'-'], 'LineWidth',lineWidthStart, ...
         'DisplayName', labels{k});
    lineWidthStart = lineWidthStart - decrement;
end
legend('Location','best');
ax = gca;
ax.FontSize = 20;

subplot(3,1,2); hold on; grid on;
ylabel('Velocity (AU/yr)');

lineWidthStart = startValue;
for k = 1:numToPlot
    plot(all_tSol{k}, all_vrSol{k}, [colorList{k},'--'], 'LineWidth',lineWidthStart, ...
         'DisplayName',[labels{k},' v_r']);
    plot(all_tSol{k}, all_vthSol{k}, [colorList{k},'-'], 'LineWidth',lineWidthStart, ...
         'DisplayName',[labels{k},' v_\theta']);
    lineWidthStart = lineWidthStart - decrement;
end
legend('Location','best');
ax = gca;
ax.FontSize = 20;

subplot(3,1,3); hold on; grid on;
xlabel('Time (yr)'); ylabel('\theta (rad)');

lineWidthStart = startValue;
for k = 1:numToPlot
    plot(all_tSol{k}, all_thetaSol{k}, [colorList{k},'-'], 'LineWidth',lineWidthStart, ...
         'DisplayName',[labels{k},' \theta']);
    lineWidthStart = lineWidthStart - decrement;
end
legend('Location','best');
ax = gca;
ax.FontSize = 20;

% Figure 3: Hamiltonian difference
figure('Color','white','Position',[0 0 1000 1000]); hold on; grid on;
title(['H(t)-H(0) for Last ', num2str(numToPlot),' Converged Solutions']);
xlabel('Time (yr)'); ylabel('H(t)-H(0)');

lineWidthStart = startValue;
for k = 1:numToPlot
    plot(all_tSol{k}, all_dH{k}, [colorList{k},'-'], 'LineWidth',lineWidthStart, ...
         'DisplayName', labels{k});
    lineWidthStart = lineWidthStart - decrement;
end
legend('Location','best');
ax = gca;
ax.FontSize = 20;

disp('Done.');

%%

%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%

function res = boundaryConditions_fixedTime(X, mu, Tfixed, ...
    r0, theta0, vr0, vtheta0, rMars, nMars, vthetaMars)
    % X = [lambda_r(0), lambda_theta(0), lambda_vr(0), lambda_vtheta(0)]
    lr0  = X(1);
    lth0 = X(2);
    lvr0 = X(3);
    lvth0= X(4);

    % Construct initial condition for state+costate
    Z_init = [r0; theta0; vr0; vtheta0; lr0; lth0; lvr0; lvth0];

    % Integrate from t=0 to t=Tfixed
    ode_opts_inner = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.01);
    [~, Z_out] = ode45(@(t,Z) stateCostateOde(t, Z, mu), ...
                       [0 Tfixed], Z_init, ode_opts_inner);

    % Final state
    rEnd     = Z_out(end,1);
    thetaEnd = Z_out(end,2);
    vrEnd    = Z_out(end,3);
    vthEnd   = Z_out(end,4);

    % Boundary conditions for fixed final time Tfixed:
    % (1) r(T) = rMars
    % (2) theta(T) = π + nMars*T
    % (3) v_r(T) = 0
    % (4) v_theta(T) = vthetaMars
    res = zeros(4,1);
    res(1) = rEnd - rMars;
    res(2) = thetaEnd - (pi + nMars*Tfixed);
    res(3) = vrEnd;
    res(4) = vthEnd - vthetaMars;
end

function dZdt = stateCostateOde(~, Z, mu)
    % Unpack the state variables and costates:
    % State: [r, theta, vr, vtheta]
    % Costate: [lambda_r, lambda_theta, lambda_vr, lambda_vtheta]
    r      = Z(1);
    theta  = Z(2);
    vr     = Z(3);
    vtheta = Z(4);
    
    lr     = Z(5);
    lth    = Z(6);
    lvr    = Z(7);
    lvth   = Z(8);
    
    % Optimal control (from ∂H/∂u = 0)
    u1 = -lvr / 2;
    u2 = -lvth / 2;
    
    % State dynamics
    drdt = vr;
    dthetadt = vtheta / r;
    dvrdt = (vtheta^2)/r - mu/(r^2) + u1;
    dvthetadt = -(vr*vtheta)/r + u2;
    
    % Costate dynamics: dot(lambda) = -∂H/∂x.
    dHdr = - lth*(vtheta/(r^2)) ...
           + lvr*( - (vtheta^2)/(r^2) + 2*mu/(r^3)) ...
           + lvth*( (vr*vtheta)/(r^2));
    dlr_dt = - dHdr;
    
    dHdtheta = 0;
    dlth_dt = - dHdtheta;
    
    dHdvr = lr - lvth*(vtheta/r);
    dlvr_dt = - dHdvr;
    
    dHdvtheta = lth*(1/r) + lvr*(2*vtheta/r) - lvth*(vr/r);
    dlvth_dt = - dHdvtheta;
    
    dZdt = [drdt; dthetadt; dvrdt; dvthetadt; dlr_dt; dlth_dt; dlvr_dt; dlvth_dt];
end

function Hval = computeHamiltonian(r,th,vr,vth, lr,ltheta,lvr,lvth, mu)
    % Optimal control for min-energy => u1 = -lvr/2, u2 = -lvth/2
    u1 = -lvr/2;
    u2 = -lvth/2;

    % L = u1^2 + u2^2
    L = u1^2 + u2^2;

    % State derivatives f(x,u):
    fr   = vr;
    fth  = vth / r;
    fvr  = (vth^2)/r - mu/(r^2) + u1;
    fvth = -(vr*vth)/r + u2;

    lamTf = lr*fr + ltheta*fth + lvr*fvr + lvth*fvth;

    % Hamiltonian = L + λ^T f
    Hval = L + lamTf;
end



