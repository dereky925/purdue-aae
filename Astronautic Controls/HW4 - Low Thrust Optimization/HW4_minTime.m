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
clear; clc; close all;

% 1) Problem Setup
mu = 1.0;  % nondimensional gravitational parameter

% Earth (departure) initial conditions
r0 = 1.0;
theta0 = 0;
vr0 = 0;
vtheta0 = sqrt(mu/r0);  % = 1 if r0 = 1, mu = 1

% Mars orbit conditions
rMars = 1.524;
vthetaMars = sqrt(mu/rMars);      % Circular orbit speed at Mars
nMars = sqrt(mu/(rMars^3));         % Mars' mean motion

% 2) Shooting Unknowns & Initial Guesses
% Unknowns: X = [ T, lambda_r(0), lambda_theta(0), lambda_vr(0), lambda_vtheta(0)]
T_guess     = 6.5;   
lr0_guess   = 5.0;
lth0_guess  = 2.0;
lvr0_guess  = 2.0;
lvth0_guess = 7.0;

X_guess = [T_guess; lr0_guess; lth0_guess; lvr0_guess; lvth0_guess];

% 3) Solve the BVP using fsolve
options = optimoptions('fsolve','Display','iter','MaxFunctionEvaluations',1e5,...
                       'TolFun',1e-8,'TolX',1e-8);

% Note: boundaryConditions now takes nMars and vthetaMars as additional arguments.
[X_sol, ~, exitflag, output] = fsolve(@(X) boundaryConditions(X, mu, ...
    r0, theta0, vr0, vtheta0, rMars, nMars, vthetaMars), X_guess, options);

disp('--------------------------------------------------');
disp(' RESULTS FROM INDIRECT METHOD (Time-Minimization with Clamped Control)');
disp('--------------------------------------------------');
if exitflag > 0
    disp('fsolve converged to a solution!');
else
    disp('fsolve did NOT converge properly; try different guesses or tolerances.');
end
disp(output.message);

% Extract solution
T_opt       = X_sol(1);
lambda_r0   = X_sol(2);
lambda_th0  = X_sol(3);
lambda_vr0  = X_sol(4);
lambda_vth0 = X_sol(5);

disp(['Optimal Transfer Time T* = ', num2str(T_opt), ' years']);
disp('Costates at t=0:');
disp(['  lambda_r(0)      = ', num2str(lambda_r0)]);
disp(['  lambda_theta(0)  = ', num2str(lambda_th0)]);
disp(['  lambda_vr(0)     = ', num2str(lambda_vr0)]);
disp(['  lambda_vtheta(0) = ', num2str(lambda_vth0)]);

% 4) Re-integrate the ODE with the Found Solution
X_init = [r0; theta0; vr0; vtheta0; lambda_r0; lambda_th0; lambda_vr0; lambda_vth0];
ode_opts = odeset('RelTol',1e-12, 'AbsTol',1e-12);
[tSol, X_all] = ode45(@(t,Z) stateCostateOde(t, Z, mu), [0 T_opt], X_init,ode_opts);

r_sol     = X_all(:,1);
theta_sol = X_all(:,2);
vr_sol    = X_all(:,3);
vtheta_sol= X_all(:,4);

% (Costates are in columns 5-8, not used in plotting here)

% 5) Convert to Cartesian & Determine Key Positions
% Spacecraft trajectory in Cartesian
x_sc = r_sol .* cos(theta_sol);
y_sc = r_sol .* sin(theta_sol);

% Mars trajectory (Mars starts at π)
theta_Mars = pi + nMars * tSol;
x_mars = rMars * cos(theta_Mars);
y_mars = rMars * sin(theta_Mars);

% Earth start position (at t = 0): (1, 0)
x_earthStart = r0 * cos(0);
y_earthStart = r0 * sin(0);

% Mars start position (at t = 0): (rMars*cos(pi), rMars*sin(pi))
x_marsStart = rMars * cos(pi);
y_marsStart = rMars * sin(pi);

% Intercept point (final spacecraft state)
x_intercept = x_sc(end);
y_intercept = y_sc(end);

% 6) Plot Trajectory & Markers
figure('Color','white','Position',[0 0 1000 1000]);
plot(x_sc, y_sc, 'b-', 'LineWidth',2, 'DisplayName','Spacecraft Trajectory'); hold on;
plot(x_mars, y_mars, 'r--', 'LineWidth',1.5, 'DisplayName','Mars Trajectory');
% Earth orbit for reference (circle with radius = 1)
plot(r0*cos(linspace(0,2*pi,200)), r0*sin(linspace(0,2*pi,200)), 'k--','DisplayName','Earth Orbit');
% Mark Earth start
plot(x_earthStart, y_earthStart, 'ko', 'MarkerFaceColor','b','MarkerSize',8, 'DisplayName','Earth Start');
% Mark Mars start
plot(x_marsStart, y_marsStart, 'ro', 'MarkerFaceColor','r','MarkerSize',8, 'DisplayName','Mars Start');
% Mark intercept (final state) with a yellow star with black edges
plot(x_intercept, y_intercept, 'p', 'MarkerFaceColor','yellow','MarkerEdgeColor','k','MarkerSize',50, 'DisplayName','Intercept');
axis equal; grid on;
legend('Location','Best');
xlabel('x (AU)'); ylabel('y (AU)');
title('Minimum-Time Transfer with Mars Motion, T^* = 6.3374 years');
ax = gca;
ax.FontSize = 20;

% 7) Plot States vs. Time
figure('Color','white','Position',[1000 0 1000 700]);
subplot(3,1,1);
plot(tSol, r_sol, 'LineWidth',2); grid on;
ylabel('r (AU)'); title('States vs. Time');
ax = gca;
ax.FontSize = 20;

subplot(3,1,2);
plot(tSol, vr_sol, 'LineWidth',2); hold on;
plot(tSol, vtheta_sol, 'LineWidth',2);
legend('v_r','v_\theta','Location','Best'); grid on;
ylabel('Velocity (AU/yr)');
ax = gca;
ax.FontSize = 20;

subplot(3,1,3);
plot(tSol, theta_sol, 'LineWidth',2); grid on;
xlabel('Time (yr)'); ylabel('\theta (rad)');
ax = gca;
ax.FontSize = 20;

% 8) Compute and Plot Hamiltonian Difference H(t)-H(0)
H_sol = zeros(size(tSol));
for i = 1:length(tSol)
    H_sol(i) = computeHamiltonian(r_sol(i), theta_sol(i), vr_sol(i), vtheta_sol(i), ...
                                   X_all(i,5), X_all(i,6), X_all(i,7), X_all(i,8), mu);
end
dH = H_sol - H_sol(1);
figure('Color','white','Position',[1000 500 1000 700]);
plot(tSol, dH, 'LineWidth',2); grid on;
xlabel('Time (yr)'); ylabel('H(t) - H(0)');
title('Time History of [H(t) - H(0)]');
ax = gca;
ax.FontSize = 20;

% % Figure 4: Control Input Time History
% % Compute the control inputs from the costates (using u1 = -lvr/2, u2 = -lvth/2)
% u1 = -X_all(:,7) / 2;
% u2 = -X_all(:,8) / 2;
% figure('Color','white','Position',[1300 500 1000 700]); hold on; grid on;
% plot(tSol, u1, 'b-', 'LineWidth',2, 'DisplayName','u_1');
% plot(tSol, u2, 'r-', 'LineWidth',2, 'DisplayName','u_2');
% legend('Location','best');
% xlabel('Time (yr)');
% ylabel('Control Input');
% title('Control Input Time History');
% ax = gca; ax.FontSize = 20;

% Figure 5: Costate Time History
% Extract costate trajectories from the integration (columns 5 to 8)
lambda_r_sol = X_all(:,5);
lambda_th_sol = X_all(:,6);
lambda_vr_sol = X_all(:,7);
lambda_vth_sol = X_all(:,8);

figure('Color','white','Position',[1300 0 1000 700]);
subplot(2,1,1); hold on; grid on;
plot(tSol, lambda_r_sol, 'b-', 'LineWidth',2, 'DisplayName','\lambda_r');
plot(tSol, lambda_th_sol, 'r-', 'LineWidth',2, 'DisplayName','\lambda_{\theta}');
ylabel('Costate (Position)');
title('Costate Time History');
legend('Location','best');
ax = gca; ax.FontSize = 20;

subplot(2,1,2); hold on; grid on;
plot(tSol, lambda_vr_sol, 'b-', 'LineWidth',2, 'DisplayName','\lambda_{vr}');
plot(tSol, lambda_vth_sol, 'r-', 'LineWidth',2, 'DisplayName','\lambda_{v\theta}');
xlabel('Time (yr)');
ylabel('Costate (Velocity)');
legend('Location','best');
ax = gca; ax.FontSize = 20;

% Compute control input over time using bangBangControl (u* = [u1, u2])
u1 = zeros(length(tSol), 1);
u2 = zeros(length(tSol), 1);
for i = 1:length(tSol)
    % Columns 7 and 8 of X_all are the costates associated with vr and vtheta.
    [u1(i), u2(i)] = bangBangControl(X_all(i,7), X_all(i,8));
end

figure('Color','white','Position',[1000 500 1000 700]); hold on; grid on;
plot(tSol, u1, 'b-', 'LineWidth', 2, 'DisplayName', 'u_1^*');
plot(tSol, u2, 'r-', 'LineWidth', 2, 'DisplayName', 'u_2^*');
plot(tSol, vecnorm([u1,u2],2,2), 'k-', 'LineWidth', 2, 'DisplayName', 'u_{mag}^*');
legend('Location','best');
xlabel('Time (yr)');
ylabel('Control Input (u*)');
title('Control Input Time History');
ax = gca;
ax.FontSize = 20;

%%
clear; clc; close all;
mu = 1.0;  % gravitational parameter (AU^3/year^2)

% Earth initial conditions
r0 = 1.0;
theta0 = 0;
vr0 = 0;
vtheta0 = sqrt(mu/r0);  % = 1

% Mars orbit info
rMars = 1.524;
nMars = sqrt(mu/(rMars^3));  % mean motion
vthetaMars = sqrt(mu/rMars); % circular orbit speed at Mars

% 2) Random Trials for Minimum-Time Problem
% Unknowns: X = [T, lambda_r(0), lambda_theta(0), lambda_vr(0), lambda_vtheta(0)]
% fixedT_guess = 6.5;  % typical guess for transfer time
N = 20;  % number of random trials
rng('default');  % for reproducibility

% Preallocate storage
X_sol_all = zeros(N,5);  % each row: [T, lr(0), lth(0), lvr(0), lvth(0)]
exitflag_all = zeros(N,1);
iterations_all = zeros(N,1);

fs_opts = optimoptions('fsolve','Display','none','MaxFunctionEvaluations',1e5,...
    'TolFun',1e-8,'TolX',1e-8);

for i = 1:N

    fixedT_guess = unifrnd(5,10);
    costateGuess = unifrnd(-10,10,[4,1]);  % random from U(-10,10)

    X_guess = [fixedT_guess; costateGuess];  % [T; lambda_r(0); lambda_theta(0); lambda_vr(0); lambda_vtheta(0)]
    [X_sol_i, ~, exitflag_i, output_i] = fsolve(@(X) boundaryConditions(X, mu, ...
        r0, theta0, vr0, vtheta0, rMars, nMars, vthetaMars), X_guess, fs_opts);
    X_sol_all(i,:) = X_sol_i(:)';
    exitflag_all(i) = exitflag_i;
    iterations_all(i) = output_i.iterations;
end

num_converged = sum(exitflag_all > 0);
disp('------------------------------------------');
disp(['Random Trials: N = ', num2str(N)]);
disp(['Number of converged solutions: ', num2str(num_converged)]);
if num_converged > 0
    disp('Iteration counts for converged samples:');
    disp(iterations_all(exitflag_all>0)');
else
    disp('No solutions converged with these random guesses.');
end

% 3) Extract Last 3 Converged Solutions and Store Data
idx_conv = find(exitflag_all > 0);
if isempty(idx_conv)
    disp('No converged solutions to plot.'); return;
end

numToPlot = min(3, length(idx_conv));
plotIndices = idx_conv(end-numToPlot+1 : end);

% Preallocate cell arrays for plotting data
all_tSol = cell(numToPlot,1);
all_rSol = cell(numToPlot,1);
all_thetaSol = cell(numToPlot,1);
all_vrSol = cell(numToPlot,1);
all_vthSol = cell(numToPlot,1);
all_xSc = cell(numToPlot,1);
all_ySc = cell(numToPlot,1);
all_dH = cell(numToPlot,1);
labels = cell(numToPlot,1);

% Also preallocate costate cell arrays (columns 5 to 8)
all_lambda_r = cell(numToPlot,1);
all_lambda_theta = cell(numToPlot,1);
all_lambda_vr = cell(numToPlot,1);
all_lambda_vth = cell(numToPlot,1);

% Define colors and line widths for clarity
colorList = {'b','g','m'};  % blue, green, magenta
lineWidthList = [8, 5, 2];

for k = 1:numToPlot
    idx = plotIndices(k);
    X_sol_final = X_sol_all(idx,:);
    T_opt = X_sol_final(1);
    lambda_r0 = X_sol_final(2);
    lambda_th0 = X_sol_final(3);
    lambda_vr0 = X_sol_final(4);
    lambda_vth0 = X_sol_final(5);
    
    labels{k} = ['Sol#', num2str(idx)];
    
    % Re-integrate the ODE for this solution
    X_init = [r0; theta0; vr0; vtheta0; lambda_r0; lambda_th0; lambda_vr0; lambda_vth0];
    ode_opts = odeset('MaxStep',0.01,'Refine',4,'RelTol',1e-10,'AbsTol',1e-12);
    [tSol, X_all_i] = ode45(@(t,Z) stateCostateOde(t,Z,mu), [0 T_opt], X_init, ode_opts);
    
    r_sol = X_all_i(:,1);
    theta_sol = X_all_i(:,2);
    vr_sol = X_all_i(:,3);
    vth_sol = X_all_i(:,4);
    lambda_r_sol = X_all_i(:,5);
    lambda_th_sol = X_all_i(:,6);
    lambda_vr_sol = X_all_i(:,7);
    lambda_vth_sol = X_all_i(:,8);
    
    all_tSol{k} = tSol;
    all_rSol{k} = r_sol;
    all_thetaSol{k} = theta_sol;
    all_vrSol{k} = vr_sol;
    all_vthSol{k} = vth_sol;
    all_lambda_r{k} = lambda_r_sol;
    all_lambda_theta{k} = lambda_th_sol;
    all_lambda_vr{k} = lambda_vr_sol;
    all_lambda_vth{k} = lambda_vth_sol;
    
    % Convert polar to Cartesian for plotting trajectory
    x_sc = r_sol .* cos(theta_sol);
    y_sc = r_sol .* sin(theta_sol);
    all_xSc{k} = x_sc;
    all_ySc{k} = y_sc;
    
    % Compute Hamiltonian difference H(t)-H(0)
    H_sol = zeros(length(tSol),1);
    for i = 1:length(tSol)
        H_sol(i) = computeHamiltonian(X_all_i(i,1), X_all_i(i,2), X_all_i(i,3), X_all_i(i,4), ...
                                       X_all_i(i,5), X_all_i(i,6), X_all_i(i,7), X_all_i(i,8), mu);
    end
    dH = H_sol - H_sol(1);
    all_dH{k} = dH;
end

% (Existing Figures for Trajectories, States, and Hamiltonian)
figure('Color','white','Position',[0 0 1000 1000]); hold on;
plot(r0*cos(linspace(0,2*pi,200)), r0*sin(linspace(0,2*pi,200)), 'k--', 'DisplayName','Earth Orbit');
plot(rMars*cos(linspace(0,2*pi,200)), rMars*sin(linspace(0,2*pi,200)), 'r--', 'DisplayName','Mars Orbit');
plot(r0, 0, 'ko','MarkerFaceColor','b','MarkerSize',20, 'DisplayName','Earth Start');
plot(rMars*cos(pi), rMars*sin(pi), 'ro','MarkerFaceColor','r','MarkerSize',20, 'DisplayName','Mars Start');
for k = 1:numToPlot
    lw = lineWidthList(k);
    col = colorList{k};
    plot(all_xSc{k}, all_ySc{k}, [col,'-'], 'LineWidth', lw, 'DisplayName', [labels{k},' Traj']);
    scatter(all_xSc{k}(end), all_ySc{k}(end), lw*600, [col,'p'], 'MarkerFaceColor', col, 'MarkerEdgeColor','k', 'DisplayName', [labels{k},' Intercept']);
end
axis equal; grid on;
legend('Location','Best');
xlabel('x (AU)'); ylabel('y (AU)');
title(['Trajectories of Last ', num2str(numToPlot), ' Converged Solutions']);
ax = gca; ax.FontSize = 20;

figure('Color','white','Position',[1000 0 1000 1000])
subplot(3,1,1); hold on; grid on;
title(['States vs. Time for Last ', num2str(numToPlot), ' Converged Solutions']);
ylabel('r (AU)');
lineWidthStart = lineWidthList(1);
for k = 1:numToPlot
    plot(all_tSol{k}, all_rSol{k}, colorList{k}, 'LineWidth', lineWidthList(k), 'DisplayName', labels{k});
end
legend('Location','Best'); ax = gca; ax.FontSize = 20;
subplot(3,1,2); hold on; grid on;
ylabel('Velocity (AU/yr)');
for k = 1:numToPlot
    plot(all_tSol{k}, all_vrSol{k}, [colorList{k},'--'], 'LineWidth', lineWidthList(k), 'DisplayName', [labels{k},' v_r']);
    plot(all_tSol{k}, all_vthSol{k}, colorList{k}, 'LineWidth', lineWidthList(k), 'DisplayName', [labels{k},' v_\theta']);
end
legend('Location','Best'); ax = gca; ax.FontSize = 20;
subplot(3,1,3); hold on; grid on;
xlabel('Time (yr)'); ylabel('\theta (rad)');
for k = 1:numToPlot
    plot(all_tSol{k}, all_thetaSol{k}, colorList{k}, 'LineWidth', lineWidthList(k), 'DisplayName', [labels{k},' \theta']);
end
legend('Location','Best'); ax = gca; ax.FontSize = 20;

figure('Color','white','Position',[0 1000 1000 1000]); hold on; grid on;
title(['H(t)-H(0) for Last ', num2str(numToPlot), ' Converged Solutions']);
xlabel('Time (yr)'); ylabel('H(t)-H(0)');
for k = 1:numToPlot
    plot(all_tSol{k}, all_dH{k}, colorList{k}, 'LineWidth', lineWidthList(k), 'DisplayName', labels{k});
end
legend('Location','Best'); ax = gca; ax.FontSize = 20;

% Costate Time History for Last 3 Solutions
figure('Color','white','Position',[1300 0 1000 700]);
subplot(2,1,1); hold on; grid on;
title(['Costate Time History for Last ', num2str(numToPlot), ' Solutions']);
for k = 1:numToPlot
    plot(all_tSol{k}, all_lambda_r{k}, colorList{k}, 'LineWidth', lineWidthList(k), 'DisplayName', [labels{k},' \lambda_r']);
    plot(all_tSol{k}, all_lambda_theta{k}, [colorList{k},'--'], 'LineWidth', lineWidthList(k), 'DisplayName', [labels{k},' \lambda_{\theta}']);
end
ylabel('Costate (Position)');
legend('Location','Best'); ax = gca; ax.FontSize = 20;
subplot(2,1,2); hold on; grid on;
for k = 1:numToPlot
    plot(all_tSol{k}, all_lambda_vr{k}, colorList{k}, 'LineWidth', lineWidthList(k), 'DisplayName', [labels{k},' \lambda_{vr}']);
    plot(all_tSol{k}, all_lambda_vth{k}, [colorList{k},'--'], 'LineWidth', lineWidthList(k), 'DisplayName', [labels{k},' \lambda_{v\theta}']);
end
xlabel('Time (yr)');
ylabel('Costate (Velocity)');
legend('Location','Best'); ax = gca; ax.FontSize = 20;

% Control Input Time History for Last 3 Solutions
figure('Color','white','Position',[1300 500 1000 700]); hold on; grid on;
for k = 1:numToPlot
    t_current = all_tSol{k};
    lambda_vr = all_lambda_vr{k};
    lambda_vth = all_lambda_vth{k};
    u1 = zeros(length(t_current),1);
    u2 = zeros(length(t_current),1);
    for i = 1:length(t_current)
        [u1(i), u2(i)] = bangBangControl(lambda_vr(i), lambda_vth(i));
    end
    u_mag = vecnorm([u1,u2],2,2);
    plot(t_current, u1, colorList{k}, 'LineWidth', lineWidthList(k), 'DisplayName', [labels{k},' u_1^*']);
    plot(t_current, u2, [colorList{k},'--'], 'LineWidth', lineWidthList(k), 'DisplayName', [labels{k},' u_2^*']);
    plot(t_current, u_mag, [colorList{k},':'], 'LineWidth', lineWidthList(k), 'DisplayName', [labels{k},' |u|^*']);
end
legend('Location','Best');
xlabel('Time (yr)');
ylabel('Control Input (u*)');
title('Control Input Time History for Last 3 Solutions');
ax = gca; ax.FontSize = 20;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = boundaryConditions(X, mu, r0, th0, vr0, vtheta0, rMars, nMars, vthetaMars)
    % X = [T, lambda_r(0), lambda_theta(0), lambda_vr(0), lambda_vtheta(0)]
    T = X(1);
    lr0 = X(2);
    lth0 = X(3);
    lvr0 = X(4);
    lvth0 = X(5);
    
    if T <= 0
        res = 1e6 * ones(5,1);
        return;
    end
    
    % Initial conditions for spacecraft and costates
    Xinit = [r0; th0; vr0; vtheta0; lr0; lth0; lvr0; lvth0];
    [~, Z_out] = ode45(@(t,Z) stateCostateOde(t,Z,mu), [0 T], Xinit);
    
    % Final states and costates
    rEnd     = Z_out(end,1);
    thetaEnd = Z_out(end,2);
    vrEnd    = Z_out(end,3);
    vthEnd   = Z_out(end,4);
    
    lrEnd    = Z_out(end,5);
    lthEnd   = Z_out(end,6);
    lvrEnd   = Z_out(end,7);
    lvthEnd  = Z_out(end,8);
    
    % Compute the optimal control at final time with clamping
    [u1End, u2End] = bangBangControl(lvrEnd, lvthEnd);
    
    % Compute Hamiltonian at final time
    Hfinal = 1 ...
             + lrEnd * vrEnd ...
             + lthEnd * (vthEnd / rEnd) ...
             + lvrEnd * ((vthEnd^2)/rEnd - mu/(rEnd^2) + u1End) ...
             + lvthEnd * (-(vrEnd*vthEnd)/rEnd + u2End);
    
    % Final boundary conditions:
    % (1) r(T) = rMars,
    % (2) theta(T) = π + nMars*T,
    % (3) v_r(T) = 0,
    % (4) v_theta(T) = vthetaMars,
    % (5) H(T) = 0.
    res = zeros(5,1);
    res(1) = rEnd - rMars;
    res(2) = thetaEnd - (pi + nMars * T);
    res(3) = vrEnd;
    res(4) = vthEnd - vthetaMars;
    res(5) = Hfinal;
end

function dZdt = stateCostateOde(~, Z, mu)
    % Z = [r, theta, vr, vtheta, lambda_r, lambda_theta, lambda_vr, lambda_vtheta]
    r = Z(1);
    theta = Z(2);
    vr = Z(3);
    vtheta = Z(4);
    
    lr = Z(5);
    lth = Z(6);
    lvr = Z(7);
    lvth = Z(8);
    
    % Compute clamped control (if unconstrained u exceeds 0.1, clamp to 0.1)
    [u1, u2] = bangBangControl(lvr, lvth);
    
    % State dynamics
    drdt = vr;
    dthetadt = vtheta / r;
    dvrdt = (vtheta^2)/r - mu/(r^2) + u1;
    dvthetadt = -(vr*vtheta)/r + u2;
    
    % Costate dynamics: dot(lambda) = -∂H/∂x.
    dHdr = - lth*(vtheta/(r^2)) + lvr*( - (vtheta^2)/(r^2) + 2*mu/(r^3)) + lvth*( (vr*vtheta)/(r^2));
    dlr_dt = - dHdr;
    
    dHdtheta = 0;  % no explicit dependence on theta
    dlth_dt = - dHdtheta;
    
    dHdvr = lr - lvth*(vtheta / r);
    dlvr_dt = - dHdvr;
    
    dHdvtheta = lth*(1/r) + lvr*(2*vtheta/r) - lvth*(vr/r);
    dlvth_dt = - dHdvtheta;
    
    dZdt = [drdt; dthetadt; dvrdt; dvthetadt; dlr_dt; dlth_dt; dlvr_dt; dlvth_dt];
end

function [u1, u2] = bangBangControl(lvr, lvth)
    % Compute the unconstrained control as: u_unc = -0.5*[lvr, lvth]
    u_unc = -0.5 * [lvr, lvth];
    u_max = 0.1;
    norm_u = norm(u_unc);

    if norm_u > u_max
        % Clamp the control: scale it to have magnitude u_max
        u1 = u_max * u_unc(1) / norm_u;
        u2 = u_max * u_unc(2) / norm_u;
    else
        % Use the unconstrained control
        u1 = u_unc(1);
        u2 = u_unc(2);
    end
end

function Hval = computeHamiltonian(r, th, vr, vth, lr, lth, lvr, lvth, mu)
    % Compute the Hamiltonian for the time-minimization problem
    [u1, u2] = bangBangControl(lvr, lvth);
    Hval = 1 ...
           + lr * vr ...
           + lth * (vth / r) ...
           + lvr * ((vth^2)/r - mu/(r^2) + u1) ...
           + lvth * (-(vr*vth)/r + u2);
end