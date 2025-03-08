%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: earthMarsIndirect_fixedT.m
% Purpose: Solve the Earth-to-Mars minimum-energy transfer problem with
%          a FIXED final time (t_f = 8 years). We do multiple random trials
%          for the costates, check convergence, and then plot the converged
%          solutions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

% 1) Problem Setup
mu = 1.0;  % gravitational parameter (AU^3/year^2)

% Earth initial conditions
r0      = 1.0;    % radius (AU)
theta0  = 0;      % angular position (rad)
vr0     = 0;      % radial velocity (AU/yr)
vtheta0 = sqrt(mu / r0);  % tangential velocity (AU/yr) => 1 if r0=1 and mu=1

% Mars orbital info
rMars       = 1.524;                  % radius of Mars orbit (AU)
nMars       = sqrt(mu / (rMars^3));   % mean motion (rad/yr)
vthetaMars  = sqrt(mu / rMars);       % tangential velocity for Mars (AU/yr)

% FIXED final time
Tfixed = 8.0;   % years

% 2) Random Trials
N = 10;               % number of random samples
rng(1);       % for reproducibility (optional)

% We'll store the costate solutions, exit flags, iteration counts
costate_sol_all = zeros(N,4);  % [lr(0), lth(0), lvr(0), lvth(0)]
exitflag_all    = zeros(N,1);
iterations_all  = zeros(N,1);

% fsolve options -- the ONLY change here is 'Display','iter-detailed'
options = optimoptions('fsolve', ...
    'Display','iter-detailed', ...  % <--- changed from 'none' to show progress
    'MaxFunctionEvaluations',1e5, ...
    'TolFun',1e-8, ...
    'TolX',1e-8);

for i = 1:N
    % Each costate is uniform(-10,10)
    costateGuess = unifrnd(-10,10,[4,1]);
    
    % Solve with fsolve
    [X_sol_i, ~, exitflag_i, output_i] = fsolve(@(X) ...
        boundaryConditions_fixedTime(X, mu, Tfixed, ...
                                     r0, theta0, vr0, vtheta0, ...
                                     rMars, nMars, vthetaMars), ...
                                     costateGuess, options);
    
    % Store results
    costate_sol_all(i,:) = X_sol_i(:).';
    exitflag_all(i)      = exitflag_i;
    iterations_all(i)    = output_i.iterations;
end

% 3) Summarize Convergence
num_converged = sum(exitflag_all > 0);
disp('------------------------------------------');
disp(['Random Trials: N = ', num2str(N)]);
disp(['Number of solutions that converged = ', num2str(num_converged)]);
if num_converged > 0
    disp('Iteration counts for converged samples:');
    disp(iterations_all(exitflag_all>0).');
else
    disp('No solutions converged with these random guesses.');
end

% Indices of converged solutions
idx_converged = find(exitflag_all>0);
if isempty(idx_converged)
    disp('No converged samples to plot. Done.');
    return;
end

% 4) Re-integrate and Plot up to the last 3 converged solutions
numToPlot = min(3, length(idx_converged));
indicesToPlot = idx_converged(end-numToPlot+1 : end);

% Preallocate cell arrays to store data for each solution we plot
all_tSol      = cell(numToPlot,1);
all_rSol      = cell(numToPlot,1);
all_thetaSol  = cell(numToPlot,1);
all_vrSol     = cell(numToPlot,1);
all_vthSol    = cell(numToPlot,1);
all_xSc       = cell(numToPlot,1);
all_ySc       = cell(numToPlot,1);
all_dH        = cell(numToPlot,1);
labels        = cell(numToPlot,1);

% We'll use some distinct line colors
colorList = {'b','g','m'};  % for up to 3 solutions

for k = 1:numToPlot
    conv_idx = indicesToPlot(k);
    lambda_r0   = costate_sol_all(conv_idx,1);
    lambda_th0  = costate_sol_all(conv_idx,2);
    lambda_vr0  = costate_sol_all(conv_idx,3);
    lambda_vth0 = costate_sol_all(conv_idx,4);
    
    disp(' ');
    disp(['--- Converged sample #', num2str(conv_idx), ...
          ' (one of the last ', num2str(numToPlot), ' ) ---']);
    disp(['  lambda_r(0)      = ', num2str(lambda_r0)]);
    disp(['  lambda_theta(0)  = ', num2str(lambda_th0)]);
    disp(['  lambda_vr(0)     = ', num2str(lambda_vr0)]);
    disp(['  lambda_vtheta(0) = ', num2str(lambda_vth0)]);
    
    % Re-integrate state+costate from t=0 to t=Tfixed
    Z_init_final = [r0; theta0; vr0; vtheta0; ...
                    lambda_r0; lambda_th0; lambda_vr0; lambda_vth0];

    ode_opts = odeset('MaxStep',0.01,'Refine',4,'RelTol',1e-10,'AbsTol',1e-12);
    [tSol, Z_out] = ode45(@(t,Z) stateCostateOde(t, Z, mu), ...
                          [0 Tfixed], Z_init_final,ode_opts);
    
    r_sol     = Z_out(:,1);
    theta_sol = Z_out(:,2);
    vr_sol    = Z_out(:,3);
    vth_sol   = Z_out(:,4);
    
    % Convert polar to Cartesian
    x_sc = r_sol .* cos(theta_sol);
    y_sc = r_sol .* sin(theta_sol);

    % Compute Hamiltonian difference H(t) - H(0)
    H_sol = zeros(size(tSol));
    for j = 1:length(tSol)
        H_sol(j) = computeHamiltonian(Z_out(j,1), Z_out(j,2), ...
                                      Z_out(j,3), Z_out(j,4), ...
                                      Z_out(j,5), Z_out(j,6), ...
                                      Z_out(j,7), Z_out(j,8), mu);
    end
    dH = H_sol - H_sol(1);

    % Store in cell arrays
    all_tSol{k}      = tSol;
    all_rSol{k}      = r_sol;
    all_thetaSol{k}  = theta_sol;
    all_vrSol{k}     = vr_sol;
    all_vthSol{k}    = vth_sol;
    all_xSc{k}       = x_sc;
    all_ySc{k}       = y_sc;
    all_dH{k}        = dH;
    labels{k}        = ['Sol#', num2str(conv_idx)];
end

% ------------------ PLOTTING ------------------ %
% We define a decreasing linewidth for clarity
startValue = 8;
decrement  = 3;

% Figure 1: Trajectories
figure('Color','white','Position',[1000 500 1000 1000]); hold on;
% Plot Earth orbit
plot(r0*cos(linspace(0,2*pi,200)), r0*sin(linspace(0,2*pi,200)), ...
     'k--', 'DisplayName','Earth Orbit');
% Plot Mars orbit
plot(rMars*cos(linspace(0,2*pi,200)), rMars*sin(linspace(0,2*pi,200)), ...
     'r--', 'DisplayName','Mars Orbit');
% Mark Earth start
plot(r0, 0, 'ko','MarkerFaceColor','b','MarkerSize',12,'DisplayName','Earth Start');
% Mark Mars start (theta = pi => x = -rMars, y=0)
plot(-rMars, 0, 'ro','MarkerFaceColor','r','MarkerSize',12,'DisplayName','Mars Start');

lineWidthStart = startValue;
for k = 1:numToPlot
    plot(all_xSc{k}, all_ySc{k}, [colorList{k},'-'], ...
         'LineWidth', lineWidthStart, 'DisplayName', [labels{k},' Traj']);
    % Mark final intercept
    scatter(all_xSc{k}(end), all_ySc{k}(end), lineWidthStart*600, ...
            [colorList{k},'p'], 'filled','MarkerEdgeColor','k', ...
            'DisplayName',[labels{k},' Intercept']);
    lineWidthStart = lineWidthStart - decrement;
end
axis equal; grid on; legend('Location','best');
xlabel('x (AU)'); ylabel('y (AU)');
title(['Trajectories of Last ', num2str(numToPlot),' Converged Solutions']);
set(gca, 'FontSize', 20);

% Figure 2: States vs. Time
figure('Color','white','Position',[1000 0 1000 1000]);

% r(t)
subplot(3,1,1); hold on; grid on;
title(['States vs. Time for Last ', num2str(numToPlot),' Converged Solutions']);
ylabel('r (AU)');
lineWidthStart = startValue;
for k = 1:numToPlot
    plot(all_tSol{k}, all_rSol{k}, [colorList{k},'-'], ...
         'LineWidth', lineWidthStart, 'DisplayName', labels{k});
    lineWidthStart = lineWidthStart - decrement;
end
legend('Location','best');
set(gca, 'FontSize', 20);

% v_r(t) and v_theta(t)
subplot(3,1,2); hold on; grid on;
ylabel('Velocity (AU/yr)');
lineWidthStart = startValue;
for k = 1:numToPlot
    plot(all_tSol{k}, all_vrSol{k},  [colorList{k},'--'], ...
         'LineWidth', lineWidthStart, 'DisplayName',[labels{k},' v_r']);
    plot(all_tSol{k}, all_vthSol{k}, [colorList{k},'-'], ...
         'LineWidth', lineWidthStart, 'DisplayName',[labels{k},' v_\theta']);
    lineWidthStart = lineWidthStart - decrement;
end
legend('Location','best');
set(gca, 'FontSize', 20);

% theta(t)
subplot(3,1,3); hold on; grid on;
xlabel('Time (yr)'); ylabel('\theta (rad)');
lineWidthStart = startValue;
for k = 1:numToPlot
    plot(all_tSol{k}, all_thetaSol{k}, [colorList{k},'-'], ...
         'LineWidth', lineWidthStart, 'DisplayName',[labels{k},' \theta']);
    lineWidthStart = lineWidthStart - decrement;
end
legend('Location','best');
set(gca, 'FontSize', 20);

% Figure 3: Hamiltonian difference
figure('Color','white','Position',[0 0 1000 1000]); hold on; grid on;
title(['H(t) - H(0) for Last ', num2str(numToPlot),' Converged Solutions']);
xlabel('Time (yr)'); ylabel('H(t) - H(0)');
lineWidthStart = startValue;
for k = 1:numToPlot
    plot(all_tSol{k}, all_dH{k}, [colorList{k},'-'], ...
         'LineWidth', lineWidthStart, 'DisplayName', labels{k});
    lineWidthStart = lineWidthStart - decrement;
end
legend('Location','best');
set(gca, 'FontSize', 20);

disp('Done.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = boundaryConditions_fixedTime(X, mu, Tfixed, ...
                                            r0, theta0, vr0, vtheta0, ...
                                            rMars, nMars, vthetaMars)
% ------------------------------------------------------------------------
% X = [ lambda_r(0), lambda_theta(0), lambda_vr(0), lambda_vtheta(0) ]
%
% For a FIXED final time Tfixed, the boundary conditions are:
%   1) r(T) = rMars
%   2) theta(T) = pi + nMars*Tfixed
%   3) v_r(T) = 0
%   4) v_theta(T) = vthetaMars
% ------------------------------------------------------------------------

    lr0  = X(1);
    lth0 = X(2);
    lvr0 = X(3);
    lvth0= X(4);

    % Construct initial condition for state+costate
    Z_init = [r0; theta0; vr0; vtheta0; lr0; lth0; lvr0; lvth0];

    % Integrate from t=0 to t=Tfixed
    ode_opts_inner = odeset('RelTol',1e-10, 'AbsTol',1e-12, 'MaxStep',0.01);
    [~, Z_out] = ode45(@(t,Z) stateCostateOde(t, Z, mu), [0 Tfixed], Z_init, ode_opts_inner);

    % Final state
    rEnd     = Z_out(end,1);
    thetaEnd = Z_out(end,2);
    vrEnd    = Z_out(end,3);
    vthEnd   = Z_out(end,4);

    % The 4 boundary conditions:
    % (1) r(T) = rMars
    % (2) theta(T) = pi + nMars*Tfixed
    % (3) v_r(T) = 0
    % (4) v_theta(T) = vthetaMars
    res = zeros(4,1);
    res(1) = rEnd - rMars;
    res(2) = thetaEnd - (pi + nMars*Tfixed);
    res(3) = vrEnd;
    res(4) = vthEnd - vthetaMars;
end

function dZdt = stateCostateOde(~, Z, mu)
% ------------------------------------------------------------------------
% ODE for the state+costate system under the cost function J = \int (u^2 dt)
%
%   State:   [r, theta, v_r, v_theta]
%   Costate: [\lambda_r, \lambda_\theta, \lambda_{v_r}, \lambda_{v_\theta}]
%
%   Control: u = (u1, u2) = -0.5 * [\lambda_{v_r}, \lambda_{v_\theta}]
% ------------------------------------------------------------------------
    % Unpack states
    r      = Z(1);
    theta  = Z(2);
    vr     = Z(3);
    vth    = Z(4);

    % Unpack costates
    lr     = Z(5);
    lth    = Z(6);
    lvr    = Z(7);
    lvth   = Z(8);

    % Optimal control (from ∂H/∂u = 0 => u = -1/2 * lambda_v)
    u1 = -0.5 * lvr;
    u2 = -0.5 * lvth;

    % State dynamics
    drdt     = vr;
    dthetadt = vth / r;
    dvrdt    = (vth^2)/r - mu/(r^2) + u1;
    dvthdt   = -(vr*vth)/r + u2;

    % Costate dynamics: dot(\lambda) = -∂H/∂x
    %
    % Hamiltonian H = u^2 + λ_r * f_r + λ_θ * f_θ + λ_{v_r} * f_{v_r} + λ_{v_θ} * f_{v_θ}
    % with u^2 = (u1^2 + u2^2).
    %
    % We'll compute partial derivatives carefully:

    % dH/dr
    dHdr = ...  
         - lth * (vth / (r^2)) ...
         + lvr * ( - (vth^2)/(r^2) + 2*mu/(r^3) ) ...
         + lvth * ( (vr*vth)/(r^2) );
    dlr_dt = -dHdr;

    % dH/dθ = 0 in this problem (no explicit dependence except costates),
    % but let's keep it for clarity:
    dHdtheta = 0;
    dlth_dt  = -dHdtheta;

    % dH/d(v_r)
    %   f_r = v_r => ∂f_r/∂v_r = 1
    %   f_{v_θ} = -(v_r v_θ)/r + u2 => partial wrt v_r = -(v_θ)/r
    dHdvr = lr - lvth*(vth/r);
    dlvr_dt = -dHdvr;

    % dH/d(v_θ)
    %   f_θ = v_θ / r => ∂f_θ/∂v_θ = 1/r
    %   f_{v_r} = (v_θ^2)/r - ... => partial wrt v_θ = 2*v_θ / r
    %   f_{v_θ} = -(v_r v_θ)/r + u2 => partial wrt v_θ = -v_r / r
    dHdvtheta = lth*(1/r) + lvr*(2*vth / r) - lvth*(vr / r);
    dlvth_dt  = -dHdvtheta;

    dZdt = [drdt; dthetadt; dvrdt; dvthdt; dlr_dt; dlth_dt; dlvr_dt; dlvth_dt];
end

function Hval = computeHamiltonian(r,th,vr,vth, lr,ltheta,lvr,lvth, mu)
% ------------------------------------------------------------------------
% Compute the Hamiltonian:
%    H = u^2 + λ^T f(x,u)
%
% where u = (u1, u2) = -0.5*[lvr, lvth].
% ------------------------------------------------------------------------
    % Optimal control
    u1 = -0.5 * lvr;
    u2 = -0.5 * lvth;

    % Cost function integrand L = u1^2 + u2^2
    L = u1^2 + u2^2;

    % State dynamics
    fr   = vr;
    fth  = vth / r;
    fvr  = (vth^2)/r - mu/(r^2) + u1;
    fvth = -(vr*vth)/r + u2;

    % λ^T f
    lamTf = lr*fr + ltheta*fth + lvr*fvr + lvth*fvth;

    % Hamiltonian
    Hval = L + lamTf;
end