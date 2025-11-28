% =========================================================================
% 
% Filename:       HW6.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE590 - Applied Control in Astronautics
% Professor:      Dr. Kenshiro Oguri
% Assignment:     HW 6
% Semester:       Spring 2025
% 
% Description: Homework 6
%
% =========================================================================

%% Problem 6

clear;clc;close all

% Define the continuous-time matrices A, B, c as in your figure:

alpha = 0.5086; % s/km
g = [-3.7114E-3, 0, 0]';
r0 = [1.5, 0, 2.0]'; % Initial position [km]
v0 = [-0.075, 0.03, 0.1]'; % Initial velocity [km/s]

m0 = 2000; % Initial mass [kg]
z0 = log(m0);

A = [
  0 0 0 1 0 0 0
  0 0 0 0 1 0 0
  0 0 0 0 0 1 0
  0 0 0 0 0 0 0
  0 0 0 0 0 0 0
  0 0 0 0 0 0 0
  0 0 0 0 0 0 0
];
B = [
  0  0  0  0
  0  0  0  0
  0  0  0  0
  1  0  0  0
  0  1  0  0
  0  0  1  0
  0  0  0 -alpha  
];
c = [
  0
  0
  0
  g
  0
];

% Build the function handle for the ODE:
fEoM = @(t,x,u) A*x + B*u + c;

% Choose a time grid tspan for discretization (N steps):
N = 1;            % number of intervals
t0 = 0; tF = 10;   % start at t=0, end at t=5
tspan = linspace(t0, tF, N+1);  % e.g. [0,1,2,3,4,5]

% Initial state x0 and a piecewise-constant control sequence
x0 = [r0; v0; z0];  % 7 states: (r1,r2,r3,v1,v2,v3,z)
% Suppose we keep the same control (u_k) = (sigma, a1, a2, a3) = (0,0,0,0)
% for demonstration. In practice, you would define actual thrust inputs.
uSeq = cell(1, N);
for k = 1:N
    uSeq{k} = [0; 0; 0; 0];  % constant control in each interval
end

options = odeset('AbsTol',1E-12,'RelTol',1E-12);

% Call the function to get discrete A_k, B_k, c_k for each interval
[Acell, Bcell, ccell, xnext] = computeDiscreteMatricesZOH(fEoM, tspan, x0, uSeq,options);


% Display results for each step:
for k = 1:N
    fprintf('--- Step k=%d ---\n', k);
    fprintf('Acell{k}:\n'), disp(Acell{k});
    fprintf('Bcell{k}:\n'), disp(Bcell{k});
    fprintf('ccell{k}:\n'), disp(ccell{k});
    % fprintf('xnext{k}:\n'), disp(xnext{k});
end

%% 

clear; clc; close all

% Problem Parameters
alpha = 0.5086;           % fuel consumption coeff [s/km]
g = [-3.7114e-3; 0; 0];   % gravity in km/s^2 (acting in -r1 direction)
T_min = 4.97;             % kg·km/s^2
T_max = 13.26;            % kg·km/s^2
gamma_min_deg = 3;        % glide slope angle in degrees

m0  = 2000;               % initial mass [kg]
z0  = log(m0);            % log of initial mass

% Discretization
dt = 1.0;                 
N  = 78;                  
t0 = 0;                   
tF = N*dt;                

% Initial state (r1=vertical)
r0 = [1.5; 0; 2.0];       % [km]
v0 = [-0.075; 0.03; 0.1]; % [km/s]

% Continuous-time system: x = [r1; r2; r3; v1; v2; v3; z]
A_ct = zeros(7);
A_ct(1,4) = 1;  
A_ct(2,5) = 1;  
A_ct(3,6) = 1;

B_ct = zeros(7,4);
B_ct(4,1) = 1;       % a1 -> v1
B_ct(5,2) = 1;       % a2 -> v2
B_ct(6,3) = 1;       % a3 -> v3
B_ct(7,4) = -alpha;  % sigma -> z
c_ct = zeros(7,1);
c_ct(4:6) = g;

yalmip('clear');

% YALMIP Setup
X = cell(N+1,1);
U = cell(N,1);
constraints = [];

% Initial condition
X{1} = sdpvar(7,1,'full');
constraints = [constraints, X{1} == [r0; v0; z0]];

for k = 2:(N+1)
   X{k} = sdpvar(7,1,'full'); 
end
for k = 1:N
   U{k} = sdpvar(4,1,'full');  % [a1; a2; a3; sigma]
end

% Discrete-Time Dynamics
[Ad, Bd, cd] = zohOneStep(A_ct, B_ct, c_ct, dt);
for k = 1:N
    xk = X{k};
    uk = U{k};
    xk1 = X{k+1};
    constraints = [constraints, xk1 == Ad*xk + Bd*uk + cd];
end

% Thrust / Norm Constraints
z_lb = zeros(N+1,1);
for i = 0:N
   mass_lb = m0 - alpha*T_max*(i*dt);
   z_lb(i+1) = log(mass_lb);
end

for k = 1:N
   xk = X{k};
   uk = U{k};
   z_k = xk(7);

   % Now a_k is the first three entries, sigma_k is the 4th
   a_k = uk(1:3);
   sigma_k = uk(4);

   e_lbk      = exp(-z_lb(k));
   lin_approx = e_lbk*(1 - (z_k - z_lb(k))); 

   constraints = [constraints, sigma_k >= T_min * e_lbk];
   constraints = [constraints, sigma_k <= T_max * lin_approx];
   constraints = [constraints, norm(a_k,2) <= sigma_k];
end

% Boundary Conditions
constraints = [constraints, X{N+1}(1:3)==0, X{N+1}(4:6)==0];

% Glide-Slope Constraint
gamma_min_rad = gamma_min_deg*(pi/180);
for k = 1:N+1
    rk = X{k}(1:3);
    constraints = [constraints, rk(1) >= tan(gamma_min_rad)*norm(rk(2:3),2)];
end

% Objective: sum of sigma_k * dt
sum_sigma = 0;
for k = 1:N
    sum_sigma = sum_sigma + U{k}(4);
end
Objective = sum_sigma * dt;

% Solve
ops = sdpsettings('solver','mosek','verbose',1);
sol = optimize(constraints, Objective, ops);

if sol.problem == 0
    disp('======================================================');
    disp('    Feasible and optimal solution found! ');
    disp('======================================================');
    disp(['  Minimum sum of sigma_k*dt = ', num2str(value(Objective),'%.4f')]);

    xFinal = value(X{N+1});
    finalMass = exp(xFinal(7));
    consumedFuel = m0 - finalMass;
    disp(['  Final mass: ', num2str(finalMass,'%.4f'), ' kg']);
    disp(['  Fuel consumed: ', num2str(consumedFuel,'%.4f'), ' kg']);

    % Plot difference: sigma_k - ||a_k||
    plotLosslessProperty(U, N, dt);

    % Plot the rest
    plotResults(X, U, N, dt, t0, tF);
else
    warning(['Solver reported problem code = ', num2str(sol.problem)]);
    disp(sol.info);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZOH discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ad, Bd, cd] = zohOneStep(A,B,c,dt)
    n = size(A,1);
    M = expm(A*dt);
    if rank(A) == n
       Ainv = inv(A);
       IntA = Ainv*(M - eye(n));
    else
       steps = 50;
       dTau  = dt/steps;
       IntA = zeros(n);
       for i = 1:steps
           tau = (i-1)*dTau;
           IntA = IntA + expm(A*tau)*dTau;
       end
    end
    Ad = M;
    Bd = IntA*B;
    cd = IntA*c;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the "lossless" difference: sigma_k - ||a_k||_2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotLosslessProperty(U, N, dt)
    diffVals = zeros(1,N);
    for k = 1:N
        % Now a_k is U{k}(1:3), sigma_k is U{k}(4)
        a_val    = value(U{k}(1:3));
        sigma_val= value(U{k}(4));
        diffVals(k) = sigma_val - norm(a_val,2);
    end
    t_ctrl = linspace(0, dt*(N-1), N);

    figure('Color','white','Position',[0,100,1000,600]);
    plot(t_ctrl, diffVals, 'k-o','LineWidth',2);
    grid on; hold on;
    xlabel('Time [s]', 'FontSize',20);
    ylabel('\sigma_k - ||a_k||_2', 'FontSize',20);
    title('Lossless Property: \sigma_k - ||a_k||_2 vs. Time', 'FontSize',20);
    set(gca,'FontSize',20);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotResults(X, U, N, dt, t0, tF)
    Xval = zeros(7, N+1);
    for k = 1:N+1
        Xval(:,k) = value(X{k});
    end
    
    rVals = Xval(1:3, :);  
    vVals = Xval(4:6, :);  
    zVals = Xval(7, :);
    massVals = exp(zVals);

    for k = 1:N
        aVals(:,k) = (vVals(:,k+1) - vVals(:,k))/dt;
    end
    
    for k = 1:N
        sigmaVals(k) = value(U{k}(4));
    end
    
    t_state = linspace(t0, tF, N+1);
    t_acc   = t_state(1:end-1) + dt/2;
    t_ctrl  = linspace(t0, tF-dt, N);
    
    % 3D Position Plot
    figure('Color','white','Position',[500,500,1000,800]);
    plot3(-rVals(2,:), rVals(3,:), rVals(1,:), 'b-o', 'LineWidth',2);
    grid on; hold on;
    xlabel('r2 [km]', 'FontSize',20);
    ylabel('r3 [km]', 'FontSize',20);
    zlabel('r1 [km] (vertical)', 'FontSize',20);
    title('3D Position Trajectory', 'FontSize',20);
    set(gca,'FontSize',20);
    view(75,20)

    % 2D Position Plot
    figure('Color','white','Position',[500,500,1000,500]);
    plot3(-rVals(2,:), rVals(3,:), rVals(1,:), 'b-o', 'LineWidth',2);
    grid on; hold on;
    xlabel('r2 [km]', 'FontSize',20);
    ylabel('r3 [km]', 'FontSize',20);
    zlabel('r1 [km] (vertical)', 'FontSize',20);
    set(gca,'FontSize',20);
    view(90,0)
    
    figure('Color','white','Position',[500,0,1000,500]);
    plot(t_state, rVals(1,:), 'r-o', t_state, rVals(2,:), 'b-o', t_state, rVals(3,:), 'g-o','LineWidth',2);
    grid on; hold on;
    xlabel('Time [s]', 'FontSize',20);
    ylabel('Position [km]', 'FontSize',20);
    legend('r1 (vertical)','r2','r3','FontSize',20);
    title('Position Components vs. Time', 'FontSize',20);
    set(gca,'FontSize',20);
    
    figure('Color','white','Position',[3000,1000,1000,500]);
    plot(t_state, vVals(1,:), 'r-o', t_state, vVals(2,:), 'b-o', t_state, vVals(3,:), 'g-o','LineWidth',2);
    grid on; hold on;
    xlabel('Time [s]', 'FontSize',20);
    ylabel('Velocity [km/s]', 'FontSize',20);
    legend('v1 (vertical)','v2','v3','FontSize',20);
    title('Velocity Components vs. Time', 'FontSize',20);
    set(gca,'FontSize',20);
    
    figure('Color','white','Position',[3000,500,1000,500]);
    plot(t_acc, aVals(1,:), 'r-o', t_acc, aVals(2,:), 'b-o', t_acc, aVals(3,:), 'g-o','LineWidth',2);
    grid on; hold on;
    xlabel('Time [s]', 'FontSize',20);
    ylabel('Acceleration [km/s^2]', 'FontSize',20);
    legend('a1 (vertical)','a2','a3','FontSize',20);
    title('Acceleration Components vs. Time', 'FontSize',20);
    set(gca,'FontSize',20);
    
    throttleLevel = (sigmaVals / (13.26e-3)) * 100;  % T_max=0.01326 in kg·km/s^2
    figure('Color','white','Position',[3000,0,1000,500]);
    plot(t_ctrl, throttleLevel, 'm-o','LineWidth',2);
    grid on; hold on;
    xlabel('Time [s]','FontSize',20);
    ylabel('Throttle Level [%]','FontSize',20);
    title('Throttle Level vs. Time','FontSize',20);
    set(gca,'FontSize',20);
end

function [Acell, Bcell, ccell, xnext] = computeDiscreteMatricesZOH(fEoM, ...
                                     tspan, x0, uSeq, options)
%COMPUTEDISCRETEMATRICESZOH 
%  Computes (A_k, B_k, c_k) for each discrete step under ZOH
%  using numeric integration of an augmented system.
%
% Inputs:
%   fEoM   = handle to continuous dynamics: fEoM(t, x, u) -> xdot
%   tspan  = [t0, t1, ..., tN], time nodes
%   x0     = initial state at t0
%   uSeq   = cell array or matrix of control inputs [u0, u1, ..., u_{N-1}]
%   options= ODE solver options (optional)
%
% Outputs:
%   Acell{k} = A_k
%   Bcell{k} = B_k
%   ccell{k} = c_k
%   xnext{k} = x_{k+1} from the actual integration
%
% Each step integrates from t_k to t_{k+1} with u(t)=u_k (ZOH).

N = length(tspan) - 1;  % number of intervals
nx = length(x0);        % dimension of the state
nu = length(uSeq{1});   % dimension of the control

Acell = cell(1, N);
Bcell = cell(1, N);
ccell = cell(1, N);
xnext = cell(1, N);

xk = x0;  % current state
for k = 1:N
    
    tk     = tspan(k);
    tkp1   = tspan(k+1);
    uk     = uSeq{k};   % constant over [tk, tk+1]
    dt     = tkp1 - tk;
    
    % Build augmented initial condition:
    % Y = [x; reshape(Ix, nx^2,1); reshape(Iu, nx*nu,1)] 
    % where Ix=eye(nx), Iu= zeros(...) if you want partial w.r.t u, etc.
    
    Ix = eye(nx);
    Iu = zeros(nx, nu); 
    % or if your system is affine in u, you might keep partial derivatives in another structure
    Y0 = [xk; reshape(Ix, nx^2, 1); reshape(Iu, nx*nu, 1)];
    
    % Integrate from tk to tk+1:
    sol = ode45(@(t, Y) augEoM(t, Y, uk, fEoM, nx, nu), [tk tkp1], Y0, options);
    
    Yfinal = deval(sol, tkp1);  % Evaluate at end
    xkp1   = Yfinal(1:nx);      % This is x_{k+1}
    
    % Parse out partial derivatives wrt x:
    offset1 = nx; 
    PHI_x   = reshape(Yfinal(offset1+1 : offset1+nx^2), nx, nx);
    
    % partial derivatives wrt u (if you included them):
    offset2 = offset1 + nx^2;
    PHI_u   = reshape(Yfinal(offset2+1 : offset2+nx*nu), nx, nu);
    
    % Store results:
    xnext{k} = xkp1;
    
    % Now define A_k, B_k, c_k:
    % Suppose the continuous system is: dot{x} = A(t)*x + B(t)*u + c(t) (affine).
    % Then from the augmented integration, PHI_x ~ A_k, PHI_u ~ B_k, 
    % but there's also a constant shift (the integral of c(t)).
    
    Acell{k} = PHI_x;    % partial wrt x
    Bcell{k} = PHI_u;    % partial wrt u
    ccell{k} = xkp1 - PHI_x*xk - PHI_u*uk;
    
    % Update xk for next step
    xk = xkp1;
end

end

%% ---------------------------------------------------------------------
function dYdt = augEoM(t, Y, uk, fEoM, nx, nu)
% augEoM: builds an augmented system to track partial derivatives wrt x,u.
%
% Y = [x; vec(PHI_x); vec(PHI_u)]
% We'll differentiate each piece w.r.t time.

% Extract x from Y:
x    = Y(1:nx);
ix1  = nx;
ix2  = nx + nx^2;
ix3  = nx + nx^2 + nx*nu; 

% The continuous dynamics:
xdot = fEoM(t, x, uk);  % must be a function handle

% If fEoM is affine: xdot = A(t)*x + B(t)*u + c(t), 
% then partial_x = A(t), partial_u = B(t).
% Otherwise, you'd linearize or compute Jacobians numerically here.

% For demonstration, let's just do a generic approach:
A_t  = approximateJacobianX(@(xx) fEoM(t, xx, uk), x);   % Nx-by-Nx
B_t  = approximateJacobianU(@(uu) fEoM(t, x, uu), uk);   % Nx-by-Nu
% c(t) is included in xdot if the system is not purely linear

% Extract the PHI_x, PHI_u from Y:
PHI_x = reshape(Y(ix1+1:ix2), nx, nx);
PHI_u = reshape(Y(ix2+1:ix3), nx, nu);

% Now define derivatives:
dxdT    = xdot;
dPHI_xdT= A_t * PHI_x;  % chain rule
dPHI_udT= A_t * PHI_u + B_t;  % chain rule

% Pack up:
dYdt = [dxdT;
        reshape(dPHI_xdT, nx*nx, 1);
        reshape(dPHI_udT, nx*nu, 1)];
end

%% -------------- (example numeric Jacobian approximations) --------------
function A = approximateJacobianX(fun, x)
% Simple finite-difference for partial w.r.t x
nx = length(x);
A  = zeros(nx, nx);
epsVal = 1e-6;
f0 = fun(x);
for i = 1:nx
   dx      = zeros(nx,1);
   dx(i)   = epsVal;
   fPlus   = fun(x+dx);
   A(:, i) = (fPlus - f0)/epsVal;
end
end

function B = approximateJacobianU(fun, u)
% Simple finite-difference for partial w.r.t u
nu = length(u);
f0 = fun(u);
nx = length(f0);
B  = zeros(nx, nu);
epsVal = 1e-6;
for j = 1:nu
   du      = zeros(nu,1);
   du(j)   = epsVal;
   fPlus   = fun(u+du);
   B(:, j) = (fPlus - f0)/epsVal;
end
end