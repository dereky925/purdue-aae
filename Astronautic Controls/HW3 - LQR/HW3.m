% =========================================================================
% 
% Filename:       HW3.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE590 - Applied Control in Astronautics
% Professor:      Dr. Kenshiro Oguri
% Assignment:     HW 2
% Semester:       Spring 2025
% 
% Description: Homework 3
%
% =========================================================================

%% Part c Notation from Class
tic
% Professor is defining feedback gain as S
% u = Sx
% S = -inv(R)*B'*K

clear; clc; close all;

A = [ 1  -2  -2  -1;
      1   1   3   0;
      0  -1   2   0;
      1   3   0  -1 ];
B = [ 1   0;
     -1   0;
      1   2;
      3   2 ];

% State weighting matrix
Q = eye(4); 

% Control input matrix
R = eye(2);   

P = 0;

% For a finite-horizon LQR, we might have a terminal cost matrix K_f. 
% Here, let's just set it to zero for simplicity:
K_f = zeros(4);

t0 = 0;
tf = 5;

% Matlab verification
K = lqr(A,B,Q,R);
disp(K);

% Initial condition for the state
x0_vec = [-5; 2; -1; 3];

% Compare with infinite-horizon LQR for reference
% This is just to see the constant-gain solution from MATLAB's lqr:
K_lqr = lqr(A, B, Q, R);
disp('Infinite-horizon LQR gain (for reference):');
disp(K_lqr);

% ---------------------------------------------------------------
% Solve the "K(t)" Riccati equation backward in time
%    dK/dt = -(A^T K + K A - K B R^{-1} B^T K + Q),  K(tf) = K_f
% ---------------------------------------------------------------
K_final = reshape(K_f, 16,1);  % vectorized version of K_f

% ODE solver options
options = odeset('AbsTol',1e-12,'RelTol',1e-12);

% Integrate from t = 5 down to t = 0
[tSolBackward, KSolBackward] = ode45(@(t,Kvec) kriccati_ode(t, Kvec, A, B, Q, R), ...
                                     [tf t0], K_final, options);

% Flip solution to be forward in time
tSol = flipud(tSolBackward);
KSol = flipud(KSolBackward);

% ---------------------------------------------------------------
% Define interpolation for K(t), then define S(t)
% ---------------------------------------------------------------
K_of_t = @(t) reshape(interp1(tSol, KSol, t, 'linear'), 4,4);

S_of_t = @(t) -inv(R)*(B'*K_of_t(t)+P');

% Evaluate S at t=0 explicitly
S0 = S_of_t(t0);
disp('Optimal feedback gain S(t) at t=0 is:');
disp(S0);

tPlot = linspace(t0, tf, 100);          % 100 sample points
S_of_t_values = zeros(2,4,length(tPlot));  % preallocate for 2x4 S(t)

for i = 1:length(tPlot)
    % Get the 2x4 matrix at time tPlot(i)
    Smat = S_of_t(tPlot(i)); 
    S_of_t_values(:,:,i) = Smat;
end

% Simulate the closed-loop system forward in time:
[tForward, xForward] = ode45(@(t,x) closed_loop_ode_S(t, x, A, B, S_of_t), ...
                             [t0 tf], x0_vec);

toc

% Plot each element of S(t) vs time
figure('Name','Time-Varying Gain S(t)','Color','white','Position',[0 0 1500 1000]);
subplot(2,1,1); hold on; grid minor;
plot(tPlot, squeeze(S_of_t_values(1,1,:)), 'DisplayName','S_{1,1}', 'LineWidth', 2);
plot(tPlot, squeeze(S_of_t_values(1,2,:)), 'DisplayName','S_{1,2}', 'LineWidth', 2);
plot(tPlot, squeeze(S_of_t_values(1,3,:)), 'DisplayName','S_{1,3}', 'LineWidth', 2);
plot(tPlot, squeeze(S_of_t_values(1,4,:)), 'DisplayName','S_{1,4}', 'LineWidth', 2);
legend('Location','best');
xlabel('Time (s)'); ylabel('Row 1 Values');
title('Row 1 of S(t)');

subplot(2,1,2); hold on; grid minor;
plot(tPlot, squeeze(S_of_t_values(2,1,:)), 'DisplayName','S_{2,1}', 'LineWidth', 2);
plot(tPlot, squeeze(S_of_t_values(2,2,:)), 'DisplayName','S_{2,2}', 'LineWidth', 2);
plot(tPlot, squeeze(S_of_t_values(2,3,:)), 'DisplayName','S_{2,3}', 'LineWidth', 2);
plot(tPlot, squeeze(S_of_t_values(2,4,:)), 'DisplayName','S_{2,4}', 'LineWidth', 2);
legend('Location','best');
xlabel('Time (s)'); ylabel('Row 2 Values');
title('Row 2 of S(t)');

set(findall(gcf,'Type','axes'),'FontSize',16);

% Plot the results
figure('Color','white','Position',[1500 0 1500 1000]);
plot(tForward, xForward, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('States');
legend('x_1','x_2','x_3','x_4','Location','Best');
grid on;
title('State trajectories under finite-horizon LQR (Prof. Notation)');
set(gca,'FontSize',20);


for i = 1:length(tForward)
    ti  = tForward(i);
    xi  = xForward(i,:)';       % column vector of the state at time ti
    Sti = S_of_t(ti);          
    ui = Sti * xi;
    
    u_vals(i,:) = ui';          % store as row in the array
end

figure('Color','white','Position',[500 0 1500 1000]);
plot(tForward, u_vals(:,1), 'LineWidth',2, 'DisplayName','u_1(t)');
hold on; grid on;
plot(tForward, u_vals(:,2), 'LineWidth',2, 'DisplayName','u_2(t)');
xlabel('Time (s)');
ylabel('Control Inputs');
legend('Location','best');
title('Control Inputs vs Time');
set(gca,'FontSize',20);


% -------------------------------------------------------------------------
% kriccati_ode: Returns dK/dt for the matrix differential Riccati equation.
%               We pass K as a 16-vector; reshape to 4x4 inside the function.
% -------------------------------------------------------------------------
function dKvec = kriccati_ode(~, Kvec, A, B, Q, R)

    Kmat = reshape(Kvec, 4,4);
    dK = -Kmat*A + Kmat*B*inv(R)*B'*Kmat - Q - A'*Kmat;

    dKvec = dK(:);
end

% -------------------------------------------------------------------------
% closed_loop_ode: Returns dx/dt for the closed-loop system 
%                  x'(t) = (A + B S(t)) x(t).
% -------------------------------------------------------------------------
function dx = closed_loop_ode_S(t, x, A, B, S_of_t)

    St = S_of_t(t);               % Evaluate S at the current time

    dx = (A + B*St)*x;           

end


%% Part c Alternative notation

% Feedback gain = K
% Typically u = -Kx

clear;clc;close all
% Given
A = [ 1  -2  -2  -1;
      1   1   3   0;
      0  -1   2   0;
      1   3   0  -1 ];
B = [ 1   0;
     -1   0;
      1   2;
      3   2 ];
Q = eye(4);       % Q = I_4
R = eye(2);       % R = I_2
P = 0;

t0 = 0;
tf = 5;

x0_vec = [-5 2 -1 3]';   % Initial condition

% Matlab verification
K = lqr(A,B,Q,R);
disp(K);

% Solve the Riccati differential equation backward in time
S_final = zeros(16,1);

% ODE Solver Options
options = odeset('AbsTol',1e-12,'RelTol',1e-12);

% Integrate from t = 5 down to t = 0
[tSolBackward, SSolBackward] = ode45(@(t,Svec) riccati_ode(t, Svec, A, B, Q, R), [tf t0], S_final, options);

% Flip solution to be foward in time
tSol = flipud(tSolBackward);
SSol = flipud(SSolBackward);

% Extract K(t) = R^{-1} B^T S(t) at each time. 
% In particular, we want K(0) = K(tSol(1)).
% For easy simulation, we will interpolate S(t) at arbitrary times.
% Then at each time during the forward simulation, we compute K(t).
%
% Define an interpolation function that returns S(t) in matrix form.

S_of_t = @(t) reshape(interp1(tSol, SSol, t, 'linear'), 4,4);
K_of_t = @(t) R \ (B' * S_of_t(t));

% Evaluate K at t=0 explicitly
K0 = K_of_t(t0);
disp('Optimal feedback gain K at t=0 is:')
disp(K0)

% We'll use ode45 again, but now forward from 0 to 5. 
% The ODE depends on the time-varying gain K(t).
[tForward, xForward] = ode45(@(t,x) closed_loop_ode(t, x, A, B, K_of_t), [t0 tf], x0_vec);

figure('Color','white','Position',[1500 0 1500 1000]);
plot(tForward, xForward, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('States');
legend('x_1','x_2','x_3','x_4','Location','Best');
grid on;
title('State trajectories under finite-horizon LQR control');
ax = gca;
ax.FontSize = 20;


function dSvec = riccati_ode(~, Svec, A, B, Q, R)

    % Reshape S from 16x1 -> 4x4
    S = reshape(Svec, 4, 4);

    % Riccati equation:
    dS = -(A'*S + S*A - S*B*(R\(B'*S)) + Q);
    
    % Return as a 16x1 vector
    dSvec = dS(:);

end

function dx = closed_loop_ode(t, x, A, B, K_of_t)

    % Evaluate K at the current time
    Kt = K_of_t(t);           
    dx = (A - B*Kt) * x;   

end


%% Part d New weighting matrices

% Professor is defining feedback gain as S
% u = Sx
% S = -inv(R)*B'*K

clear; clc; close all;

A = [ 1  -2  -2  -1;
      1   1   3   0;
      0  -1   2   0;
      1   3   0  -1 ];
B = [ 1   0;
     -1   0;
      1   2;
      3   2 ];
Q = eye(4);           
R = 1000*eye(2);      

% For a finite-horizon LQR, we might have a terminal cost matrix K_f. 
% Here, let's just set it to zero for simplicity:
K_f = zeros(4);

t0 = 0;
tf = 5;

% Matlab verification
K = lqr(A,B,Q,R);
disp(K);

% Initial condition for the state
x0_vec = [-5; 2; -1; 3];

% Compare with infinite-horizon LQR for reference
% This is just to see the constant-gain solution from MATLAB's lqr:
K_lqr = lqr(A, B, Q, R);
disp('Infinite-horizon LQR gain (for reference):');
disp(K_lqr);


% Solve the "K(t)" Riccati equation backward in time
% dK/dt = -(A^T K + K A - K B R^{-1} B^T K + Q),  K(tf) = K_f
K_final = reshape(K_f, 16,1);  

% ODE solver options
options = odeset('AbsTol',1e-12,'RelTol',1e-12);

% Integrate from t = 5 down to t = 0
[tSolBackward, KSolBackward] = ode45(@(t,Kvec) kriccati_ode(t, Kvec, A, B, Q, R), ...
                                     [tf t0], K_final, options);

% Flip solution to be forward in time
tSol = flipud(tSolBackward);
KSol = flipud(KSolBackward);


% Define interpolation for K(t), then define S(t)
K_of_t = @(t) reshape(interp1(tSol, KSol, t, 'linear'), 4,4);

S_of_t = @(t) -inv(R)*B'*K_of_t(t);

% Evaluate S at t=0 explicitly
S0 = S_of_t(t0);
disp('Optimal feedback gain S(t) at t=0 is:');
disp(S0);

tPlot = linspace(t0, tf, 100);          % 100 sample points
S_of_t_values = zeros(2,4,length(tPlot));  % preallocate for 2x4 S(t)

for i = 1:length(tPlot)
    % Get the 2x4 matrix at time tPlot(i)
    Smat = S_of_t(tPlot(i)); 
    S_of_t_values(:,:,i) = Smat;
end

% Plot each element of S(t) vs time
figure('Name','Time-Varying Gain S(t)','Color','white','Position',[0 0 1500 1000]);
subplot(2,1,1); hold on; grid minor;
plot(tPlot, squeeze(S_of_t_values(1,1,:)), 'DisplayName','S_{1,1}', 'LineWidth', 2);
plot(tPlot, squeeze(S_of_t_values(1,2,:)), 'DisplayName','S_{1,2}', 'LineWidth', 2);
plot(tPlot, squeeze(S_of_t_values(1,3,:)), 'DisplayName','S_{1,3}', 'LineWidth', 2);
plot(tPlot, squeeze(S_of_t_values(1,4,:)), 'DisplayName','S_{1,4}', 'LineWidth', 2);
legend('Location','best');
xlabel('Time (s)'); ylabel('Row 1 Values');
title('Row 1 of S(t)');

subplot(2,1,2); hold on; grid minor;
plot(tPlot, squeeze(S_of_t_values(2,1,:)), 'DisplayName','S_{2,1}', 'LineWidth', 2);
plot(tPlot, squeeze(S_of_t_values(2,2,:)), 'DisplayName','S_{2,2}', 'LineWidth', 2);
plot(tPlot, squeeze(S_of_t_values(2,3,:)), 'DisplayName','S_{2,3}', 'LineWidth', 2);
plot(tPlot, squeeze(S_of_t_values(2,4,:)), 'DisplayName','S_{2,4}', 'LineWidth', 2);
legend('Location','best');
xlabel('Time (s)'); ylabel('Row 2 Values');
title('Row 2 of S(t)');

set(findall(gcf,'Type','axes'),'FontSize',16);

% Simulate the closed-loop system forward in time:
[tForward, xForward] = ode45(@(t,x) closed_loop_ode_S(t, x, A, B, S_of_t), ...
                             [t0 tf], x0_vec);


% Plot the results
figure('Color','white','Position',[1500 0 1500 1000]);
plot(tForward, xForward, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('States');
legend('x_1','x_2','x_3','x_4','Location','Best');
grid on;
title('State trajectories under finite-horizon LQR (Prof. Notation)');
set(gca,'FontSize',20);


for i = 1:length(tForward)
    ti  = tForward(i);
    xi  = xForward(i,:)';       % column vector of the state at time ti
    Sti = S_of_t(ti);           
    ui = Sti * xi;
    
    u_vals(i,:) = ui';          % store as row in the array
end

figure('Color','white','Position',[500 0 1500 1000]);
plot(tForward, u_vals(:,1), 'LineWidth',2, 'DisplayName','u_1(t)');
hold on; grid on;
plot(tForward, u_vals(:,2), 'LineWidth',2, 'DisplayName','u_2(t)');
xlabel('Time (s)');
ylabel('Control Inputs');
legend('Location','best');
title('Control Inputs vs Time');
set(gca,'FontSize',20);

%% Part e Analytical solution
tic
% Finite-Horizon LQR via Two-Point Boundary-Value Approach
% -------------------------------------------------------------
% We DO NOT solve the Riccati ODE. Instead, we form an augmented
% system in (x, lambda) and use matrix exponentials to enforce:
%   x(0) = x0   (given)
%   lambda(T) = 0   (no terminal cost, P=0).
%
% Then we compute lambda(0) from these boundary conditions,
% and get x(t), lambda(t), and the optimal control 
%   u(t) = -R^(-1) B^T lambda(t).
%
% This solves the same problem as in (c)/(d), and for (f) we
% verify that x*(t), lambda*(t), u*(t) match the Riccati-based
% solution. 
%
% -------------------------------------------------------------

clear; clc; close all;

% 1) Define system matrices (from the problem statement)
A = [ 1  -2  -2  -1;
      1   1   3   0;
      0  -1   2   0;
      1   3   0  -1 ];
B = [ 1   0;
     -1   0;
      1   2;
      3   2 ];

% Dimensions
n = size(A,1);   % n=4
m = size(B,2);   % m=2

% Define cost weighting matrices
Q = eye(n);  
R = 1000*eye(m);  
% P = 0 => no terminal cost, so lambda(T)=0

% Time horizon and initial condition
t0 = 0;
tf = 5;
x0 = [-5; 2; -1; 3];   % given initial state

% Form the augmented system dot(y)= A_aug  where y=[x; lambda]
%   We have from Pontryagin's principle:
%    Control law: u(t) = - R^-1 B^T lambda(t)
%     dot(x)      = A x + B u   = (A x) - (B R^-1 B^T lambda)
%     dot(lambda) = -Q x - A^T lambda

A_aug = [ A , -B*inv(R)*B';
         -Q ,    -A' ];         

% 5) Enforce boundary conditions: x(0)=x0, lambda(tf)=0
%    y(t) = exp(A_aug*t)*z(0). Let z(0)=[x0; lambda(0)].
%    Then y(tf)=[ x(tf); lambda(tf) ]= exp(A_aug*tf)*z(0).
%    We want lambda(tf)=0 => the last n components of z(tf)=0.
%    Partition exp(A_aug*tf) into 4 blocks, each n x n:
%       exp(A_aug*tf) = [ Phi_11,  Phi_12
%                     Phi_21,  Phi_22 ]
%    Then  0 = lambda(tf) = Phi_21*x0 + Phi_22*lambda(0).
%    => lambda(0) = -Phi_22^-1 * Phi_21 * x0

exp_A_aug_tf = expm(A_aug*(tf-t0)); 

Phi_11_tf = exp_A_aug_tf(1:n,       1:n     ); % Upper left quadrant
Phi_12_tf = exp_A_aug_tf(1:n,       n+1:2*n ); % Upper right quadrant
Phi_21_tf = exp_A_aug_tf(n+1:2*n,   1:n     ); % Lower left quadrant
Phi_22_tf = exp_A_aug_tf(n+1:2*n,   n+1:2*n ); % Lower right quadrant

% lam(T) = Phi_21_tf * x0 + Phi_22_tf * lam(0)
% Want to find lam(0)
% lam(T) = 0 from boundary condition (no terminal cost, P = 0)
% Pontryagin's Minimum Principle => lam(T) = 2Px(T) or Px(T)
% => 0 = Phi_21_tf * x0 + Phi_22_tf * lam(0)
lambda0 = - Phi_22_tf \ (Phi_21_tf * x0);

% Sample times, compute x(t), lambda(t), and u(t)
Npoints = 200;
tvec = linspace(t0, tf, Npoints);

for i = 1:Npoints

    ti = tvec(i);

    % Create STMs up to ti
    exp_A_aug_ti = expm(A_aug*(ti - t0));
    
    % Partition again into 4 blocks (n x n)
    Phi_11_ti = exp_A_aug_ti(1:n,       1:n     ); % Upper left quadrant
    Phi_12_ti = exp_A_aug_ti(1:n,       n+1:2*n ); % Upper right quadrant
    Phi_21_ti = exp_A_aug_ti(n+1:2*n,   1:n     ); % Lower left quadrant
    Phi_22_ti = exp_A_aug_ti(n+1:2*n,   n+1:2*n ); % Lower right quadrant
    
    x_ti      = Phi_11_ti*x0 + Phi_12_ti*lambda0;
    lambda_ti = Phi_21_ti*x0 + Phi_22_ti*lambda0;
    
    xvals(:,i)   = x_ti;
    lamvals(:,i) = lambda_ti;
    
    % Control law: u(t) = - R^-1 B^T lambda(t)
    u_ti = - ( R \ (B' * lambda_ti) );
    uvals(:,i) = u_ti;
end

toc

% Plot x(t), lambda(t), u(t)
figure('Color','white','Position',[1000 100 1400 900]);
subplot(3,1,1)
hold on
plot(tvec, xvals, 'LineWidth',2);
xlabel('Time (s)'); ylabel('States x(t)');
legend('x_1','x_2','x_3','x_4','Location','best');
grid minor; title('State Trajectories via STM');


subplot(3,1,2)
hold on
plot(tvec, lamvals, 'LineWidth',2);
xlabel('Time (s)'); ylabel('\lambda(t)');
legend('\lambda_1','\lambda_2','\lambda_3','\lambda_4','Location','best');
grid minor; title('Costate Trajectories');

subplot(3,1,3)
hold on
plot(tvec, uvals, 'LineWidth',2);
xlabel('Time (s)'); ylabel('Control Inputs u(t)');
legend('u_1','u_2','Location','best');
grid minor; title('Control Inputs vs Time');

set(findall(gcf,'Type','axes'),'FontSize',16);

