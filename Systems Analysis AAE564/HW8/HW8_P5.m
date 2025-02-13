
%% Part A Linearize the system
clear;clc

% m0 = 2;    % Mass of the cart
% m = 1;    % Mass of the pendulum
% l = 1;    % Length of the pendulum
% g = 1; % Gravitational acceleration
% theta_e = 0;
% 
% M = [m0 + m, -m*l*cos(theta_e);
%      -m*l*cos(theta_e), m*l^2];
% 
% % Invert the matrix
% M_inv = inv(M);
% M_ = [0   0 ;
%       0 -m*l*g*cos(theta_e)];
% 
% M_final = M_inv * M_;
% 
% b = inv(M) * [1; 0];
% 
% A_top = [zeros(2), eye(2)];
% A_bottom = [M_final, zeros(2)];
% A = [A_top; A_bottom]
% 
% B = [0; 0; b];

c = 5;
k = 0.5;

% Define symbolic variables
syms m0 m1 m2 l1 l2 theta1_e theta2_e g real
% Define the matrix M 
M = [m0 + m1 + m2, -m1*l1*cos(theta1_e), -m2*l2*cos(theta2_e);
     -m1*l1*cos(theta1_e), m1*l1^2, 0;
     -m2*l2*cos(theta2_e), 0, m2*l2^2];

% Invert the matrix
K = [-k 0 0;
      0 -m1*l1*g*cos(theta1_e) 0;
      0 0 -m2*l2*g*cos(theta2_e)];

C = [-c 0 0; 0 0 0; 0 0 0];

Phi = [1; 0; 0];


A_top = [zeros(3), eye(3)];
A_bottom = [M\K, M\C];
A = [A_top; A_bottom];

B = [0; 0; 0; M\Phi];

A = double ( subs(A, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 1, 1, 0, 0})  );
% A = double ( subs(A, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 1, 1, pi, pi}) );

B = double( subs(B, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 1, 1, 0, 0}));

% P4 E4
% A = double ( subs(A, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 1, 0.5, 0, 0})  )
% B = double( subs(B, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 1, 0.5, 0, 0}))

% A = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1;
%      -0.25 -0.5 -0.5 -2.5 0 0;
%      -0.25 -1.5 -0.5 -2.5 0 0;
%      -0.5  -1 -3 -5 0 0];
% 
% B = [0;0;0;0.5;0.5;1];

C = [1 0 0 0 0 0];
D = [0];

sys = ss(A,B,C,D);

[zeros, poles, ~] = zpkdata(sys);

celldisp(zeros)

celldisp(poles)


%% Part B Simulate non linear system

clear;clc;close all

% Initial conditions
y0 = 0;               % Initial position of the cart
y_dot0 = 0;           % Initial velocity of the cart
theta0 = 0 * pi/180;      % Initial angle of the pendulum (45 degrees)
theta_dot0 = 0;       % Initial angular velocity of the pendulum
initial_conditions = [y0; y_dot0; theta0; theta_dot0];

% Time span
t_span = [0 60]; 

options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
[t, X] = ode45(@(t, x) pendulum_cart_system(t, x), t_span, initial_conditions);

a = figure();
subplot(2,1,1);
hold on
grid minor
plot(t, X(:,1),'LineWidth',2); % Cart position
title('Cart Position Over Time');
xlabel('Time (s)');
ylabel('Position (y)');
ax = gca;
ax.FontSize = 16;
subplot(2,1,2);
hold on
grid minor
plot(t, X(:,3)*180/pi,'LineWidth',2); % Pendulum angle
title('Pendulum Angle Over Time');
xlabel('Time (s)');
ylabel('Angle [degrees]');
a.Position = [100 100 1400 1000];
ax = gca;
ax.FontSize = 16;


%% Part B Simulate Linearized system

clear;clc;close all

% Initial conditions
y0 = 0;               % Initial position of the cart
y_dot0 = 0;           % Initial velocity of the cart
theta0 = 0 * pi/180;      % Initial angle of the pendulum (45 degrees)
theta_dot0 = 0;       % Initial angular velocity of the pendulum
initial_conditions = [y0; y_dot0; theta0; theta_dot0];

% Time span
t_span = [0 60]; % Simulate for 10 seconds

options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
[t, X] = ode45(@(t, x) linearizedSinglePendulum(t, x), t_span, initial_conditions);

a = figure();
subplot(2,1,1);
hold on
grid minor
plot(t, X(:,1),'LineWidth',2); % Cart position
title('Cart Position Over Time');
xlabel('Time (s)');
ylabel('Position (y)');
ax = gca;
ax.FontSize = 16;
subplot(2,1,2);
hold on
grid minor
plot(t, X(:,3)*180/pi,'LineWidth',2); % Pendulum angle
title('Pendulum Angle Over Time');
xlabel('Time (s)');
ylabel('Angle [degrees]');
a.Position = [100 100 1400 1000];
ax = gca;
ax.FontSize = 16;
