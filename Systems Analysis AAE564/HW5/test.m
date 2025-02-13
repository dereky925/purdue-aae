clear;clc

m0 = 2; % Mass of the cart
m1 = 1; % Mass of pendulum 1
m2 = 1; % Mass of pendulum 2
l1 = 1; % Length of pendulum 1
l2 = 0.5; % Length of pendulum 2
g = 1; % Gravity

% Set equilibrium angles to pi (pendulums are upside down)
theta_e1 = pi; % Equilibrium configuration for pendulum 1
theta_e2 = pi; % Equilibrium configuration for pendulum 2

M = [m0 + m1 + m2, -m1*l1*cos(theta_e1), -m2*l2*cos(theta_e2);
    -m1*l1*cos(theta_e1), m1*l1^2, 0;
    -m2*l2*cos(theta_e2), 0, m2*l2^2];

K = [0, 0, 0;
     0, -m1*l1*g*cos(theta_e1), 0;
     0, 0, -m2*l2*g*cos(theta_e2)];

% Input matrix (for control input u)
Phi = [1; 0; 0];

% Inverse of M
M_inv = inv(M);

% Define A and B matrices for state-space representation
A_top = [zeros(3), eye(3)];
A_bottom = [M_inv * K, zeros(3)];
A = [A_top; A_bottom];

B = [zeros(3, 1); M_inv * Phi];

% Initial conditions
x0 = [0; 0; pi+0.01; 0; pi+0.01; 0];

tfinal = 60;
dt = 0.1;
t = 0:dt:tfinal;
u = 0; 

[t, x] = ode45(@(t, x) pendulum_dynamics(t, x, A, B, u), t, x0);
%
name = 'L8';

a = figure(1);
hold on
sgtitle(name,'FontSize',30)
subplot(3,2,1)
hold on
grid minor
plot(t,x(:,1))
xlabel('Time [s]','FontSize',20)
ylabel('Cart Position','FontSize',20)
title('Cart Position vs. Time','FontSize',20)

subplot(3,2,2)
hold on
grid minor
plot(t,x(:,2))
xlabel('Time [s]','FontSize',20)
ylabel('Cart Velocity','FontSize',20)
title('Cart Velocity vs. Time','FontSize',20)

subplot(3,2,3)
hold on
grid minor
plot(t,x(:,3).*180/pi)
xlabel('Time [s]','FontSize',20)
ylabel('Theta [deg]','FontSize',20)
title('\theta_1 Angular Position vs. Time','FontSize',20)

subplot(3,2,4)
hold on
grid minor
plot(t,x(:,4).*180/pi)
xlabel('Time [s]','FontSize',20)
ylabel('Angular Velocity [deg/s]','FontSize',20)
title('\theta_1 Angular Velocity vs. Time','FontSize',20)

subplot(3,2,5)
hold on
grid minor
plot(t,x(:,5).*180/pi)
xlabel('Time [s]','FontSize',20)
ylabel('Theta [deg]','FontSize',20)
title('\theta_2 Angular Position vs. Time','FontSize',20)

subplot(3,2,6)
hold on
grid minor
plot(t,x(:,6).*180/pi)
xlabel('Time [s]','FontSize',20)
ylabel('Angular Velocity [deg/s]','FontSize',20)
title('\theta_2 Angular Velocity vs. Time','FontSize',20)

a.Position = [100 100 1400 1000];
exportgraphics(a, string(name) + '.pdf', 'ContentType', 'vector');

function dx = pendulum_dynamics(t, x, A, B, u)
    % State-space dynamics: dx/dt = A*x + B*u
    dx = A * x + B * u;
end