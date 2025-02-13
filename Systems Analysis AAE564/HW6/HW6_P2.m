clear;clc;close all

function dxdt = odesystem(t, x)

    Omega = 1;
    w_sqrd = 2*Omega^2;

    % dxdt = zeros(4,1);

    dxdt(1,1) = x(3);  % dx1/dt = x3
    dxdt(2,1) = x(4);  % dx2/dt = x4
    
    dxdt(3,1) = (Omega^2 - w_sqrd/2)*x(1) + w_sqrd/2*x(2);  % dx3/dt
    dxdt(4,1) = (Omega^2 - w_sqrd/2)*x(2) + w_sqrd/2*x(1);  % dx4/dt
    
end

% Initial conditions
x0 = [1; -1; 0; 0];  % delta_q1, delta_q2, q1_dot, q2_dot

tspan = [0 20];

options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
[t, x] = ode45(@odesystem, tspan, x0, options);

a = figure;
hold on
grid minor
plot(t, x(:,1), 'r','LineWidth',2)
plot(t, -x(:,2), 'b--','LineWidth',2);
legend('\delta q_1(t)', '-\delta q_2(t)','FontSize',20);
xlabel('Time [sec]','FontSize',20);
ylabel('Displacement','FontSize',20);
title('\delta q_1(t) and -\delta q_2(t) vs Time','FontSize',20);
a.Position = [100 100 1000 800];

b = figure;
hold on
grid minor
plot(x(:,1), x(:,2), 'k','LineWidth',2);
xlabel('\delta q_1(t)','FontSize',20);
ylabel('\delta q_2(t)','FontSize',20);
title('\delta q_2(t) vs \delta q_1(t)','FontSize',20);
b.Position = [1200 100 1000 800];

%%
% Initial conditions
x0 = [1; 1; -1; -1];  % delta_q1, delta_q2, q1_dot, q2_dot

tspan = [0 20];

options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
[t, x] = ode45(@odesystem, tspan, x0, options);

a = figure;
hold on
grid minor
plot(t, x(:,1), 'r','LineWidth',2)
plot(t, x(:,2), 'b--','LineWidth',2);
legend('\delta q_1(t)', '\delta q_2(t)','FontSize',20);
xlabel('Time [sec]','FontSize',20);
ylabel('Displacement','FontSize',20);
title('\delta q_1(t) and -\delta q_2(t) vs Time','FontSize',20);
a.Position = [100 100 1000 800];

b = figure;
hold on
grid minor
plot(x(:,1), x(:,2), 'k','LineWidth',2);
xlabel('\delta q_1(t)','FontSize',20);
ylabel('\delta q_2(t)','FontSize',20);
title('\delta q_2(t) vs \delta q_1(t)','FontSize',20);
b.Position = [1200 100 1000 800];

%%
% Initial conditions
x0 = [1; 1; 1; 1];  % delta_q1, delta_q2, q1_dot, q2_dot

tspan = [0 20];

options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
[t, x] = ode45(@odesystem, tspan, x0, options);

a = figure;
hold on
grid minor
plot(t, x(:,1), 'r','LineWidth',2)
plot(t, x(:,2), 'b--','LineWidth',2);
legend('\delta q_1(t)', '\delta q_2(t)','FontSize',20);
xlabel('Time [sec]','FontSize',20);
ylabel('Displacement','FontSize',20);
title('\delta q_1(t) and -\delta q_2(t) vs Time','FontSize',20);
a.Position = [100 100 1000 800];

b = figure;
hold on
grid minor
plot(x(:,1), x(:,2), 'k','LineWidth',2);
xlabel('\delta q_1(t)','FontSize',20);
ylabel('\delta q_2(t)','FontSize',20);
title('\delta q_2(t) vs \delta q_1(t)','FontSize',20);
b.Position = [1200 100 1000 800];









