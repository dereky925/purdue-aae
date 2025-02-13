%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 564 HW1 Double Pendulum simulation
% Author: Derek Yu (Aug 2024)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;close all;

tfinal = 60;
dt = 0.1;
t = 0:dt:tfinal;

name = 'L8';
th1_init = 1.8113;
th2_init = 1.8113;

x0 = [0 0 th1_init 0 th2_init 0];

options = odeset('RelTol', 1E-12, 'AbsTol', 1E-12, 'InitialStep', dt, 'MaxStep', dt);
[t,x] = ode45(@doublePendulumSS, t, x0,options);

%
close all;

a = figure();
sgtitle(name,'FontSize',30)
subplot(3,2,1)
hold on
plot(t,x(:,1))
xlabel('Time [s]','FontSize',20)
ylabel('Cart Position','FontSize',20)
title('Cart Position vs. Time','FontSize',20)

subplot(3,2,2)
hold on
plot(t,x(:,2))
xlabel('Time [s]','FontSize',20)
ylabel('Cart Velocity','FontSize',20)
title('Cart Velocity vs. Time','FontSize',20)

subplot(3,2,3)
hold on
plot(t,x(:,3).*180/pi)
xlabel('Time [s]','FontSize',20)
ylabel('Theta [deg]','FontSize',20)
title('\theta_1 Angular Position vs. Time','FontSize',20)

subplot(3,2,4)
hold on
plot(t,x(:,4).*180/pi)
xlabel('Time [s]','FontSize',20)
ylabel('Angular Velocity [deg/s]','FontSize',20)
title('\theta_1 Angular Velocity vs. Time','FontSize',20)

subplot(3,2,5)
hold on
plot(t,x(:,5).*180/pi)
xlabel('Time [s]','FontSize',20)
ylabel('Theta [deg]','FontSize',20)
title('\theta_2 Angular Position vs. Time','FontSize',20)

subplot(3,2,6)
hold on
plot(t,x(:,6).*180/pi)
xlabel('Time [s]','FontSize',20)
ylabel('Angular Velocity [deg/s]','FontSize',20)
title('\theta_2 Angular Velocity vs. Time','FontSize',20)

a.Position = [100 100 1400 1000];
exportgraphics(a, string(name) + '.pdf', 'ContentType', 'vector');

