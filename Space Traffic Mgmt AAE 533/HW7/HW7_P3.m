% =========================================================================
% 
% Filename:       HW7_P3.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE533 - Space Traffic Management
% Professor:      Dr. Carolin Frueh
% Contact:        cfrueh@purdue.edu
% Assignment:     HW 7
% Semester:       Fall 2024
% 
% Description:
% Propogate mean and covariance
% 
%
% =========================================================================
clear;clc;close all

mu_Earth = 3.9860044e14; %m^3/s^2

% Set initial mean
m0 = [153446.180; 41874155.872; 0;  3066.875;   -11.374;   0];

% Set initial covariance
P0 = [ 6494.080,   -376.139,    0,     0.0160,     -0.494,        0;
       -376.139,    22.560,     0,    -9.883e-4,    0.0286,       0;
           0,         0,     1.205,       0,          0,      -6.071e-5;
        0.0160,   -9.883e-4,    0,     4.437e-8,  -1.212e-6,      0;
        -0.494,     0.0286,     0,    -1.212e-6,   3.762e-5,      0;
           0,         0,     -6.071e-5,   0,          0,      3.390e-9]; 

[d,e] = eig(P0);

% Make covariance semi-positive definite
[V, D] = eig(P0);        % Compute eigenvalues (D) and eigenvectors (V)
D(D < 0) = 0;            % Set negative eigenvalues to zero
P0_psd = V * D * V';     % Reconstruct the matrix

for i = 1:1000
    x0 = m0 + chol(P0_psd)'*randn(length(m0),1);
end
% x0 - m0;

% Mapping matrix
M = [zeros(3,3); eye(3)];
Q = 0.0*eye(3);

num_orbits = 2;
[a,ecc,incl,RAAN,argp,nu,~,~,~] = ijk2keplerian(m0(1:3), m0(4:6)); 
P_coast = num_orbits * (2*pi)/(sqrt(mu_Earth/(a^3))); %Coast time

% Time
t0 = 0;
dt = 10;
tf = P_coast;
tv = t0:dt:tf;

options = odeset('RelTol', 1e-10,'AbsTol',1e-15);

% Code from Lesson 19 
for i = 1:2

    if i == 2
        m0 = x0;
        linespec1 = "k--";
        linespec2 = "k--";
    else
        linespec1 = "r";
        linespec2 = "b";
    end


    % Storage for sigmas and determinant
    sigma_plot = zeros(6,length(tv));
    sigma_plot(:,1) = sqrt(diag(P0_psd));
    detP_plot = zeros(1,length(tv));
    detP_plot(:,1) = det(P0_psd);
    error_plot = zeros(6,length(tv));
   
    
    tkm1 = t0;
    mkm1 = m0;
    Pkm1 = P0_psd;
    
    for k=2:length(tv)
    
        tk = tv(k);
    
        % Propagate from tkm1 to tk
        [~,XX] = ode45(@eom_tbp_ekf,[tkm1 tk], [mkm1;Pkm1(:)], options, mu_Earth,M,Q);
    
        % Determine propogated mean and covariance
        mkm = XX(end,1:6)';
        Pkm = reshape(XX(end,7:42)',6,6);
    
        % Save the sigmas and determinant of P for plotting
        sigma_plot(:,k) = sqrt(diag(Pkm));
        detP_plot(:,k) = det(Pkm);
        error_plot(:,k) = mkm;
    
        % Cycle to the next time
        tkm1 = tk;
        mkm1 = mkm;
        Pkm1 = Pkm;
    
    end

    % Save off mean + covariance of no process noise propogation 
    if i < 2
        mkm0 = mkm;
        Pkm0 = Pkm;
        err0 = error_plot;
    end

    % Compare mean and covariance 
    if i > 1
        disp('Diff:')
        disp(mkm - mkm0)
        disp(Pkm - Pkm0)
       
        disp('Mean Random Draw:')
        disp(mkm)
        disp('Covariance Random Draw:')
        disp(Pkm)
    end

    sigma_plot(1,end)
    
    figA = figure(1);
    subplot(3,1,1)
    hold on
    grid minor
    grid on
    plot(NaN,NaN,'r','LineWidth', 2)
    plot(NaN,NaN,'k--','LineWidth', 2)
    plot(tv/P_coast*2,sigma_plot(1,:)*3,linespec1, 'LineWidth', 2)
    plot(tv/P_coast*2,-sigma_plot(1,:)*3,linespec1, 'LineWidth', 2)
    title('Propogated Position Covariance')
    xlabel('Orbits Periods')
    ylabel('X Position [m]')
    ax = gca;   
    ax.FontSize = 14;
    legend('No Random Draws','Monte Carlo','Location','best','FontSize',18)
    subplot(3,1,2)
    hold on
    grid minor
    grid on
    plot(tv/P_coast*2,sigma_plot(2,:)*3,linespec1, 'LineWidth', 2)
    plot(tv/P_coast*2,-sigma_plot(2,:)*3,linespec1, 'LineWidth', 2)
    xlabel('Orbits Periods')
    ylabel('Y Position [m]')
    ax = gca;   
    ax.FontSize = 14;
    subplot(3,1,3)
    hold on
    grid minor
    grid on
    plot(tv/P_coast*2,sigma_plot(3,:)*3,linespec1, 'LineWidth', 2)
    plot(tv/P_coast*2,-sigma_plot(3,:)*3,linespec1, 'LineWidth', 2)
    xlabel('Orbits Periods')
    ylabel('Z Position [m]')
    ax = gca;   
    ax.FontSize = 14;
    figA.Position = [50 50 1000 800];
    
    
    figB = figure(2);
    subplot(3,1,1)
    hold on
    grid minor
    grid on
    plot(NaN,NaN,'b','LineWidth', 2)
    plot(NaN,NaN,'k--','LineWidth', 2)
    plot(tv/P_coast*2,sigma_plot(4,:)*3,linespec2, 'LineWidth', 2)
    plot(tv/P_coast*2,-sigma_plot(4,:)*3,linespec2, 'LineWidth', 2)
    title('Propogated Velocity Covariance')
    xlabel('Orbits Periods')
    ylabel('X Velocity [m/s]')
    legend('No Random Draws','Monte Carlo','Location','best','FontSize',18)
    ax = gca;   
    ax.FontSize = 14;
    subplot(3,1,2)
    hold on
    grid minor
    grid on
    plot(tv/P_coast*2,sigma_plot(5,:)*3,linespec2, 'LineWidth', 2)
    plot(tv/P_coast*2,-sigma_plot(5,:)*3,linespec2, 'LineWidth', 2)
    subplot(3,1,3)
    xlabel('Orbits Periods')
    ylabel('Y Velocity [m/s]')
    ax = gca;   
    ax.FontSize = 14;
    hold on
    grid minor
    grid on
    plot(tv/P_coast*2,sigma_plot(6,:)*3,linespec2, 'LineWidth', 2)
    plot(tv/P_coast*2,-sigma_plot(6,:)*3,linespec2, 'LineWidth', 2)
    figB.Position = [1000 50 1000 800];
    xlabel('Orbits Periods')
    ylabel('Z Velocity [m/s]')
    ax = gca;   
    ax.FontSize = 14;



end

% Mean Diff plot
figC = figure(3);
subplot(3,1,1)
hold on
grid minor
grid on
plot(NaN,NaN,'r','LineWidth', 2)
plot(tv/P_coast*2,error_plot(1,:) - err0(1,:),'r', 'LineWidth', 2)
title('Propogated Position Mean')
xlabel('Orbits Periods')
ylabel('X Position [m]')
ax = gca;   
ax.FontSize = 14;
% legend('No Random Draws','Location','best','FontSize',18)
subplot(3,1,2)
hold on
grid minor
grid on
plot(tv/P_coast*2,error_plot(2,:) - err0(2,:),'r', 'LineWidth', 2)
xlabel('Orbits Periods')
ylabel('Y Position [m]')
ax = gca;   
ax.FontSize = 14;
subplot(3,1,3)
hold on
grid minor
grid on
plot(tv/P_coast*2,error_plot(3,:) - err0(3,:),'r', 'LineWidth', 2)
xlabel('Orbits Periods')
ylabel('Z Position [m]')
ax = gca;   
ax.FontSize = 14;
figC.Position = [50 50 1000 800];


figD = figure(4);
subplot(3,1,1)
hold on
grid minor
grid on
plot(NaN,NaN,'b','LineWidth', 2)
plot(tv/P_coast*2,error_plot(4,:) - err0(4,:),'b', 'LineWidth', 2)
title('Propogated Velocity Mean')
xlabel('Orbits Periods')
ylabel('X Velocity [m/s]')
% legend('Difference','Location','best','FontSize',18)
ax = gca;   
ax.FontSize = 14;
subplot(3,1,2)
hold on
grid minor
grid on
plot(tv/P_coast*2,error_plot(5,:) - err0(5,:),'b', 'LineWidth', 2)
subplot(3,1,3)
xlabel('Orbits Periods')
ylabel('Y Velocity [m/s]')
ax = gca;   
ax.FontSize = 14;
hold on
grid minor
grid on
plot(tv/P_coast*2,error_plot(6,:) - err0(6,:),'b', 'LineWidth', 2)
figD.Position = [1000 50 1000 800];
xlabel('Orbits Periods')
ylabel('Z Velocity [m/s]')
ax = gca;   
ax.FontSize = 14;


