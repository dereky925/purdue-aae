clear; clc; close all

a = 8000; %[km]
ecc = 0.15; 
incl = 80; % [°]
raan = 30; % [°]
argp = 20; % [°]
nu = 0; % [°]

numOrbits = 5;
tfinal = numOrbits * 2*pi*sqrt(a^3/MU); %tfinal
t = [1:tfinal];

[pos_0,vel_0] = keplerian2eci(a,ecc,incl,raan,argp,nu);

options = odeset('AbsTol',1E-12,'RelTol',1E-12);

% regular propagation
[~,Z] = ode45(@(t,x) two_body_ode(t,x),t,[pos_0 ; vel_0],options);

% Define initial STM as a 6x6 identity matrix
Phi0 = eye(6);
X0 = [pos_0; vel_0; Phi0(:)]; % Append STM as a flattened column vector

[t, X] = ode45(@STM_two_body, tspan, X0, options);
% Extract final STM
Phi_final = reshape(X(end, 7:end), 6, 6);


orb_fig = OrbitPlotSetup();

orb_fig.Position = [0 0 1000 1000];
plot3(Z(:,1),Z(:,2),Z(:,3),'r','LineWidth',4)
plot3(states(:,1),states(:,2),states(:,3),'b','LineWidth',4)


