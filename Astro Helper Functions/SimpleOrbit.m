clear;clc;close all

[r_eci, v_eci, a, e, i, RAAN, w, nu] = getISSState();

num_orbits = 50;
T = num_orbits * (2*pi)*(sqrt(a^3/MU)); % Orbit period
t = 1:T;

options = odeset('RelTol',1E-12,'AbsTol',1E-12);

% Propagate orbital elements
[tvec, pos] = ode45(@(t,x) two_body_ode(t,x),t, [r_eci,v_eci], options);

r = vecnorm(pos,2,2);

orb = OrbitPlotSetup;
orb.Position = [0 0 1000 1000];
plot3(pos(:,1),pos(:,2),pos(:,3),'r','LineWidth',4);

figure('Color','white','Position',[500 100 1500 1000])
hold on;grid minor
plot(tvec/3600,r)

