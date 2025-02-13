close all;clear;clc


load('RIOD_5min.mat');
load('RIOD_15min.mat')
load('Gibbs_5min.mat')
load('Gibbs_15min.mat')

size = 300;


a = figure(1);
hold on
grid minor

% [a,ecc,incl,RAAN,argp,nu]

scatter(1,RIOD_5min(1), size,'r', 'filled', 'MarkerEdgeColor', 'k')
scatter(2,RIOD_5min(2), size,'g', 'filled', 'MarkerEdgeColor', 'k')
scatter(3,RIOD_5min(3), size,'b', 'filled', 'MarkerEdgeColor', 'k')
scatter(4,RIOD_5min(4), size,'c', 'filled', 'MarkerEdgeColor', 'k')
scatter(5,RIOD_5min(5), size,'m', 'filled', 'MarkerEdgeColor', 'k')
scatter(6,RIOD_5min(6), size,'y', 'filled', 'MarkerEdgeColor', 'k')

scatter(1,RIOD_15min(1), size,'r', 'd', 'filled', 'MarkerEdgeColor', 'k')
scatter(2,RIOD_15min(2), size,'g', 'd', 'filled', 'MarkerEdgeColor', 'k')
scatter(3,RIOD_15min(3), size,'b', 'd', 'filled', 'MarkerEdgeColor', 'k')
scatter(4,RIOD_15min(4), size,'c', 'd', 'filled', 'MarkerEdgeColor', 'k')
scatter(5,RIOD_15min(5), size,'m', 'd', 'filled', 'MarkerEdgeColor', 'k')
scatter(6,RIOD_15min(6), size,'y', 'd', 'filled', 'MarkerEdgeColor', 'k')

xlabel('Orbital Elements','FontSize',15)
ylabel('% Change in Orbital Elements','FontSize',15)
title('% Change in Restricted IOD Orbital Elements vs Truth','FontSize',25)
legend('Semi-major axis','Eccentricity','Inclination','RAAN','Arg of Perigee','True Anomaly','FontSize',25)
% axis([0.5 3.5 -0.2 0.6])

a.Position = [100 100 1000 1000];

b = figure(2);
hold on
grid minor
scatter(1,Gibbs_5min(1), size,'r', 'filled', 'MarkerEdgeColor', 'k')
scatter(2,Gibbs_5min(2), size,'g', 'filled', 'MarkerEdgeColor', 'k')
scatter(3,Gibbs_5min(3), size,'b', 'filled', 'MarkerEdgeColor', 'k')
scatter(4,Gibbs_5min(4), size,'c', 'filled', 'MarkerEdgeColor', 'k')
scatter(5,Gibbs_5min(5), size,'m', 'filled', 'MarkerEdgeColor', 'k')
scatter(6,Gibbs_5min(6), size,'y', 'filled', 'MarkerEdgeColor', 'k')

scatter(1,Gibbs_15min(1), size,'r', 'd', 'filled', 'MarkerEdgeColor', 'k')
scatter(2,Gibbs_15min(2), size,'g', 'd', 'filled', 'MarkerEdgeColor', 'k')
scatter(3,Gibbs_15min(3), size,'b', 'd', 'filled', 'MarkerEdgeColor', 'k')
scatter(4,Gibbs_15min(4), size,'c', 'd', 'filled', 'MarkerEdgeColor', 'k')
scatter(5,Gibbs_15min(5), size,'m', 'd', 'filled', 'MarkerEdgeColor', 'k')
scatter(6,Gibbs_15min(6), size,'y', 'd', 'filled', 'MarkerEdgeColor', 'k')

xlabel('Orbital Elements','FontSize',15)
ylabel('% Change in Orbital Elements','FontSize',15)
title('% Change in Restricted Gibbs Orbital Elements vs Truth','FontSize',25)
legend('Semi-major axis','Eccentricity','Inclination','RAAN','Arg of Perigee','True Anomaly','FontSize',25)

b.Position = [100 100 1000 1000];
