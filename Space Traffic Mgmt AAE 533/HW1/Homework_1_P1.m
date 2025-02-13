clear;clc;close all

mu = 3.9860044e14; %m^3/s^2

% LEO Parameters =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

a = 6771000; % Semi-major axis [m] LEO
e = 0.01; % Eccentricity
i = 30; % Inclination [deg]
RAAN = 200; % Right Ascension of Ascending Node [deg]
w = 0; % Argument of Perigee [deg]
v = 0; % True Anomaly [deg]

num_orbits = 2;
P_coast = num_orbits * (2*pi)/(sqrt(mu/(a^3))); %Coast time
t = 1:1:P_coast;

[r_ijk,v_ijk] = keplerian2ijk(a,e,i,RAAN,w,v);

options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
[Tout, Y] = ode45(@two_body_ode,t,[r_ijk v_ijk],options);

% GEO Parameters =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

a = 42378137; % Semi-major axis [m] GEO
e = 0.01; % Eccentricity
i = 1; % Inclination [deg]
RAAN = 0; % Right Ascension of Ascending Node [deg]
w = 0; % Argument of Perigee [deg]
v = 0; % True Anomaly [deg]

num_orbits = 2;
P_coast = num_orbits * (2*pi)/(sqrt(mu/(a^3))); %Coast time
t = 1:1:P_coast;

[r_ijk,v_ijk] = keplerian2ijk(a,e,i,RAAN,w,v);

options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
[Tout, Z] = ode45(@two_body_ode,t,[r_ijk v_ijk],options);


% Orbit Plot =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
imData = imread('2_no_clouds_4k.jpg');
[xS,yS,zS] = sphere(50);
earth_radius = 6378.1370*1000;  % meters
xSE = earth_radius*xS;
ySE = earth_radius*yS;
zSE = earth_radius*zS;
surface(xSE,ySE,zSE,'FaceColor','texturemap','cdata',flipud(imData),'edgecolor','none');
axis equal
a = view(3);
grid on
xlabel('Inertial x (m)')
ylabel('Inertial y (m)')
zlabel('Inertial z (m)')
ch = get(gca,'children');
hold on 
plot3(Y(:,1), Y(:,2), Y(:,3), 'r','Linewidth', 4);
plot3(Z(:,1), Z(:,2), Z(:,3), 'b','Linewidth', 4);
legend('','LEO','GEO','FontSize',15)


%% Calc Error

clear;clc;close all

mu = 3.9860044e14; %m^3/s^2

% LEO

a = 6771000; % Semi-major axis [m] LEO
e = 0.01; % Eccentricity
i = 30; % Inclination [deg]
RAAN = 200; % Right Ascension of Ascending Node [deg]
w = 90; % Argument of Perigee [deg]
v = 90; % True Anomaly [deg]

num_orbits = 2;
P_coast = num_orbits * (2*pi)/(sqrt(mu/(a^3))); %Coast time
t = 1:1:P_coast;

[r_ijk,v_ijk] = keplerian2ijk(a,e,i,RAAN,w,v);

relTolList = [1e-4, 1e-7, 1e-10];
absTolList = [1e-5, 1e-10, 1e-15];

j=1;
for k=1:3
    options = odeset('RelTol', relTolList(k),'AbsTol',absTolList(k));
    [Tout, X(:,j:k*6)] = ode45(@two_body_ode,t,[r_ijk v_ijk],options);
    
    j=j+6; % ode45 X output increment
end


leoFinalPosList_ijk = X(end,:);

j=1;
for k=1:3
    [a,ecc,incl,RAAN,argp,nu,truelon,arglat,lonper] = ...
    ijk2keplerian(leoFinalPosList_ijk(:,j:j+2), leoFinalPosList_ijk(:,j+3:j+5));
    
    leoFinalPosList_kep(k,:) = [a,ecc,incl,RAAN,argp,nu,truelon,arglat,lonper];
    j=j+6;
end


% Redefine Overwritten Original Orbital Elements
a = 6771000; % Semi-major axis [m] LEO
e = 0.01; % Eccentricity
i = 30; % Inclination [deg]
RAAN = 200; % Right Ascension of Ascending Node [deg]
w = 90; % Argument of Perigee [deg]
v = 90; % True Anomaly [deg]

% Calculate % change
leoPercentChangeVector(:,1) =  (leoFinalPosList_kep(:,1) - a)./a .* 100;
leoPercentChangeVector(:,2) =  (leoFinalPosList_kep(:,2) - e)./e .* 100;
leoPercentChangeVector(:,3) =  (leoFinalPosList_kep(:,3) - i)./i .* 100;
leoPercentChangeVector(:,4) =  (leoFinalPosList_kep(:,4) - RAAN)./RAAN .* 100;
leoPercentChangeVector(:,5) =  (leoFinalPosList_kep(:,5) - w)./w .* 100;
leoPercentChangeVector(:,6) =  (leoFinalPosList_kep(:,6) - v)./v .* 100;

x = 1:3;
a = figure();
hold on
grid minor
scatter(NaN,NaN,400,'r','filled','MarkerEdgeColor','k','LineWidth',1.5)
scatter(NaN,NaN,400,'b','filled','MarkerEdgeColor','k','LineWidth',1.5)
scatter(NaN,NaN,400,'g','filled','MarkerEdgeColor','k','LineWidth',1.5)

scatter(1,leoPercentChangeVector(1,:),200,'r','filled','MarkerEdgeColor','k','LineWidth',1.5)

scatter(2,leoPercentChangeVector(2,1),280,'b','filled','MarkerEdgeColor','k','LineWidth',1.5)
scatter(2,leoPercentChangeVector(2,2),260,'b','filled','MarkerEdgeColor','k','LineWidth',1.5)
scatter(2,leoPercentChangeVector(2,3),240,'b','filled','MarkerEdgeColor','k','LineWidth',1.5)
scatter(2,leoPercentChangeVector(2,4),220,'b','filled','MarkerEdgeColor','k','LineWidth',1.5)
scatter(2,leoPercentChangeVector(2,5),200,'b','filled','MarkerEdgeColor','k','LineWidth',1.5)
scatter(2,leoPercentChangeVector(2,6),200,'b','filled','MarkerEdgeColor','k','LineWidth',1.5)

scatter(3,leoPercentChangeVector(3,1),280,'g','filled','MarkerEdgeColor','k','LineWidth',1.5)
scatter(3,leoPercentChangeVector(3,2),260,'g','filled','MarkerEdgeColor','k','LineWidth',1.5)
scatter(3,leoPercentChangeVector(3,3),240,'g','filled','MarkerEdgeColor','k','LineWidth',1.5)
scatter(3,leoPercentChangeVector(3,4),220,'g','filled','MarkerEdgeColor','k','LineWidth',1.5)
scatter(3,leoPercentChangeVector(3,5),200,'g','filled','MarkerEdgeColor','k','LineWidth',1.5)
scatter(3,leoPercentChangeVector(3,6),200,'g','filled','MarkerEdgeColor','k','LineWidth',1.5)
a.Position = [100 100 1000 1000];

xlabel('Decreasing Abs and Rel Tol','FontSize',15)
ylabel('% Change in Orbital Elements','FontSize',15)
title('% Change in Orbital Elements after 2 Orbits','FontSize',15)
legend('RelTol: 1e-4, AbsTol:1e-5', 'RelTol: 1e-7, AbsTol:1e-10', 'RelTol: 1e-10, AbsTol:1e-15','FontSize',25)
axis([0.5 3.5 -0.2 0.6])

