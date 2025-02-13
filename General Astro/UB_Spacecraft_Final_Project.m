%% MAE 425 SPACECRAFT PROJECT - Group 5
clear;clc;close all

R = 6378.137e3; %Radius of earth
a = 6778e3; % meters
e = 0.02; %eccentricity
i = 29; %degrees
w = 30; %degrees
O = 0; %degress
v = 0; %degrees
mu = 3.9860044e14; %m^3/s^2

% LOITER =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
num_orbits = 10;
P_coast = num_orbits * (2*pi)/(sqrt(mu/(a^3))); %Coast time
t_total = P_coast;
total_days = floor(t_total/86400);
total_hours = floor(mod(t_total,86400)/3600);
total_min = floor(mod(t_total,3600)/60);
total_sec = mod(t_total,60);
fprintf('Number of orbits completed: %.0f\n',num_orbits)
fprintf('Completed on: Feb %.0f 2022, %.0f Hr, %.0f min, %.3f sec\n',1 + total_days, total_hours, total_min, total_sec)
disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')

% TRANSFER ORBIT =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

ra = a*(1+e);
rp = a*(1-e);
r_target = 42164e3;

%Initial Velocity at Perigee
Vp = sqrt(mu*((2/rp)-(1/a)));
% Transfer Orbit
aT1 = (rp + r_target)/2;
% Transfer orbit perigee velocity
VpT = sqrt(mu*((2/rp)-(1/aT1)));
D_v1 = VpT - Vp;
% Transfer orbit apogee velocity
VaT = sqrt(mu*((2/r_target)-(1/aT1)));
%Target orbit circulat velocity
V_circ = sqrt(mu/r_target);
D_v2 =  V_circ - VaT;
%Calculate delta V
P1DV_Total = D_v1 + D_v2;

dv_transfer = 3.8294e+3; %m/s

% Period of Path 1
P = (2*pi)/(sqrt(mu/(aT1^3)));
t_transfer = (1/2)*P;
%Calculate running date
t_total = t_total + t_transfer;
total_days = floor(t_total/86400);
total_hours = floor(mod(t_total,86400)/3600);
total_min = floor(mod(t_total,3600)/60);
total_sec = mod(t_total,60);
fprintf('Transfer orbit delta V: %.4f m/s\n',P1DV_Total)
fprintf('Transfer orbit time: %.4f sec\n',t_transfer)
fprintf('Completed on: Feb %.0f 2022, %.0f Hr, %.0f min, %.3f sec\n',1 + total_days, total_hours, total_min, total_sec)

%Calculate date and JD
hour = 0 + floor(t_transfer/24);
minute = 0 + floor(t_transfer/(24*60));
sec = 0;
year = 2022;
month = 2;
day = 1 + floor(t_transfer/86400); %86400 seconds in a day

J_transfer = 367 * year - floor( 7/4 * (year+floor((month+9)/12))) + floor(275 * month...
    /9) + day + 1721013.5 + ((sec/60+minute)/60+hour)/24;

disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')

% PLANE CHANGE =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
e = 0.02;
i = 29; %degrees
w = 30; %degrees
O = 0; %degress
v = 0; %degrees
a = 42164e3; %a of target sat

V_circ = sqrt(mu/a);
dv_planechange = 2*V_circ*sind(i/2);

w_circ = sqrt(mu/a^3); %angular velocity of circular orbit
P_circ = (2*pi)/(sqrt(mu/(a^3))); %period of circular orbit

t_node = 150*pi/180 / w_circ; %wait time to get to line of nodes

%Calculate running date
t_total = t_node+t_total;
total_days = floor(t_total/86400);
total_hours = floor(mod(t_total,86400)/3600);
total_min = floor(mod(t_total,3600)/60);
total_sec = mod(t_total,60);
fprintf('Wait time till line of nodes: %.3f sec\n',t_node)
fprintf('Plane change delta V: %.3f m/s\n',dv_planechange)
fprintf('Completed on: Feb %.0f 2022, %.0f Hr, %.0f min, %.3f sec\n',1 + total_days, total_hours, total_min, total_sec)
disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')

% Calculate solar align time =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
ra = 42164e3; rp = 42164e3; i = 0; RAAN = 0; w = 0;  nu = 0; P_0 = P_circ; 
e = abs((ra-rp)/(rp+ra)); a = ra/(1+e); r = rp;

co = cosd(RAAN);   so = sind(RAAN);
cw = cosd(w);      sw = sind(w);
ci = cosd(i);      si = sind(i);

e_sat_sun = [0.67038, -0.68081, -0.29513]';

A_P_N = [co*cw-so*sw*ci -co*sw-so*cw*ci so*si;
         so*cw+co*sw*ci -so*sw+co*cw*ci -co*si;
         sw*si cw*si ci];

r_p = [r*cosd(nu) r*sind(nu) 0];
v_p = [-sqrt(mu/(a*(1-e^2)))*sind(nu) sqrt(mu/(a*(1-e^2)))*(e+cosd(nu)) 0];

r_N = A_P_N * r_p';
v_N = A_P_N * v_p';
%integrate to solve for sensor alignment time
options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
[Tout, Y] = ode45(@two_body_ode,[0,P_0],[r_N' v_N'],options);
V_vec = Y(:,4:6);
for i = 1:length(V_vec)
    b2(i,:) = V_vec(i,:) / norm(V_vec(i,:));
    theta_sun(i) = acosd(dot(-b2(i,:),e_sat_sun));
end
[min_ang,t_ind] = min(theta_sun);
t_sun = Tout(t_ind);
inertial_angle = t_sun/max(Tout)*360; %degrees
% figure()
% plot(Tout,theta_sun,'LineWidth',1.5)

pos_tgt = 261.99 + w_circ*180/pi*(P_coast + t_transfer + t_node + t_sun);
pos_chase = 360 + t_sun*w_circ*180/pi; %chase position after waiting for solar

%Calculate running date
t_total = t_total + t_sun;
total_days = floor(t_total/86400);
total_hours = floor(mod(t_total,86400)/3600);
total_min = floor(mod(t_total,3600)/60);
total_sec = mod(t_total,60);
fprintf('Wait time to align solar panels: %.3f sec\n',t_node)
fprintf('Completed on: Feb %.0f 2022, %.0f Hr, %.0f min, %.3f sec\n',1 + total_days, total_hours, total_min, total_sec)
disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')

% RENDEZVOUS =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
ang = pos_chase - pos_tgt;
if ang > 180
    ang = ang-360;
elseif ang < -180
    ang = ang+360;
end
if ang > 0
    ang = ang - 0.5;
else
    ang = ang + 0.5;
end
theta = ang;
n_tgt = sqrt(mu/a^3); %mean motion of target (angular velocity)
w_tgt = n_tgt; %angular velocity of circular orbit
P_tgt = (2*pi)/w_tgt; %period of target (same as period of circle orbit before)
k = 1; %number of orbits to rendezvous
delta_t = (theta*(pi/180))/(k*w_tgt);
Pchase = P_tgt + delta_t; %required chaser period to rendezvous

t_rend = k*Pchase; %rendezvous time

achase = ((mu*Pchase^2)/(4*(pi^2)))^(1/3); %rendezvous semi latus

Vchase = sqrt(mu*((2/a)-(1/achase))); 
Vcirc = sqrt(mu/a);
dv_rend = abs(Vchase - Vcirc)*2; %double delta v to uncircularize then recicularize

%Calculate running date
t_total = t_total+t_rend;
total_days = floor(t_total/86400);
total_hours = floor(mod(t_total,86400)/3600);
total_min = floor(mod(t_total,3600)/60);
total_sec = mod(t_total,60);
% fprintf('Chase: %.2f°\n',pos_chase)
% fprintf('Target: %.2f°\n',pos_tgt)
fprintf('Theta: %.2f°\n',theta)
fprintf('Rendezvous time: %.3f sec\n',t_rend)
fprintf('Rendezvous delta V: %.3f m/s\n',dv_rend)
fprintf('Completed on: Feb %.0f 2022, %.0f Hr, %.0f min, %.3f sec\n',1 + total_days, total_hours, total_min, total_sec)

final_pos_tgt = mod(pos_tgt*w_circ*t_total*180/pi,360);
final_pos_chase = mod(pos_chase, 360);

disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')

dv_total = dv_planechange + dv_transfer + dv_rend; %Total delta V

T_total = P_coast + t_transfer + t_node + t_rend + t_sun + P_tgt; %Total time

fprintf('Total delta V: %.3f km/s\n',dv_total/1000)



%% Plotting =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
close all;clearvars -except P_coast t_transfer t_node t_rend T_total P_tgt t_sun final_pos_chase t_ind V_vec
R = 6378.137e3;
mu = 3.9860044e14;

%        [loiter,    Transfer,    coast,    solar panel,    rendezvous,    MEV,       Galaxy]
ra_vec = [6913560,    6913560,   6913560,    6913560,        6913560,    6913560,   6913560,             6913560,    6913560,    6913560,    6913560,    6913560,    6913560,    6913560,    6913560,    6913560,    6913560,    6913560,    6913560,       ];
rp_vec = [6642440,    6913560,    6913560,    6913560,      6913560,    6913560,   6913560,              6913560,    6913560,    6913560,    6913560,    6913560,    6913560,    6913560,    6913560,    6913560,    6913560,    6913560,    6913560,       ];
i_vec =    [75,          95,        85,         70,               100,          80,          90,              75,         80,         84,         93,         98,         96,         71,         77,         80,         82,         98,         81,       ];
RAAN_vec = [25,           90,         150,          200,               250,          300,         340,       120,        10,        60,           80,        140,        160,        220,        250,        290,        320,        340,        115,       ];
w_vec =    [30,          30,        30,         0,               0,          0,          0                     0,          0,          0,          0,          0,          0,          0,          0,          0,          0,          0,          0,       ]; 
nu_vec =   [0,           0,         180,        0,        final_pos_chase,final_pos_chase, 261.99,             0,          0,          0,          0,          0,          0,          0,          0,          0,          0,          0,          0,       ];
P_0_vec = [P_coast,  t_transfer,   t_node,     t_sun,         t_rend,      P_tgt,   T_total              P_coast,    P_coast,    P_coast,    P_coast,    P_coast,    P_coast,    P_coast,    P_coast,    P_coast,    P_coast,    P_coast,    P_coast,       ];

%preallocate - if not large enough, plot values will appear at 0,0
x = NaN(150000,length(i_vec)*3+3);
phi = NaN(50000,length(i_vec));
lam = NaN(50000,length(i_vec));
az = NaN(50000,length(i_vec));
el = NaN(50000,length(i_vec));

hour = 0;
minute = 0;
sec = 0;
year = 2022;
month = 2;
day = 1;

JD = 367 * year - floor( 7/4 * (year+floor((month+9)/12))) + floor(275 * month...
    /9) + day + 1721013.5 + ((sec/60+minute)/60+hour)/24;

t_length = 1;
for k = 1:length(i_vec)
    
    ra = ra_vec(k);
    rp = rp_vec(k);
    e = abs((ra-rp)/(rp+ra));
    a = ra/(1+e);
    i = i_vec(k);
    RAAN = RAAN_vec(k);
    w = w_vec(k);
    nu = nu_vec(k);
    P_0 = P_0_vec(k);
    
    r = rp;
    % a fix for strange rendezous behavior
    if k == 5
        r = ra;
        e = -e;
    end

    r_p = [r*cosd(nu) r*sind(nu) 0];
    
    co = cosd(RAAN);
    cw = cosd(w);
    ci = cosd(i);
    so = sind(RAAN);
    sw = sind(w);
    si = sind(i);
    
    A_P_N = [co*cw-so*sw*ci -co*sw-so*cw*ci so*si;
             so*cw+co*sw*ci -so*sw+co*cw*ci -co*si;
             sw*si cw*si ci];
    
    v_p = [-sqrt(mu/(a*(1-e^2)))*sind(nu) sqrt(mu/(a*(1-e^2)))*(e+cosd(nu)) 0];
    
    %rotate vectors to inertial frame
    r_N = A_P_N * r_p';
    v_N = A_P_N * v_p';
    
    options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
    % propogate orbits
    [Tout, Y] = ode45(@two_body_ode,[0,P_0],[r_N' v_N'],options); %transpose initial to row vectors
    
    % a fix for last case, galaxy 15. begin orbit at 0 seconds instead of chaining times together
    if k==length(i_vec)
        x(1:length(Y),k*3-2:k*3) = Y(:,1:3);
    else
        %chain orbit times
        x(t_length:t_length + length(Y)-1,k*3-2:k*3) = Y(:,1:3);
        %store time vector and positions
        t_vec(t_length:t_length + length(Y)-1,1) = Tout;
        x_vec(t_length:t_length + length(Y)-1,1:3) = Y(:,1:3);
        t_length = t_length+length(Y);
     
    end
    % =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    r_N_int = Y(:,1:3);
    
    % reset julian date for last case and store time vector
    if k == length(i_vec)
        hour = 0;
        minute = 0;
        sec = 0;
        year = 2022;
        month = 2;
        day = 1;
        
        JD = 367 * year - floor( 7/4 * (year+floor((month+9)/12))) + floor(275 * month...
            /9) + day + 1721013.5 + ((sec/60+minute)/60+hour)/24;
    end
    T = JD + Tout/86400;

    %store case 6 time for movie
    if k == 6
        T_f = T;
    end
    
    %site info
    site_phi = 37.174722;
    site_lam = -5.615833;
    r_E_site = [R*cosd(site_phi)*cosd(site_lam);
                R*cosd(site_phi)*sind(site_lam);
                R*sind(site_phi)];
    
    s_phi = sind(site_phi);
    s_lam = sind(site_lam);
    c_phi = cosd(site_phi);
    c_lam = cosd(site_lam);
    
    A_SEZ = [s_phi*c_lam s_phi*s_lam -c_phi;
             -s_lam c_lam 0;
             c_phi*c_lam c_phi*s_lam s_phi];
    
    
    for i = 1:length(T)
        T_ERA_x = 280.46061837504 + 360.985612288808 * (T(i) - 2451545.0);
        T_ERA(i) = mod(T_ERA_x,360);
    
        rot3 = [cosd(T_ERA(i)) sind(T_ERA(i)) 0;
                -sind(T_ERA(i)) cosd(T_ERA(i)) 0;
                0 0 1];
        
        r_ECEF(i,:) = rot3 * r_N_int(i,:)';
    
        phi(i,k) = asind( r_ECEF(i,3)/norm(r_ECEF(i,:)) );
        lam(i,k) = atan2d(r_ECEF(i,2),r_ECEF(i,1));
        alt(i,k) = norm(r_ECEF(i,:)) - R;
    
        p_E(i,:) = r_ECEF(i,:) - r_E_site';
        p_SEZ(i,:) = A_SEZ * p_E(i,:)';
    
        p(i,k) = norm(p_SEZ(i,:));
        el(i,k) = asind(p_SEZ(i,3)/p(i,k));
        az(i,k) = atan2d(p_SEZ(i,2),-1*p_SEZ(i,1));

        JD = T(end);
        
    end

end

% Find end of mission index
for i = length(lam):-1:1
    if isnan(lam(i,6)) == 0
        ind = i;
        break
    end
end
% Find end of galaxy 15 index
for i = length(x):-1:1
    if isnan(x(i,end)) == 0
        ind_g = i;
        break
    end
end
% Find end of galaxy 15 index
for i = length(lam):-1:1
    if isnan(lam(i,7)) == 0
        ind_g_g = i;
        break
    end
end
% Find end of mev index
for i = length(x):-1:1
    if isnan(x(i,18)) == 0
        ind_g_m = i;
        break
    end
end

% Orbit Plot =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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
% a = plot3(x(:,1), x(:,2), x(:,3), 'm', 'Linewidth', 2); %loiter
% b = plot3(x(:,4), x(:,5), x(:,6), 'r', 'Linewidth', 2); %transfer
% c = plot3(x(:,7), x(:,8), x(:,9), 'b', 'Linewidth', 2); %coast
% d = plot3(x(:,10), x(:,11), x(:,12), 'Color', [0.9290 0.6940 0.1250] , 'Linewidth', 2.5); %solar
% e = plot3(x(:,13), x(:,14), x(:,15), 'g', 'Linewidth', 2.5); %rendezvous
% f = plot3(x(:,16), x(:,17), x(:,18), '--c', 'Linewidth', 2.5); %MEV
% g = plot3(x(:,19), x(:,20), x(:,21), 'k', 'Linewidth', 2); %Galaxy
% 
% plot3(x(:,22), x(:,23), x(:,24), '--k', 'Linewidth', 2); 
% plot3(x(:,25), x(:,26), x(:,27), '--r', 'Linewidth', 2); 
% plot3(x(:,28), x(:,29), x(:,30), '--b', 'Linewidth', 2); 
% plot3(x(:,31), x(:,32), x(:,33), '--m', 'Linewidth', 2); 
% plot3(x(:,34), x(:,35), x(:,36), '--g', 'Linewidth', 2); 
% plot3(x(:,37), x(:,38), x(:,39), '--c', 'Linewidth', 2); 
% plot3(x(:,40), x(:,41), x(:,42), '--y', 'Linewidth', 2); 
% plot3(x(:,43), x(:,44), x(:,45), '--y', 'Linewidth', 2); 
% plot3(x(:,46), x(:,47), x(:,48), '--r', 'Linewidth', 2); 
% plot3(x(:,49), x(:,50), x(:,51), '--b', 'Linewidth', 2); 
% plot3(x(:,52), x(:,53), x(:,54), '--m', 'Linewidth', 2); 
% plot3(x(:,55), x(:,56), x(:,57), '--c', 'Linewidth', 2); 




% plot3([0,100000000],[0,0],[0,0])
% plot3([0,0],[0,100000000],[0,0])
% plot3([0,.67038]*100000000,[0,-.68081]*100000000,[0,-.29513]*100000000)
% plot3([0,V_vec(t_ind,1)],[0,V_vec(t_ind,2)],[0,V_vec(t_ind,3)])
% h = scatter3(x(ind_g,19), x(ind_g,20), x(ind_g,21), 150,'kp'); %Galaxy Check
% i = scatter3(x(ind_g_m,16), x(ind_g_m,17), x(ind_g_m,18), 4500,'kp'); %MEV Check
%legend([a b c d e f g],'Loiter','Transfer Orbit','Coast to line of nodes', 'Align solar panels', 'Rendezvous', 'MEV Final','Galaxy 15')
movegui(a,[0 0]);

%% 

% Ground Trace =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
be = figure();
grid
hold on
a = plot(lam(1:5:end,1),phi(1:5:end,1),'.m','MarkerSize',6); %loiter
b = plot(lam(1:5:end,2),phi(1:5:end,2),'.r','MarkerSize',6); %transfer
c = plot(lam(1:5:end,3),phi(1:5:end,3),'.b','MarkerSize',10); %coast
d = plot(lam(1:5:end,4),phi(1:5:end,4),'.', 'Color',[0.9290 0.6940 0.1250],'MarkerSize',55); %solar
e = plot(lam(1:5:end,5),phi(1:5:end,5),'.g','MarkerSize',10); %rendezvous
f = scatter(lam(1:5:end,6),phi(1:5:end,6),400,'p','MarkerEdgeColor','k','MarkerFaceColor','c'); %MEV
g = scatter(lam(1:5:end,7),phi(1:5:end,7),100,'p','MarkerEdgeColor','k','MarkerFaceColor','k'); %Galaxy
h = scatter(site_lam,site_phi,200,'MarkerEdgeColor','k','MarkerFaceColor','y','Linewidth',2);
% h = scatter(lam(ind,6),phi(ind,6),200,'py','filled','MarkerEdgeColor','k');
legend('Loiter','Transfer Orbit','Coast to line of nodes', 'Align solar panels', 'Rendezvous', 'MEV Final','Galaxy 15','MOSS Site','Location','eastoutside')
xlim([-180,180])
ylim([-90,90])
h = image(xlim,-ylim,imData); 
uistack(h,'bottom')
xlabel('Longitude °')
ylabel('Latitude°')
title('Ground Trace')
% movegui(be,[500 0]);
set(be,'Position',[200 0 1100 500])

for i=1:length(az)
    if az(i)<0
        az(i) = az(i)+360;
    end
end
% az = az + 360*(az<0);

%Observation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c = figure();
grid
hold on
plot(az(1:3:end,1),el(1:3:end,1),'.m') %loiter
plot(az(1:3:end,2),el(1:3:end,2),'.r') %transfer
plot(az(1:3:end,3),el(1:3:end,3),'.b') %coast
plot(az(1:3:end,4),el(1:3:end,4),'.', 'Color',[0.9290 0.6940 0.1250]) %solar
a = plot(az(1:3:end,5),el(1:3:end,5),'.g','MarkerSize',15); %rendezvous
scatter(az(1:3:end,6),el(1:3:end,6),500,'.c') %MEV
scatter(az(1:3:end,7),el(1:3:end,7),150,'.k') %Galaxy
scatter(az(ind,6),el(ind,6),200,'py','filled','MarkerEdgeColor','k')
xlim([0,360])
ylim([10,90])
xlabel('Azimuth °')
ylabel('Elevation °')
title('Sensor Observation')
movegui(c,[1000 0]);

%% Movie
% P_0= T_f(end);
% [Tout, Y] = ode45(@two_body_ode,[0,P_0],[r_N' v_N'],options); %Assuming galaxy is last case
% x(1:length(Y),22:24) = Y(:,1:3);
% %subtract some time off b/c ODE45 time is not linear - makes plotting annoying
% fix = 30;
% 
% close all;clc
% vidObj = VideoWriter('orbit');
% vidObj.FrameRate = 120;
% open(vidObj);
% t = 1:length(x);
% a = view(3);
% msize = 40;
% % for k = 7000:100:length(x)
% % for k = 10000:100:length(x)
% for k = 50:10:14000    
%     k1 = 15;
%     k2 = 30;
%     k3 = 45;
% 
%     scatter3(x(k+fix,19),x(k+fix,20),x(k+fix,21),msize,'k','MarkerFaceColor','k');
% 
%     hold on
%     grid on
%     surface(xSE,ySE,zSE,'FaceColor','texturemap','cdata',flipud(imData),'edgecolor','none');
%     xlabel('X [m]','FontSize',14)
%     ylabel('Y [m]','FontSize',14)
%     zlabel('Z [m]','FontSize',14)
% 
%     scatter3(x(k,1),x(k,2),x(k,3),msize,'magenta','MarkerFaceColor','m');
%     scatter3(x(k,4),x(k,5),x(k,6),msize,'red','MarkerFaceColor','r');
%     scatter3(x(k,7),x(k,8),x(k,9),msize,'blue','MarkerFaceColor','b');
%     scatter3(x(k,10),x(k,11),x(k,12),msize,'Color',[0.9290 0.6940 0.1250],'MarkerFaceColor',[0.9290 0.6940 0.1250]);
%     scatter3(x(k,13),x(k,14),x(k,15),msize,'green','MarkerFaceColor','g');
%     scatter3(x(k,16),x(k,17),x(k,18),msize,'cyan','MarkerFaceColor','c');
%     scatter3(x(k+fix,22),x(k+fix,23),x(k+fix,24),msize,'k','MarkerFaceColor','k');
%     if k >50
%         
%         scatter3(x(k-k1,1),x(k-k1,2),x(k-k1,3),msize,'magenta','MarkerFaceColor','m','MarkerFaceAlpha',.75);
%         scatter3(x(k-k1,4),x(k-k1,5),x(k-k1,6),msize,'red','MarkerFaceColor','r','MarkerFaceAlpha',.75);
%         scatter3(x(k-k1,7),x(k-k1,8),x(k-k1,9),msize,'blue','MarkerFaceColor','b','MarkerFaceAlpha',.75);
%         scatter3(x(k-k1,10),x(k-k1,11),x(k-k1,12),msize,'Color',[0.9290 0.6940 0.1250],'MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerFaceAlpha',.75);
%         scatter3(x(k-k1,13),x(k-k1,14),x(k-k1,15),msize,'green','MarkerFaceColor','g','MarkerFaceAlpha',.75);
%         scatter3(x(k-k1,16),x(k-k1,17),x(k-k1,18),msize,'cyan','MarkerFaceColor','c','MarkerFaceAlpha',.75);
% %         scatter3(x(k-k1,13),x(k-k1,14),x(k-k1,15),msize,'k','MarkerFaceColor','k','MarkerFaceAlpha',.75);
%     
%         scatter3(x(k-k2,1),x(k-k2,2),x(k-k2,3),msize,'magenta','MarkerFaceColor','m','MarkerFaceAlpha',.5);
%         scatter3(x(k-k2,4),x(k-k2,5),x(k-k2,6),msize,'red','MarkerFaceColor','r','MarkerFaceAlpha',.5);
%         scatter3(x(k-k2,7),x(k-k2,8),x(k-k2,9),msize,'blue','MarkerFaceColor','b','MarkerFaceAlpha',.5);
%         scatter3(x(k-k2,10),x(k-k2,11),x(k-k2,12),msize,'Color',[0.9290 0.6940 0.1250],'MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerFaceAlpha',.5);
%         scatter3(x(k-k2,13),x(k-k2,14),x(k-k2,15),msize,'green','MarkerFaceColor','g','MarkerFaceAlpha',.5);
%         scatter3(x(k-k2,16),x(k-k2,17),x(k-k2,18),msize,'cyan','MarkerFaceColor','c','MarkerFaceAlpha',.5);
% %         scatter3(x(k-k2,13),x(k-k2,14),x(k-k2,15),msize,'k','MarkerFaceColor','k','MarkerFaceAlpha',.5);
%     
%         scatter3(x(k-k3,1),x(k-k3,2),x(k-k3,3),msize,'magenta','MarkerFaceColor','m','MarkerFaceAlpha',.25);
%         scatter3(x(k-k3,4),x(k-k3,5),x(k-k3,6),msize,'red','MarkerFaceColor','r','MarkerFaceAlpha',.25);
%         scatter3(x(k-k3,7),x(k-k3,8),x(k-k3,9),msize,'blue','MarkerFaceColor','b','MarkerFaceAlpha',.25);
%         scatter3(x(k-k3,10),x(k-k3,11),x(k-k3,12),msize,'Color',[0.9290 0.6940 0.1250],'MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerFaceAlpha',.25);
%         scatter3(x(k-k3,13),x(k-k3,14),x(k-k3,15),msize,'green','MarkerFaceColor','g','MarkerFaceAlpha',.25);
%         scatter3(x(k-k3,16),x(k-k3,17),x(k-k3,18),msize,'cyan','MarkerFaceColor','c','MarkerFaceAlpha',.25);
% %         scatter3(x(k-k3,13),x(k-k3,14),x(k-k3,15),msize,'k','MarkerFaceColor','k','MarkerFaceAlpha',.25);
%         axis([min(min(x)) max(max(x)) min(min(x)) max(max(x)) min(min(x))/1.5 max(max(x))/1.5])
%     end
% 
%     hold off
%         
% %      view(110,10)      %-30,5 original
%     str1 = ['t = ' num2str(t(k)) 's'];
%     title(str1)
%     currframe = getframe(gcf);
%     writeVideo(vidObj,currframe);
%     drawnow
% 
% end

%% Attitude

% clear;clc;close all
% 
% 
% t = -90;
% C = [cosd(t) 0 -sind(t);
%      0 1 0;
%      sind(t) 0 cosd(t)]
% 
% d = 44.6378;
% B = [cosd(d) sind(d) 0;
%      -sind(d) cosd(d) 0;
%      0 0 1]
% A = C*B
% 
% q4 = 0.5*sqrt(trace(A)+1);
% 
% q123 = 1/(4*q4) * [A(2,3)-A(3,2); A(3,1)-A(1,3); A(1,2)-A(2,1)];
% 
% q = [q123;q4]







