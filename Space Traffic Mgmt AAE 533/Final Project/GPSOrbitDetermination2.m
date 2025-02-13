%%
clear;clc;close all

mu = 3.986e14;
RE = 6.371E6; % Earth Radius [m]
c = 299792458; % [m/s]

load('position_POE_ECI.mat');
load('velocity_POE_ECI.mat');

filePath = 'SWOT_VOR_AXVCNE20241121_171715_20241030_225923_20241101_005923.nc';

% info = ncinfo(filePath);  % Inspect the file metadata
% disp(info.Attributes);    % Display global attributes

% Read position data (ECEF position vectors)
position_POE = h5read(filePath, '/position');

% POSITION WAS FOUND TO BE ECEF!! WOW...
% Get metadata for a variable (e.g., 'x')
% xInfo = info.Variables(strcmp({info.Variables.Name}, 'position')); % Replace 'x' with the variable name
% xInfo = info.Variables(strcmp({info.Variables.Name}, 'time'))
% 
% % Display attributes
% disp('Attributes for x:');
% disp(xInfo.Attributes);

% Read velocity data (ECEF velocity vectors)
velocity_POE = h5read(filePath, '/velocity');

% Read time data (UTC) - Epoch: 'seconds since 2000-01-01 00:00:00.0'
time = h5read(filePath, '/time');
epoch = datetime(2000, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
POE_utcTime = epoch + seconds(time);

% Orbit Plot  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
plotA = figure(1);
hold on
imData = imread('2_no_clouds_4k.jpg');  % Load the texture image

[xS, yS, zS] = sphere(50);  % Generate a sphere for the globe
earth_radius = 6378.1370 * 1000;  % Earth radius in meters

xSE = earth_radius * xS;
ySE = earth_radius * yS;
zSE = earth_radius * zS;

% Create the rotation matrices
lon = 1;
Ry = R3(-lon+.6);
Rx = R1(0);

% Combine the rotations (Ry first for longitude, then Rx for latitude)
rotationMatrix = Rx * Ry;

% Apply the rotation to the coordinates
rotated_coords = rotationMatrix * [xSE(:)'; ySE(:)'; zSE(:)'];

% Reshape the rotated coordinates back to match the surface grid dimensions
xSE_rot = reshape(rotated_coords(1, :), size(xSE));
ySE_rot = reshape(rotated_coords(2, :), size(ySE));
zSE_rot = reshape(rotated_coords(3, :), size(zSE));

% Plot the Earth with texture mapping, using the rotated coordinates
surface(xSE_rot, ySE_rot, zSE_rot, 'FaceColor', 'texturemap', 'CData', ...
    flipud(imData), 'EdgeColor', 'none', 'FaceAlpha', 1);

% Adjust axes
axis equal;
view(3);  % 3D view
grid on;
xlabel('Inertial x (m)');
ylabel('Inertial y (m)');
zlabel('Inertial z (m)');

scatter(0,0,1000,'py','filled','MarkerEdgeColor','k')
plot(0, 0, 'r-', 'LineWidth', 5);

% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

gps1 = rinexread('FALK00FLK_R_20243050000_01D_GN.rnx');
gps2 = rinexread('CPVG00CPV_R_20243050000_01D_GN.rnx');
gps3 = rinexread('MCM400ATA_R_20243050000_01D_GN.rnx');
gps4 = rinexread('DYNG00GRC_R_20243050000_01D_GN.rnx');
gps5 = rinexread('HKWS00HKG_R_20243050000_01D_GN.rnx');
gps6 = rinexread('MGO200USA_R_20243050000_01D_GN.rnx');
gps7 = rinexread('PIE100USA_R_20243050000_01D_GN.rnx');

gps1 = gps1.GPS;
gps2 = gps2.GPS;
gps3 = gps3.GPS;
gps4 = gps4.GPS;
gps5 = gps5.GPS;
gps6 = gps6.GPS;
gps7 = gps7.GPS;

combinedTimetable = vertcat(gps1,gps2,gps3,gps4,gps5,gps6,gps7);

% Time might be unique'd automatically, so only check for SatelliteID after
columnsToCheck = combinedTimetable(:, {'SatelliteID'});

% Find unique rows based on the specified columns
[~, uniqueIdx] = unique(columnsToCheck, 'rows', 'stable');

% Filter the original timetable
uniqueTimetable = combinedTimetable(uniqueIdx, :);

GPS_ephemeris = sortrows(uniqueTimetable, {'Time', 'SatelliteID'});

rinexFile = 'SWOT_L1_GPSP_RINEX_1280_20241031T015930_20241031T032050_PIC2_01.rnx';
SWOT_Observation = rinexread(rinexFile);

% obsvec = [38,39,40,41,43];
% gpsvec = [63,65,68,69,72];
obsvec = [38,39,40,41,42,43,45,46,47];
gpsvec = [63,65,68,69,70,72,76,77,78];

[rcvr_pos1, rcvr_vel1, clock_bias1, sat_pos1] = compute_everything(rinexFile, GPS_ephemeris, obsvec, gpsvec);
scatter3(rcvr_pos1(1),rcvr_pos1(2),rcvr_pos1(3),600,'py','filled','MarkerEdgeColor','k')
scatter3(sat_pos1(:,1),sat_pos1(:,2),sat_pos1(:,3),200,'r.','MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.3)


rinexFile = 'SWOT_L1_GPSP_RINEX_1280_20241031T032100_20241031T051310_PIC2_01.rnx';
SWOT_Observation = rinexread(rinexFile);

obsvec = [2466,2467,2468,2470,2471,2472,2473,2474,2475,2476,2477];
gpsvec = [95,96,101,102,104,106,108,111,112,113,116];

[rcvr_pos2, rcvr_vel2, clock_bias2, sat_pos2] = compute_everything(rinexFile, GPS_ephemeris, obsvec, gpsvec);
scatter3(rcvr_pos2(1),rcvr_pos2(2),rcvr_pos2(3),100,'py','filled','MarkerEdgeColor','k')
scatter3(sat_pos2(:,1),sat_pos2(:,2),sat_pos2(:,3),200,'r.','MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.3)


rinexFile = 'SWOT_L1_GPSP_RINEX_1280_20241031T051310_20241031T070940_PIC2_01.rnx';
SWOT_Observation = rinexread(rinexFile);

obsvec = [2637,2638,2639,2640,2641,2642,2644,2645];
gpsvec = [128,132,135,137,140,144,147,149];

[rcvr_pos3, rcvr_vel3, clock_bias3, sat_pos3] = compute_everything(rinexFile, GPS_ephemeris, obsvec, gpsvec);
scatter3(rcvr_pos3(1),rcvr_pos3(2),rcvr_pos3(3),100,'py','filled','MarkerEdgeColor','k')
scatter3(sat_pos3(:,1),sat_pos3(:,2),sat_pos3(:,3),200,'r.','MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.3)


rinexFile = 'SWOT_L1_GPSP_RINEX_1280_20241031T070950_20241031T082730_PIC2_01.rnx';
SWOT_Observation = rinexread(rinexFile);

obsvec = [2894,2895,2896,2897,2898,2899,2900,2901,2903,2904,2905];
gpsvec = [159,165,167,169,170,171,172,176,177,178,184];

[rcvr_pos4, rcvr_vel4, clock_bias4, sat_pos4] = compute_everything(rinexFile, GPS_ephemeris, obsvec, gpsvec);
scatter3(rcvr_pos4(1),rcvr_pos4(2),rcvr_pos4(3),100,'py','filled','MarkerEdgeColor','k')
scatter3(sat_pos4(:,1),sat_pos4(:,2),sat_pos4(:,3),200,'r.','MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.3)


rinexFile = 'SWOT_L1_GPSP_RINEX_1280_20241031T095940_20241031T104450_PIC2_01.rnx';
SWOT_Observation = rinexread(rinexFile);

obsvec = [25,26,27,28,30,31,32,33,34,35];
gpsvec = [196,197,202,203,205,209,213,214,215,219];

[rcvr_pos5, rcvr_vel5, clock_bias5,sat_pos5] = compute_everything(rinexFile, GPS_ephemeris, obsvec, gpsvec);
scatter3(rcvr_pos5(1),rcvr_pos5(2),rcvr_pos5(3),100,'py','filled','MarkerEdgeColor','k')
scatter3(sat_pos5(:,1),sat_pos5(:,2),sat_pos5(:,3),200,'r.','MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.3)

rinexFile = 'SWOT_L1_GPSP_RINEX_1280_20241031T104450_20241031T122320_PIC2_01.rnx';
SWOT_Observation = rinexread(rinexFile);

obsvec = [4329,4330,4331,4332,4333,4334,4335,4336];
gpsvec = [230,234,237,239,242,247,248,251];

[rcvr_pos6, rcvr_vel6, clock_bias6, sat_pos6] = compute_everything(rinexFile, GPS_ephemeris, obsvec, gpsvec);
scatter3(rcvr_pos6(1),rcvr_pos6(2),rcvr_pos6(3),100,'py','filled','MarkerEdgeColor','k')
scatter3(sat_pos6(:,1),sat_pos6(:,2),sat_pos6(:,3),200,'r.','MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.3)


count = 1;
for i = 22:49
    SVID(count) = GPS_ephemeris.SatelliteID(i);
    a = GPS_ephemeris.sqrtA(i)^2;
    e = GPS_ephemeris.Eccentricity(i);
    omega = GPS_ephemeris.omega(i) * 180/pi; % Radians->deg
    RAAN = GPS_ephemeris.OMEGA0(i) * 180/pi; % Radians->deg
    raan_hist(count) = RAAN;

    M0 = GPS_ephemeris.M0(i);

    omega = mod(omega, 360);
    if omega < 0
        omega = omega + 360;
    end

    RAAN = mod(RAAN, 360);
    if RAAN < 0
        RAAN = RAAN + 360;
    end
    
    Crs = GPS_ephemeris.Crs(i);
    Crc = GPS_ephemeris.Crc(i);
    i0 = GPS_ephemeris.i0(i) * 180/pi; % Radians->deg
    Cic = GPS_ephemeris.Cic(i);
    Cis = GPS_ephemeris.Cis(i);
    SVClockBias = GPS_ephemeris.SVClockBias(i);
    
    t = 1:(14*60*60);
    
    options = odeset('RelTol', 1e-8,'AbsTol',1e-10);

    [r_ijk, v_ijk] = keplerian2ijk(a,e,i0,RAAN,omega,0);
    
    [T, Z] = ode45(@two_body_ode,t,[r_ijk, v_ijk],options);
    
    % plot3(Z(:,1), Z(:,2), Z(:,3), 'r', 'LineWidth', 2);

    patch('XData', Z(:,1), 'YData', Z(:,2), 'ZData', Z(:,3), ...
    'EdgeColor', 'red', ...
    'EdgeAlpha', 0.1, ...
    'FaceColor', 'none');

    % scatter3(Z(:,1), Z(:,2), Z(:,3), 1000,'r.','filled',...
    %                 'MarkerEdgeAlpha', 0.5, 'MarkerFaceAlpha', 0.5);

    count = count + 1;

end

% KF ====================================================================

options = odeset('RelTol', 1e-8,'AbsTol',1e-10);

x0ref = [rcvr_pos1 velocity_POE_ECI(1:3,1085)'];

% Timing parameters
t0 = 0.0;
dt = 2*60*60;
tf = 12*60*60;
tv = t0:dt:tf;

% Propogate estimated state and plot KF orbit
t = 1:tf;
[~, Z_prop] = ode45(@two_body_ode,t,x0ref,options);
plotA = figure(1);
% plot3(Z_prop(:,1), Z_prop(:,2), Z_prop(:,3),'m', 'LineWidth', 2);

rng(100);

% Specify initial mean/covariance
m0 = x0ref'; % mean
P0 = diag([50000 ; 50000 ; 50000 ; 2000; 2000; 2000 ].^2); % [m] [m/s]
x0 = x0ref'; % set x0 as initial state for Kalman Filter

% Process noise
Q = blkdiag(eye(3) * 100, eye(3) * 10).^2;

tkm1 = t0; % [s], km1 := k minus 1
xkm1 = x0;
Pkm1 = P0;
mkm1 = [m0(1:3); m0(4:6)];
% mkm1 = [Z(1:3,1); Z(4:6,1)];

STM = eye(6); % Form STM

% Measurements
posVar = 5; % [m]
velVar = 99999; % [m/s]
zk = [ zeros(3,1) rcvr_pos2' rcvr_pos3' rcvr_pos4' rcvr_pos5' rcvr_pos6' rcvr_pos1';...
       zeros(3,1) rcvr_vel2' rcvr_vel3' rcvr_vel4' rcvr_vel5' rcvr_vel6' rcvr_vel1'];
[rows,cols] = size(zk);

% "Truth orbit"
% % % % % % % % truth = Z(tv(2:end),:)';
truth = [position_POE_ECI(:,1085:720:5405); velocity_POE_ECI(:,1085:720:5405)];

% Form Covariance matrix for Kalman filter measurements
for i = 1:cols
    R(:,:,i) = blkdiag(eye(3) * posVar, eye(3) * velVar).^2;
end

sigma_plot = zeros(6,length(tv)-2);
mean_plot = zeros(6,length(tv)-2);
error_plot = zeros(6,length(tv)-2);
residuals_plot = zeros(6,length(tv)-2);

% for i = 2:length(tv)-1 
for i = 2:length(tv)
    % Extract time
    tk = tv(i);

    % Propogate state and STM from tkm1 to tk
    [~, XX] = ode45(@(t,x) two_body_STM(t,x),[tkm1,tk],[mkm1; STM(:)], options);

    xref = XX(end,1:6)'; % Update reference state
    % xref
    Phi = reshape(XX(end,7:end)',6,6);

    % Propogate Covariance and Force Symmetry
    % Pkm1 = Phi*Pkm1 + Pkm1*Phi' + Q;
    Pkm1 = Phi*Pkm1*Phi' + Q; % Slightly different than above, used for discrete filters
    Pkm1 = 0.5*(Pkm1 + Pkm1');
    Pk = Pkm1;
    Pkp = Pk;

    % Form Mapping Matrix
    Ht = eye(6);
    H = Ht*Phi;
    H = Ht;

    % Update mean and covariance with measurements ========================

    % Compute Kalman Gain
    Ck = Pkm1 * H'; % K numerator
    Wk = H * Pkm1 * H' + R(:,:,i); % K denominator
    Kk = Ck/Wk;

    if mod(i,1) == 0 % Use measurements every 4*30 seconds
        % Update mean
        mkp = xref + Kk*(zk(:,i) - xref); % zk = measurements, zref = model
    
        % Update covariance
        Pkp = Pk - Kk*H*Pk; 
        % Pkp = (eye(length(Kk)) - Kk*H) * Pk * (eye(length(Kk)) - Kk*H)' + Kk * R(:,:,i) * Kk';
        Pkp = 0.5 * (Pkp + Pkp'); % Force symmetry

    end

    % Store data for plotting
    sigma_plot(:,i) = real(sqrt(diag(Pkp)));
    mean_plot(:,i) = mkp;   
    error_plot(:,i) = mkp - truth(:,i);
    % error_plot(:,i) = mkp;
    residuals_plot(:,i) = zk(:,i) - xref;

    % Cycle mean, covariance, time
    mkm1 = mkp;
    Pkm1 = Pkp;
    tkm1 = tk;

end

tvec = tv(1:end)/3600;

cov = figure();
subplot(2,1,1)
hold on
grid minor
plot(tvec,sigma_plot(1,:),'r','LineWidth',5)
plot(tvec,sigma_plot(2,:),'b','LineWidth',4)
plot(tvec,sigma_plot(3,:),'g','LineWidth',3)
xlabel('Time [hours]')
ylabel('Position Uncertainty [m]')
title('Estimated Uncertainity vs. Time')
legend('X Position Uncertainty','Y Position Uncertainty','Z Position Uncertainty')
ax = gca();
ax.FontSize = 16;

subplot(2,1,2)
hold on
grid minor
plot(tvec,sigma_plot(4,:),'r','LineWidth',5)
plot(tvec,sigma_plot(5,:),'b','LineWidth',4)
plot(tvec,sigma_plot(6,:),'g','LineWidth',3)
xlabel('Time [hours]')
ylabel('Velocity Uncertainty [m/s]')

legend('X Velocity Uncertainty','Y Velocity Uncertainty','Z Velocity Uncertainty')
cov.Position = [1000 10 600 600];
ax = gca();
ax.FontSize = 16;

residuals = figure();
subplot(2,1,1)
hold on
grid minor
plot(tvec,residuals_plot(1,:),'r','LineWidth',5)
plot(tvec,residuals_plot(2,:),'b','LineWidth',4)
plot(tvec,residuals_plot(3,:),'g','LineWidth',3)
xlabel('Time [hours]')
ylabel('Position Residual [m]')
title('Position & Velocity Residuals vs. Time')
legend('X Position Residual','Y Position Residual','Z Position Residual')
ax = gca();
ax.FontSize = 16;

subplot(2,1,2)
hold on
grid minor
plot(tvec,residuals_plot(4,:),'r','LineWidth',5)
plot(tvec,residuals_plot(5,:),'b','LineWidth',4)
plot(tvec,residuals_plot(6,:),'g','LineWidth',3)
xlabel('Time [hours]')
ylabel('Velocity Residual [m/s]')

legend('X Velocity Residual','Y Velocity Residual','Z Velocity Residual')
residuals.Position = [1500 10 600 600];
ax = gca();
ax.FontSize = 16;


kf = figure();
subplot(2,1,1)
semilogy(tvec,abs(error_plot(1,:)),'r','LineWidth',3)
hold on;
grid minor
semilogy(tvec,abs(error_plot(2,:)),'b','LineWidth',3)
semilogy(tvec,abs(error_plot(3,:)),'g','LineWidth',3)
% plot(tvec,sigma_plot(3,:)*3,'k--','LineWidth',3)
% plot(tvec,sigma_plot(3,:)*-3,'k--','LineWidth',3)

xlabel('Time [hours]')
ylabel('Position [m]')
title('Position Error')
legend('X Position Error','Y Position Error','Z Position Error','Location','best')

ax = gca();
ax.FontSize = 16;

subplot(2,1,2)

semilogy(tvec,abs(error_plot(4,:)),'r','LineWidth',3)
hold on;
grid minor
semilogy(tvec,abs(error_plot(5,:)),'b','LineWidth',3)
semilogy(tvec,abs(error_plot(6,:)),'g','LineWidth',3)
% plot(tvec,sigma_plot(6,:)*3,'k--','LineWidth',3)
% plot(tvec,sigma_plot(6,:)*-3,'k--','LineWidth',3)

xlabel('Time [hours]')
ylabel('Velocity [m/s]')
title('Velocity Error')
legend('X Velocity Error','Y Velocity Error','Z Velocity Error','Location','best')

kf.Position = [10 10 600 600];
ax = gca();
ax.FontSize = 16;

% % Propogate estimated state and plot KF orbit
% t = 1:tf;
% [~, Z_est_KF] = ode45(@two_body_ode,t,[mean_plot(1:3) mean_plot(4:6)],options);
% 
% 
% plotA = figure(1);
% plot3(Z_est_KF(:,1), Z_est_KF(:,2), Z_est_KF(:,3),'b', 'LineWidth', 2);

% Truth orbit
% POE_utcTime_array = [year(POE_utcTime) month(POE_utcTime)...
%                     day(POE_utcTime) hour(POE_utcTime)...
%                     minute(POE_utcTime) second(POE_utcTime)];
% position_POE_ECI = zeros(3,length(position_POE));
% for i = 1:length(position_POE)
%     position_POE_ECI(:,i) = ecef2eci( POE_utcTime_array(i,:), position_POE(:,i));
%     velocity_POE_ECI(:,i) = ecef2eci( POE_utcTime_array(i,:), velocity_POE(:,i));
% end

% for i = 1:length(position_POE_ECI)
%     position_POE_ECI(:,i) = R3(0.61)*position_POE_ECI(:,i);
%     velocity_POE_ECI(:,i) = R3(0.61)*velocity_POE_ECI(:,i);
% end


load('position_POE_ECI.mat');
load('velocity_POE_ECI.mat');
plotA = figure(1);
% scatter3(position_POE_ECI(1,1085),position_POE_ECI(2,1085),position_POE_ECI(3,1085),200,'pr','filled','MarkerEdgeColor','k')


% scatter3(mean_plot(1,:),mean_plot(2,:),mean_plot(3,:),30,'pr','filled','MarkerEdgeColor','k')


% Truth orbit
plot3(position_POE_ECI(1,1085:5405), position_POE_ECI(2,1085:5405), position_POE_ECI(3,1085:5405),'r', 'LineWidth', 2);
plotA.Position = [200 200 700 700];

legend('','GPS Measurements','Truth Orbit')

ax=gca;
ax.FontSize = 20;

%%

% rinexFile = 'SWOT_L1_GPSP_RINEX_1280_20241031T155320_20241031T172120_PIC2_01.rnx';
% SWOT_Observation = rinexread(rinexFile);
% ps = SWOT_Observation.GPS;
% vec = 37:48;
% for i = vec
%     pseudos5 = [pseudos1 ps.C1C(i)];
% end


load('position_POE_ECI.mat');

scatter3(position_POE_ECI(1,1085),position_POE_ECI(2,1085),position_POE_ECI(3,1085),200,'pr','filled','MarkerEdgeColor','k')

% position_POE_ECI = ecef2eci( [2024 10 31 2 00 00], position_POE(:,1085));
% scatter3(position_POE_ECI(1),position_POE_ECI(2),position_POE_ECI(3),500,'pg','filled','MarkerEdgeColor','k')

% Truth orbit
% plot3(position_POE(1,1085:1300), position_POE(2,1085:1300), position_POE(3,1085:1300),'b', 'LineWidth', 2);
plot3(position_POE_ECI(1,1085:1300), position_POE_ECI(2,1085:1300), position_POE_ECI(3,1085:1300),'r', 'LineWidth', 2);


