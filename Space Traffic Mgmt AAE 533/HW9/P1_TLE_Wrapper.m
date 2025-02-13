% =========================================================================
% 
% Filename:       P1_TLE_Wrapper.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE533 - Space Traffic Management
% Professor:      Dr. Carolin Frueh
% Contact:        cfrueh@purdue.edu
% Assignment:     HW 9
% Semester:       Fall 2024
% 
% Description:
% EKF for orbit determination. SGP4 code taken from Dr. Frueh
%
%
% =========================================================================

clear;clc; close all;

% LLA of Armstrong Hall
lat = 40.431;
lon = -86.915;
alt = 0;

mu = 3.986e14;
RE = 6.371E6; % Earth Radius [m]


% SGP4 ====================================================================

cc=0; % set line counter

    fid = fopen('STARLINK.txt'); % load the TLE
    
    tline2='gg';
    
    while ischar(tline2)
        cc = cc+1; % counter
        name = fgets(fid);% for the ones with three lines
        tline1 = fgets(fid); % collect first line of two line elements
        tline2 = fgets(fid); % collect second line of two line elements

        if tline2>0 % stop at the end of the file
            % initialize the propagation
            [satrec, startmfe, stopmfe, deltamin] ...
            = twoline2rv(721, tline1, tline2, 'c', 'd');

            epoch = (strsplit(tline1, ' '));
            epoch = char(epoch(4));
            year = 2000 + str2double(epoch(1:2));  % Adds 2000 for 21st century (adjust for other cases)
            day_of_year = str2double(epoch(3:5));
            fractional_day = str2double(['0' epoch(6:end)]);
            date = datetime(year, 1, 1, 'TimeZone', 'UTC') + days(day_of_year - 1) + days(fractional_day);

            period = 60; % [minutes]
        
            % how far shall the TLE be propagated [minutes]
            tsince = period; 

            % extract position and velocity
            for i = 1:period
                tsince = i;
                [satrec, r(i,:), v(i,:)] = sgp4(satrec, tsince); 
            end
        
        end

    end

scatter3(r(1,1)*1000, r(1,2)*1000, r(1,3)*1000, 700, 'py', 'filled','MarkerEdgeColor','k'); 

plot3(r(:,1)*1000, r(:,2)*1000, r(:,3)*1000,'r', 'LineWidth', 2); 


Z_temp = [r,v]*1000;
for i = 1:width(Z_temp)
    x = 1:length(Z_temp(:,i));  % Original indices
    xi = linspace(1, length(Z_temp(:,i)), 60 * length(Z_temp(:,i)));  % Upsampled indices
    Z(:,i) = interp1(x, Z_temp(:,i), xi, 'linear');  % Linear interpolation
end

Tout = 1:period*60;
t = Tout;

% ====================================================================

% Put LLA into vector
lla = [lat,lon,alt];

% Compute ECI XYZ from LLA and UTC time
pos_eci = lla2eci(lla, [2024 11 11 16 47 15]);

mu = 3.9860044e14; %m^3/s^2

utc_datetime = datetime(2024, 11, 11, 16, 47, 15, 'TimeZone', 'UTC');
JD = juliandate(utc_datetime);


options = odeset('RelTol', 1e-10,'AbsTol',1e-15);

Initial_State = Z(1,:);

az = zeros(length(Z),1);
h = zeros(length(Z),1);
ra = zeros(length(Z),1);
dec = zeros(length(Z),1);

% Create measurements
for i = 1:length(Z)
    JD = juliandate(utc_datetime + seconds(Tout(i)));
    pos_eci = lla2eci(lla, datevec(utc_datetime + seconds(Tout(i))));
    pos_eci_plot(i,:) = pos_eci;
    [az(i),h(i),ra(i),dec(i),ra_geoc(i),dec_geoc(i)] = getAzElRaDec(JD,pos_eci,Z(i,1:3));
end

variance = 2; % arcsec
z1 = ra_geoc' + variance*randn();
z2 = dec_geoc' + variance*randn();


% Gibbs' method IOD =======================================================
% Take 3 position measurements from propogation

r_1_Gibbs = [Z(1,1), Z(1,2), Z(1,3)];
r_2_Gibbs = [Z(length(t)/2,1), Z(length(t)/2,2), Z(length(t)/2,3)];
r_3_Gibbs = [Z(end,1), Z(end,2), Z(end,3)];

plane_verification = (r_1_Gibbs/norm(r_1_Gibbs)) * ...
    (cross(r_2_Gibbs,r_3_Gibbs)/norm(cross(r_2_Gibbs,r_3_Gibbs)))';

n = norm(r_1_Gibbs) * cross(r_2_Gibbs,r_3_Gibbs)...
    + norm(r_2_Gibbs)*cross(r_3_Gibbs,r_1_Gibbs) + norm(r_3_Gibbs)...
    * cross(r_1_Gibbs,r_2_Gibbs);
d = cross(r_1_Gibbs,r_2_Gibbs) + cross(r_2_Gibbs,r_3_Gibbs) ...
    + cross(r_3_Gibbs,r_1_Gibbs);
s = r_1_Gibbs*(norm(r_2_Gibbs)-norm(r_3_Gibbs)) ...
    + r_2_Gibbs*(norm(r_3_Gibbs)-norm(r_1_Gibbs)) ...
    + r_3_Gibbs*(norm(r_1_Gibbs)-norm(r_2_Gibbs));

v2 = sqrt( mu/(norm(n)*norm(d)) ) ...
    * ( (cross(d,r_2_Gibbs)/norm(r_2_Gibbs)) + s);


[a_Gibbs,ecc_Gibbs,incl_Gibbs,RAAN_Gibbs,argp_Gibbs,nu_Gibbs,...
    truelon,arglat,lonper] = ijk2keplerian(r_2_Gibbs, v2);

[Tout_Gibbs, Z_Gibbs] = ode45(@two_body_ode,t,[r_2_Gibbs v2],options);

% Propogate back to t1
[Tout_IOD, Z_IOD] = ode45(@two_body_ode,t/2,[r_2_Gibbs -v2],options);



% LUMVE
t1 = Tout(1);
x1 = Z_IOD(end,:);

t0 = t1;
x0ref = x1';
x0 = zeros(6,1);

iter = 4;

T_meas = [1, length(Tout)/2, Tout(end)];
% T_meas = [1:length(Tout)];
% Tm is the time indexed at the 3 measurements
Tm = Tout(T_meas);


% Form Covariance matrix
for i = 1:length(Tm)
    Rm(:,:,i) = 2^2 * eye(2);
end

Zm = [z1(T_meas)'.*3600 ; z2(T_meas)'.*3600];
% Zm = meas(:,T_meas).*3600

obs1time = utc_datetime;
obs2time = utc_datetime + seconds(length(Tout)/2);
obs3time = utc_datetime + seconds(length(Tout));

% Vector of Station ECI coords at measurement times
Rt = [lla2eci(lla, datevec(obs1time))' lla2eci(lla, datevec(obs2time))'...
    lla2eci(lla, datevec(obs3time))'];

% Rt = [];
% for i = 1:length(Tm)
%     obstime = utc_datetime+seconds(Tm(i));
%     Rt = [Rt,lla2eci(lla, datevec(obstime))'];
% end
options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
count = 0;
for loop = 1:iter

    % initialize Lambda and lambda to 0, no prior info
    lam = zeros(6,1);
    Lam = zeros(6,6);

    %initialize time, reference state, and STM
    tkm1 = t0;
    xref = x0ref/1000;
    Phi = eye(6);

    % residuals?
    rest = zeros(length(Tm),1);
    resm = zeros(length(Tm),2);

    for k = 1:length(Tm)
        count = count + 1;
        % extract the time, measurement, covariance and station position
        % for the kth observation
        tk = T_meas(k); % Observation Time
        zk = Zm(:,k);   % Measurement value, RA dec
        Rk = Rm(:,:,k); % Covariance
        rk = Rt(:,k)/1000; % Station ECI position

        % determine LUMVE weighting matrix
        Ri = inv(Rk);

        % Propogate the reference state and STM but only if we're past t0
        if (tk>t0)
            [~, XX] = ode45(@(t,x) two_body_STM_km(t,x),[tkm1,tk],[xref; Phi(:)],options);
            xref = XX(end,1:6)'; % Update reference state
            Phi = reshape(XX(end,7:end)',6,6);
        end

        rosi = xref(1:3) - rk;
        x = rosi(1);
        y = rosi(2);
        z = rosi(3);
        wsq = x*x + y*y;
        w = sqrt(wsq);
        rhosq = wsq + z*z;
        zref = [atan2d(y,x); atan2d(z,w)];
        zref = zref*3600;

        Ht = [-y/wsq, x/wsq, 0, 0, 0, 0;
              -x*z/(w*rhosq), -y*z/(w*rhosq), w/rhosq, 0, 0, 0];
        Ht = Ht * 2.0626e+05;

        H = Ht*Phi;

        lam = lam + H'*Ri*(zk - zref);
        Lam = Lam + H'*Ri*H;

        rest(k) = tk;
        resm(k,:) = zk - zref;

        tkm1 = tk;

    end

    delx = Lam\lam;

    x0ref = x0ref + delx;
    x0 = x0 - delx;

end

xhat = x0ref';
P = inv(Lam);

error_lumve = xhat(1:3) - Z(1, 1:3)

% Propogate estimated state and plot LUMVE orbit
[Tout_est, Z_est_LUMVE] = ode45(@two_body_ode,t,[xhat(1:3) -xhat(4:6)],options);


%% Kalman Filter ===========================================================
close all

x0ref = 1.0e+06 * [-6.6590  -1.8960 0.1807 0.0014 -0.0043 0.0061 ];

rng(100);

% Timing parameters
t0 = 0.0;
dt = 1.0;
tf = 3600.0;
tv = t0:dt:tf;

% Specify initial mean/covariance
m0 = x0ref'; % mean
P0 = diag([50000 ; 50000 ; 50000 ; 2000; 2000; 2000 ].^2); % [m] [m/s]
x0 = x0ref'; % set x0 as initial state for Kalman Filter

% Process noise
Q = blkdiag(eye(3) * 10, eye(3) * 10).^2;

tkm1 = t0; % [s], km1 := k minus 1
xkm1 = x0;
Pkm1 = P0;
mkm1 = [m0(1:3); m0(4:6)];
% mkm1 = [Z(1:3,1); Z(4:6,1)];

STM = eye(6); % Form STM

% Generate measurements with the same propagator as the kalman filter
% If you use SGP4 to make the measurements, oscillations & large residuals
% will be observed
t = tv;
[~, Z] = ode45(@two_body_ode,t,[x0ref(1:3) x0ref(4:6)],options);

% Measurements
posVar = 6; % [m]
velVar = 1; % [m/s]
% zk = [-z1(tv(2:end))'.*3600 ; z2(tv(2:end))'.*3600]; % Optical
zk = [Z(tv(2:end),:)]'; % Position + velocity measurements
% Add noise
zk = [zk(1:3,:) + posVar*randn(); zk(4:6,:) + velVar*randn()];

% "Truth orbit"
truth = Z(tv(2:end),:)';

% mkm1 = [truth(1:3,1);truth(4:6,1)];

% Compute station coords at every measurement
rk = pos_eci_plot(tv(2:end),:);

% Form Covariance matrix for Kalman filter measurements
for i = 1:length(Z)
    R(:,:,i) = blkdiag(eye(3) * posVar, eye(3) * velVar).^2;
end

sigma_plot = zeros(6,length(tv)-2);
mean_plot = zeros(6,length(tv)-2);
error_plot = zeros(6,length(tv)-2);
residuals_plot = zeros(6,length(tv)-2);

for i = 2:length(tv)-1 
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
    sigma_plot(:,i-1) = real(sqrt(diag(Pkp)));
    mean_plot(:,i-1) = mkp;   
    error_plot(:,i-1) = mkp - truth(:,i);
    residuals_plot(:,i-1) = zk(:,i) - xref;

    % Cycle mean, covariance, time
    mkm1 = mkp;
    Pkm1 = Pkp;
    tkm1 = tk;

end

close all

cov = figure();
subplot(2,1,1)
hold on
grid minor
plot(tv(2:end-1),sigma_plot(1,:),'r','LineWidth',5)
plot(tv(2:end-1),sigma_plot(2,:),'b','LineWidth',4)
plot(tv(2:end-1),sigma_plot(3,:),'g','LineWidth',3)
xlabel('Time [s]')
ylabel('Position Uncertainty [m]')
title('Estimated Uncertainity vs. Time')
legend('X Position Uncertainty','Y Position Uncertainty','Z Position Uncertainty')
ax = gca();
ax.FontSize = 16;

subplot(2,1,2)
hold on
grid minor
plot(tv(2:end-1),sigma_plot(4,:),'r','LineWidth',5)
plot(tv(2:end-1),sigma_plot(5,:),'b','LineWidth',4)
plot(tv(2:end-1),sigma_plot(6,:),'g','LineWidth',3)
xlabel('Time [s]')
ylabel('Velocity Uncertainty [m/s]')

legend('X Velocity Uncertainty','Y Velocity Uncertainty','Z Velocity Uncertainty')
cov.Position = [1000 10 1000 1000];
ax = gca();
ax.FontSize = 16;

residuals = figure();
subplot(2,1,1)
hold on
grid minor
plot(tv(2:end-1),residuals_plot(1,:),'r','LineWidth',5)
plot(tv(2:end-1),residuals_plot(2,:),'b','LineWidth',4)
plot(tv(2:end-1),residuals_plot(3,:),'g','LineWidth',3)
xlabel('Time [s]')
ylabel('Position Residual [m]')
title('Position & Velocity Residuals vs. Time')
legend('X Position Residual','Y Position Residual','Z Position Residual')
ax = gca();
ax.FontSize = 16;

subplot(2,1,2)
hold on
grid minor
plot(tv(2:end-1),residuals_plot(4,:),'r','LineWidth',5)
plot(tv(2:end-1),residuals_plot(5,:),'b','LineWidth',4)
plot(tv(2:end-1),residuals_plot(6,:),'g','LineWidth',3)
xlabel('Time [s]')
ylabel('Velocity Residual [m/s]')

legend('X Velocity Residual','Y Velocity Residual','Z Velocity Residual')
residuals.Position = [1500 10 1000 1000];
ax = gca();
ax.FontSize = 16;


kf = figure();
subplot(2,1,1)
hold on;
grid minor
plot(tv(2:end-1),error_plot(1,:),'r')
% plot(sigma_plot(1,:)*3,'r--')
% plot(sigma_plot(1,:)*-3,'r--')

plot(tv(2:end-1),error_plot(2,:),'b')
% plot(sigma_plot(2,:)*3,'b--')
% plot(sigma_plot(2,:)*-3,'b--')

plot(tv(2:end-1),error_plot(3,:),'g')
plot(tv(2:end-1),sigma_plot(3,:)*3,'k--','LineWidth',3)
plot(tv(2:end-1),sigma_plot(3,:)*-3,'k--','LineWidth',3)

xlabel('Time [s]')
ylabel('Position [m]')
title('Position Error and 3\sigma Bounds')
legend('X Position Error','Y Position Error','Z Position Error','3\sigma Bound','Location','best')

ax = gca();
ax.FontSize = 16;

subplot(2,1,2)
hold on;
grid minor
plot(tv(2:end-1),error_plot(4,:),'r')
% plot(sigma_plot(4,:)*3,'r--')
% plot(sigma_plot(4,:)*-3,'r--')

plot(tv(2:end-1),error_plot(5,:),'b')
% plot(sigma_plot(5,:)*3,'b--')
% plot(sigma_plot(5,:)*-3,'b--')

plot(tv(2:end-1),error_plot(6,:),'g')
plot(tv(2:end-1),sigma_plot(6,:)*3,'k--','LineWidth',3)
plot(tv(2:end-1),sigma_plot(6,:)*-3,'k--','LineWidth',3)

xlabel('Time [s]')
ylabel('Velocity [m/s]')
title('Velocity Error and 3\sigma Bounds')
legend('X Velocity Error','Y Velocity Error','Z Velocity Error','3\sigma Bound','Location','best')

kf.Position = [10 10 1000 1000];
ax = gca();
ax.FontSize = 16;

% Propogate estimated state and plot KF orbit
[~, Z_est_KF] = ode45(@two_body_ode,t,[mean_plot(1:3) mean_plot(4:6)],options);


% Orbit Plot  ===========================================================
plotA = figure();
hold on
imData = imread('2_no_clouds_4k.jpg');  % Load the texture image

[xS, yS, zS] = sphere(50);  % Generate a sphere for the globe
earth_radius = 6378.1370 * 1000;  % Earth radius in meters

xSE = earth_radius * xS;
ySE = earth_radius * yS;
zSE = earth_radius * zS;

% Create the rotation matrices
Ry = R3(-lon+1.9);

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
    flipud(imData), 'EdgeColor', 'none');

% Adjust axes
axis equal;
view(3);  % 3D view
grid on;
xlabel('Inertial x (m)');
ylabel('Inertial y (m)');
zlabel('Inertial z (m)');


% Truth orbit
plot3(Z(:,1), Z(:,2), Z(:,3), 'r', 'LineWidth', 3);
% plot3(zk(1,:), zk(2,:), zk(3,:), 'b', 'LineWidth', 4);
% plot3(Z_IOD(:,1), Z_IOD(:,2), Z_IOD(:,3), 'b--', 'LineWidth', 4);

%Plot estimated LUMVE orbit
plot3(Z_est_LUMVE(:,1), Z_est_LUMVE(:,2), Z_est_LUMVE(:,3), 'b:', 'LineWidth', 4);

%Plot estimated KF orbit
plot3(mean_plot(1,1:end-1), mean_plot(2,1:end-1), mean_plot(3,1:end-1), 'g:', 'LineWidth', 4);

plot3(pos_eci_plot(:,1),pos_eci_plot(:,2),pos_eci_plot(:,3), 'm', 'LineWidth', 4);

legend('','Truth Orbit','LUMVE Orbit','EKF Orbit','Station Position')

plotA.Position = [1600 100 1000 1000];
ax = gca();
ax.FontSize = 16;




