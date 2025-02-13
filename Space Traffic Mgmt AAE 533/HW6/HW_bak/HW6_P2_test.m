%AAE 533 Space Traffic Management
%Homework 6 Problem 2
%Brandon Hughes

%

%clear functions
clc
clear
clf

%Constants
mu_E = 3.986004418e5; % Gravitational parameter for Earth (km^3/s^2)
R = 6371;

%Armstrong Hall Coordinates
GD_Lat = 40.4237;    %[deg] Geodetic Latitude
GD_Long = -86.9212;  %[deg] Geodetic Longitude
h = 185/1000;        %[km] altitude of Purdue
R_ECEF = getECEFfromGD (GD_Lat, GD_Long, h)
%Calculate Goecentric Coordinates
[GC_Lat ,GC_Long] = getGCFromECEF ( R_ECEF );

%TLE for the ISS
tle = {'1 25544U 98067A   24293.74605873  .00028943  00000-0  51103-3 0  9999'
       '2 25544  51.6398  58.4385 0009222  86.0952  55.2425 15.50084362477861'};

% TLE Data: Each line should be a string in a cell array
line1 = tle{1};  % First line of TLE
line2 = tle{2};  % Second line of TLE

%Extracting data from Line 2
a = 400 + R;                         % [km] Altitude             
i = str2double(line2(9:16));         % [deg] Inclination
Omega = str2double(line2(18:25));    % [deg] Right Ascension of Ascending Node
e = str2double(['0.' line2(27:33)]); % Eccentricity
w = str2double(line2(35:42));        % [deg] Argument of Perigee
v = str2double(line2(44:51));        % [deg] Mean Anomaly
mm = str2double(line2(53:63));       % [rev per day] Mean Motion

%Application of Keplers Equation
[r_ijk,v_ijk] = keplerian2ijk(a*1000,e,i,Omega,w,v);
r = r_ijk/1000;
v = v_ijk/1000;

%Integration Initialization
dt = 100;
t0 = 0;
tnext = t0+dt;
options = odeset('RelTol', 1E-12, 'InitialStep', dt, 'MaxStep', dt);
c = 1;

%Keplerian Orbit Integration
state = zeros(2500,6);
epoch = zeros(2500,1);

%Orbit Initialization
x0 = [  r(1) r(2) r(3) v(1) v(2) v(3)];
state(1,:) = x0;
tend = 2500;

%Compute Function
[tt, xt] = ode45('orbit', [t0 tend], x0, options);

%Find satellite RA and dec
[RA , dec] = getGCRAdec ( xt(:, 1:3) ); %deg

%Find postion vector of station
% TLE epoch extraction from line 1
tle_epoch = line1(19:32);  % Extract YYDDD.DDDDDDDD from line 1

% Convert TLE epoch to components
YY = str2double(tle_epoch(1:2));  % Year (last two digits)
DDD = str2double(tle_epoch(3:5)); % Day of the year
fracDay = .74605873;  % Fraction of the day

% Adjust the year (assuming TLE data is from 2000s)
year = 2000 + YY;  
startDate = datetime(year, 1, 1) + days(DDD - 1) + fracDay * days(1);  % Starting date-time

% Convert tt to date vector format
date_vectors = datevec(startDate + seconds(tt));

% Convert date vectors to Julian Date (JD)
JD = juliandate(date_vectors);

%Sidereal Time
[LMSTd , GMSTd ] = getSidereal (JD , GD_Long);

% Convert tt to date vector format
date_vectors = datevec(startDate + seconds(tt));

for i = 1:length(LMSTd)
    R3 = [cosd(GMSTd(i)) sind(GMSTd(i)) 0; -sind(GMSTd(i)) cosd(GMSTd(i)) 0; 0 0 1];
    R3_neg = R3.';
    R_ECI(i, :) = R_ECEF * R3_neg;
end

GC_Lat = ones(size(xt,1),1)*GC_Lat;
% 
% [RA_TC , dec_TC ] = getTCRADec (RA , dec , GC_Lat ,...
%                                 LMSTd , xt(:, 1:3) , R_ECI );
%      
rho = (xt(:,1:3) - R_ECI)';
w = vecnorm(rho(1:2,:));
meas = [atan2d(rho(2,:),rho(1,:)); atan2d(rho(3,:),w)];
R_ECI = R_ECI';

noise = 2;          %arcseconds
% RA = (RA * 3600) + noise;     %arcseconds
% dec = (dec * 3600) + noise;   %arcseconds

plot(dec)

v2 = herrickGibbs(xt(1, 1:3), xt(30, 1:3), xt(60, 1:3), tt(1), tt(30), tt(60));
x2 = xt(30, 1:3);

x02 = [x2 v2];

[~, XX] = ode45('orbit', [tt(30) tt(1)], x02, options);
t1 = tt(1);
x1 = XX(end,:);

error = x1(1:3) - xt(1, 1:3) ;
% End of herrick gibs orbit determination


%Set Up Time
t0 = t1;     %starting time
x0ref = x1';
x0 = zeros(6,1);
sig = 2;

for i = 1:length(xt)
    Rm(:,:,i) = sig^2 * eye(size(meas,1));
end
count = 0;
iter = 4;
for i = 1:iter
    Lam = zeros(6,6);
    lam = zeros(6,1);
    
    tkm1 = t0;
    xref = x0ref;
    Phi = eye(6);
    
    rest = zeros(length(tt),1);
    resm = zeros(length(tt),2);
    for k = 1:length(tt)
        count = count + 1;
        tk = tt(k);
        zk = meas(:,k)*3600;
        Rk = Rm(:,:,k);
        rk = R_ECI(:,k);
        
        Ri = inv(Rk);
        x = [xref(:);Phi(:)];

        if k == 1
            x;
            [tkm1, tk];
        end
        if (tk > t0)
            options = odeset('RelTol',1e-12,'AbsTol',1e-12);
            [~, XX] = ode45(@(t,x) two_body_STM(t,x), [tkm1, tk], x, options);
            xref = XX(end,1:6)';
            Phi = reshape(XX(end,7:end)',[6,6]);
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
        Ht = Ht * 2.062648062470964e+05;
        
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

    % if i == 1
    %     figure()
    % end
    % subplot(4, 2, 2 * (i - 1) + 1);
    % plot(rest, resm(:, 1), '.');
    % xlabel("Time [sec]");
    % ylabel("Error in RA [arcsec]");
    % title("Iteration " + num2str(i) + " - RA Residuals");
    % grid on;
    % 
    % % Plot Dec residuals in the right column
    % subplot(4, 2, 2 * (i - 1) + 2);
    % plot(rest, resm(:, 2), '.');
    % xlabel("Time [sec]");
    % ylabel("Error in Dec [arcsec]");
    % title("Iteration " + num2str(i) + " - Dec Residuals");
    % grid on;
   
end

xhat = x0ref';
P = inv(Lam);


error_lumve = xhat(1:3) - xt(1, 1:3) 
xt(end,1:3);


close all

% Orbit Plot  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
plotA = figure();
hold on
imData = imread('2_no_clouds_4k.jpg');  % Load the texture image

[xS, yS, zS] = sphere(50);  % Generate a sphere for the globe
earth_radius = 6378.1370 * 1000;  % Earth radius in meters

xSE = earth_radius * xS;
ySE = earth_radius * yS;
zSE = earth_radius * zS;

% Create the rotation matrices
Ry = R3(-1+2);

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

R_ECI = R_ECI*1000;
R_ECEF = R_ECEF*1000;
plot3(R_ECI(1,:),R_ECI(2,:),R_ECI(3,:), 'r', 'LineWidth', 4);
scatter3(R_ECEF(1),R_ECEF(2),R_ECEF(3),50)

plotA.Position = [1600 100 1000 1000];