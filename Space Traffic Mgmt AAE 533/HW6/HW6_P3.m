clear;clc;close all

obs = load('/Users/derekyu/Mac Documents/MATLAB/Purdue/Space Traffic Mgmt AAE 533/HW6/observations.mat');

% LLA of new mexico
lat = 32.54;
lon = -105;
alt = 7300;

% Put LLA into vector
lla = [lat,lon,alt];

% Compute ECI XYZ from LLA and UTC time, Sep 6, 0:30 AM 2024 [m]
pos_eci = lla2eci(lla, [2024 10 16 12 0 0]);

mu = 3.9860044e14; %m^3/s^2

utc_datetime = datetime(2024, 10, 16, 12, 0, 0, 'TimeZone', 'UTC');
JD = juliandate(utc_datetime);

% ISS pos vel vector 10 17 2024 12:00:00.00
r_ijk = [656.0758  6596.8028 1473.0872]' * 1000; % [m]
v_ijk = [-4.9631 -0.8127 5.7866]' * 1000; % [m/s]

[a,ecc,incl,RAAN,argp,nu,truelon,arglat,lonper] = ...
                                            ijk2keplerian(r_ijk, v_ijk);

num_orbits = 0.5;
P_coast = num_orbits * (2*pi)/(sqrt(mu/(a^3))); %Coast time
t = 1:1:P_coast;

options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
[Tout, Z] = ode45(@two_body_ode,t,[r_ijk v_ijk],options);

Initial_State = Z(1,:);

az = zeros(length(Z),1);
h = zeros(length(Z),1);
ra = zeros(length(Z),1);
dec = zeros(length(Z),1);

% Create measurements
for i = 1:length(Z)
    [az(i),h(i),ra(i),dec(i)] = getAzElRaDec(JD,pos_eci,Z(i,1:3));
end
variance = 2; % arcsec
z1 = ra + variance*randn();
z2 = dec + variance*randn();


% Gibbs' method IOD =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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
t1 = 0;
x1 = Z_IOD(end,:);

t0 = t1;
x0ref = x1';
x0 = zeros(6,1);

iter = 4;

% T_meas = [1, length(Tout)/2, Tout(end)];
T_meas = [1, 400, 800];
% T_meas = [1:length(Tout)];
% Tm is the time indexed at the 3 measurements
Tm = Tout(T_meas);

% Covariance matrix
for i = 1:length(Tm)
    Rm(:,:,i) = 2^2 * eye(2);
end

z1 = obs.ra;
z2 = obs.dec;

Zm = [z1(T_meas).*3600 ; z2(T_meas).*3600];

obs1time = utc_datetime;
obs2time = utc_datetime + seconds(length(Tout)/2);
obs3time = utc_datetime + seconds(length(Tout));

Rt = [lla2eci(lla, datevec(obs1time))' lla2eci(lla, datevec(obs2time))'...
    lla2eci(lla, datevec(obs3time))'];

% Rt = [];
% for i = 1:length(Tm)
%     obstime = utc_datetime+seconds(Tm(i));
%     Rt = [Rt,lla2eci(lla, datevec(obstime))'];
% end

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
        tk = T_meas(k);
        zk = Zm(:,k);
        Rk = Rm(:,:,k);
        rk = Rt(:,k)/1000;

        % determine LUMVE weighting matrix
        Ri = inv(Rk);

        % Propogate the reference state and STM but only if we're past t0
        if (tk>t0)
            options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
            [~, XX] = ode45(@(t,x) two_body_STM(t,x),[tkm1,tk],...
                [xref; Phi(:)],options);
            xref = XX(end,1:6)';
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


%% Orbit Plot  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
plotA = figure();
hold on
imData = imread('2_no_clouds_4k.jpg');  % Load the texture image

[xS, yS, zS] = sphere(50);  % Generate a sphere for the globe
earth_radius = 6378.1370 * 1000;  % Earth radius in meters

xSE = earth_radius * xS;
ySE = earth_radius * yS;
zSE = earth_radius * zS;

% Create the rotation matrices
Ry = R3(-lon+2);

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

plot3(Z(:,1), Z(:,2), Z(:,3), 'r', 'LineWidth', 4);
plot3(Z_IOD(:,1), Z_IOD(:,2), Z_IOD(:,3), 'b--', 'LineWidth', 4);

plotA.Position = [1600 100 1000 1000];

















