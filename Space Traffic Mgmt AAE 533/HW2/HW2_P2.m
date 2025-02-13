clear;clc;close all

load('data.mat');

% Compute polar motion matrix =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% https://datacenter.iers.org/data/207/bulletinb-353.txt

xp = 125.744 / 1000 / 3600 * pi/180; % [mas] milliarc seconds / 1000 / 3600 * 180/pi = [rad]
yp = 455.691 / 1000 / 3600 * pi/180; % [mas] milliarc seconds / 1000 / 3600 * 180/pi = [rad]

polarRot = [1   0    xp;
            0   1   -yp;
           -xp  yp   1 ];

station_ITRF = [data.StationXITRFkm data.StationYITRFkm data.StationZITRFkm]';

% Compute sidereal time/angle  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
JD = data.TimeMJD + 2400000.5; 

% Compute Julian Centuries, T1
T1 = (JD - 2451545) ./ 36525; 

% Compute JD0, keep only integer portion
JD0 = floor(JD); 

% Compute T0
T0 = (JD0 - 2451545) ./ 36525; 

% Mean sidereal time of Greenwich 0h Earth time (UT)
UT = 24110.54841 + 8640184.812866*T0 + 0.093104*T1.^2 - 0.0000062*T1.^3 + (1.0027279093*(mod(JD,1)-.5).*24).*3600;
UT = UT ./ 3600;

% Compute sidereal time
sidereal_time = mod( UT , 24);

% Compute sidereal angle [deg]
sidereal_angle = sidereal_time * 15

% Compute Gregorian time
gregorian_time = mjd2gregorian(data.TimeMJD);

station_GTOD = zeros(3,length(station_ITRF)); % Preallocate size
station_TOD = zeros(3,length(station_ITRF)); % Preallocate size
station_MOD = zeros(3,length(station_ITRF)); % Preallocate size
station_J2000 = zeros(3,length(station_ITRF)); % Preallocate size

for i = 1:length(sidereal_angle)

    % Rotate  ITRF by Polar Motion to obtain GTOD/PEF
    station_GTOD(:,i) = polarRot' * station_ITRF(:,i);

    % Rotate GTOD by Sidereal Time to obtain TOD
    R3_st_temp = R3(sidereal_angle(i) * pi/180);
    station_TOD(:,i) = R3_st_temp' * station_GTOD(:,i);

        % Rotate TOD by Nutation to obtain MOD
    N = nutation(gregorian_time(i,1:end-1));
    station_MOD(:,i) = N' * station_TOD(:,i);

    % Rotate MOD by Precession to obtain J2000
    P = precession(JD(i));
    station_J2000(:,i) = P' * station_MOD(:,i);
end

station_J2000'

validate = station_J2000' - [data.StationXJ2000km data.StationYJ2000km data.StationZJ2000km];


figure()
hold on
grid minor
plot(validate)
ylabel("Error [km]",'FontSize',15)
xlabel("Sample Number",'FontSize',15)
legend("Station X J2000", "Station Y J2000", "Station Z J2000",'FontSize',15)
title('Computed J2000 Coordinates vs. Validation Data','FontSize',15)

