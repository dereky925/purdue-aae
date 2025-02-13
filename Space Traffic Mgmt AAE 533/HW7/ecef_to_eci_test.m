clear;clc;

utc_datetime = datetime(2024, 10, 16, 12, 0, 0, 'TimeZone', 'UTC');
JD = juliandate(utc_datetime);

% Define Earth's rotation rate in rad/s
omega = 7.2921159e-5; % Earth's angular velocity in rad/s, sidereal day
% omega = 15 * pi/180 / 60 / 60 

% Define the elapsed time since the reference epoch in seconds
t = JD * 24 * 60 * 60; % define time in seconds since reference epoch;

% Calculate the Earth rotation angle
theta = omega * t;

% Example ECI position vector
r_eci = [7000; 0; 0]; % Example vector in ECI frame

% Rotate to get the ECEF position
r_ecef = R3(theta) * r_eci
r_ecef = R3(theta)' * r_ecef

[r_ecef_matlab,~,~] = eci2ecef(datevec(utc_datetime) ,r_eci',[0 0 0],[0 0 0]);


