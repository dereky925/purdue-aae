clear;clc;close all

c = 3E8; % Speed of light [m/s]

h_ = [90 75 25]; % Observed elevation [°]
z_ = 90 - h_;
p = 1013.25; % Atmospheric pressure, assume sea-level [millibar]
T = 293.15; % Local temperature, [°K]

RH = 0.72; % Average relative humidity in Indiana
A = 1.2378847E-5;
B = -1.9121316E-5;
C = 33.93711047;
D = -6.3431645E3;
e = 0.01*RH*exp(A*T^2 + B*T + C + D/T); % [millibar]

dz = 16.271.*tand(z_).*(1+0.0000394.*tand(z_).^2.*((p-0.156*e)/T)).*((p-0.156*e)/T)-0.0749*(tand(z_).^3+tand(z_))*(p/1000);
dz = dz / 3600 % Convert arcseconds to degrees

h = h_ - dz; % True elevation [degrees]

geo_alt = 36000; % GEO Altitude [km]
obs_horz = geo_alt * cosd(h_); % Observed horizontal distance [km]
obs_vert = geo_alt * sin(h_); % Observed vertical distance [km]

true_horz = geo_alt * cosd(h); % True horizontal distance [km]
true_vert = geo_alt * sin(h); % True vertical distance [km]

error_horz = abs(obs_horz - true_horz) % [km]
error_vert = abs(obs_vert - true_vert) % [km]

% 3.b.
signal_travel_time = geo_alt * 1000 / c; % GEO altitude / speed of light = [s]

vel_sat = 3000; % Typical velocity of GEO satellite [m/s]

dist_error = vel_sat * signal_travel_time / 1000 % [km]

% 3.c.
w_Earth = 15 / 3600 * pi/180; % Rotation speed of Earth [rad/s]
angular_error = signal_travel_time * w_Earth; % [rad]
absolute_error = geo_alt * angular_error   % [km]

total_error = norm(error_horz(3),error_vert(3)) + dist_error + absolute_error % [km]
