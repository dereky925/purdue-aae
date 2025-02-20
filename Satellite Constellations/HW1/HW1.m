% =========================================================================
% 
% Filename:       HW1_P1.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE590 - Satellite Constellations and Formation
% Professor:      Dr. David Arnas
% Assignment:     HW 1
% Semester:       Spring 2025
% 
% Description: Homework 1
%
% =========================================================================

%% Problem 1
% Verification of own keplerian <-> cartesian functions
clear;clc

x = [-6659000 -1896000 180700 1400 -4300 6100];

[a,ecc,incl,raan,argp,nu] = eci2keplerian(x(1:3)./1000,x(4:6)./1000);

% Verification
[a_v,ecc_v,incl_v,RAAN_v,argp_v,nu_v,...
    truelon,arglat,lonper] = ijk2keplerian(x(1:3), x(4:6));


[r_vec,rdot_vec] = keplerian2eci(a,ecc,incl,raan,argp,nu);
[r_vec_v,rdot_vec_v] = keplerian2ijk(a*1000,ecc,incl,raan,argp,nu);


M_ISS = deg2rad(318.596);
ecc_ISS = 0.1;

[nu,~] = m2nu(M_ISS,ecc_ISS)

M = nu2m(nu, ecc);
M = rad2deg(M)

%% Problem 2
% Testing and plotting of an orbit

clear;clc;close all

mu = 398600.4418; % [km^3/s^2]

% Given
w_earth = 7.2921151467E-5; % rad/s

ecc = 0.74;
incl = 63.4349; % deg
raan = -86.915798; % deg
argp = 270; % deg
nu = 0;

% Compute T (period) for 2 orbits per sidereal day
T = 180/EARTH_RATE * 60 * 60; % [s]

% Compute semi-major axis
a = ( mu*(T/(2*pi))^2 )^(1/3);

% a = 3.4806E4;
% T = 2*pi*sqrt(a^3/mu);
% ecc = 0.6;
% incl = 180-75; % deg
% nu = 0;
% argp = 160; % deg
% raan = 0; % deg

% Convert orbital elements to cartesian
[r_vec,rdot_vec] = keplerian2eci(a,ecc,incl,raan,argp,nu);
% [r_vec,rdot_vec] = keplerian2ijk(a,ecc,incl,raan+360,argp,nu);

% Create time vector starting from 0
t = linspace(0,T,100);

% Numerical propagation
options = odeset('RelTol', 1e-13,'AbsTol',1e-15);
[~,eci_coords] = ode45(@two_body_ode,t,[r_vec rdot_vec],options);

% Analytical propagation
[eci_pos, eci_vel] = AnalyticalOrbitPropagator(a,ecc,incl,raan,argp,nu,t);

% =========================================================================
% Orbit Plot
orbit_plot_figure = OrbitPlotSetup();

% Plot Orbit
plot3(eci_coords(:,1), eci_coords(:,2), eci_coords(:,3), 'r', 'LineWidth', 3);
plot3(eci_pos(:,1), eci_pos(:,2), eci_pos(:,3), 'b--', 'LineWidth', 3);
legend('','Numerical ODE45 Orbit','Analytical Orbit')
orbit_plot_figure.Position = [0 0 1000 1000];

% =========================================================================
% Compute radial distance
% vecnorm(data, 2 = norm 2, 2 = norm by rows)
% Numerical Soln
radial_distance = vecnorm(eci_coords(:,1:3),2,2);
vel_mag = vecnorm(eci_coords(:,4:6),2,2);

% Analytical Soln
radial_distance_analytical = vecnorm(eci_pos,2,2);
vel_mag_analytical = vecnorm(eci_vel,2,2);

rad_vel_plot = figure();
rad_vel_plot.Position = [100 0 1600 1000];

t_hours = t/3600;
subplot(2,2,1)
hold on
grid minor
plot(t_hours,radial_distance,'r','LineWidth',4)
plot(t_hours,radial_distance_analytical,'b--','LineWidth',4)
ax = gca;
ax.FontSize = 20;
title('Radial Distance vs. Time')
ylabel('Radial Distance [km]')
xlabel('Time [Hours]')
legend('Numerical ODE45','Analytical','Location','best')

subplot(2,2,2)
hold on
grid minor
plot(t_hours,vel_mag,'r','LineWidth',4)
plot(t_hours(1:end-1),vel_mag_analytical(1:end-1),'b--','LineWidth',4)
ax = gca;
ax.FontSize = 20;
title('Velocity Magnitude vs. Time')
ylabel('Velocity Magnitude [km/s]')
xlabel('Time [Hours]')
legend('Numerical ODE45','Analytical','Location','best')

subplot(2,2,3)
hold on
grid minor
plot(t_hours,radial_distance-radial_distance_analytical,'r','LineWidth',4)
ax = gca;
ax.FontSize = 20;
title('Radial Distance Difference vs. Time')
xlabel('Time [Hours]')
ylabel('Radial Distance [km]')


subplot(2,2,4)
hold on
grid minor
plot(t_hours(1:end-1),(vel_mag(1:end-1)-vel_mag_analytical(1:end-1)),'r','LineWidth',4)
ax = gca;
ax.FontSize = 20;
title('Velocity Magnitude Difference vs. Time')
xlabel('Time [Hours]')
ylabel('Velocity Magnitude [km/s]')


% =========================================================================
% Plot orbital elements over a period

% Compute orbital elements from cartesian
for i = 1:length(eci_coords)
    % Numeric
    [a(i),ecc(i),incl(i),raan(i),argp(i),nu(i)] = ...
    eci2keplerian(eci_coords(i,1:3), eci_coords(i,4:6));
end

for i = 1:length(eci_pos)
    % Analytical
    [a_a(i),ecc_a(i),incl_a(i),raan_a(i),argp_a(i),nu_a(i)] = ...
    eci2keplerian(eci_pos(i,:), eci_vel(i,:));
end

% Create cell arrays to index entire variables
orb_data_num = {a, ecc, incl, raan, argp, nu};
orb_data_an = {a_a, ecc_a, incl_a, raan_a, argp_a, nu_a};
titles = {'a [km]', 'e', 'i [°]', 'RAAN [°]', '\omega [°]', '\nu [°]'};
ylabels = {'', '', '', '', '', ''}; 

orbital_elements = figure(); 
sgtitle('Numerical ODE45: Orbital Elements over 1 Period','FontSize',30)
orbital_elements.Position = [200 50 1800 1000];

for i = 1:6
    subplot(2,3,i);       % pick the i-th subplot
    hold on; grid minor;  % set up axes
    plot(t_hours, orb_data_num{i},'r', 'LineWidth', 4);
    
    ax = gca;
    ax.FontSize = 15;

    title(titles{i}, 'FontSize', 20,'FontWeight', 'bold');
    ylabel(ylabels{i}, 'FontSize', 15);
    xlabel('Time [Hours]')
    
end

orbital_elements_an = figure(); 
sgtitle('Analytical Method: Orbital Elements over 1 Period','FontSize',30)
orbital_elements_an.Position = [300 50 1800 1000];

for i = 1:6
    subplot(2,3,i);       % pick the i-th subplot
    hold on; grid minor;  % set up axes
    plot(t_hours(1:end-1), real(orb_data_an{i}(1:end-1)),'b', 'LineWidth', 4);
    
    ax = gca;
    ax.FontSize = 15;

    title(titles{i}, 'FontSize', 20,'FontWeight', 'bold');
    ylabel(ylabels{i}, 'FontSize', 15);
    xlabel('Time [Hours]')
end

% =========================================================================
% Problem 3
% Ground track plots
% =========================================================================
arms_lat = 40.43094;
arms_lon = -86.915798;

% Need data in ECEF to compute lat/long
% 1. Convert UTC to JD Time
% 2. Compute Greenwich Mean Sidereal Time using JD
% 3. Rotate ECI coordinates to ECEF with GMST

% Create time vector starting from 0
t = linspace(0,2*T,800);

% Numerical propagation
options = odeset('RelTol', 1e-12,'AbsTol',1e-15);
[time_vec,eci_coords] = ode45(@two_body_ode,t,[r_vec rdot_vec],options);

% Offset time such that apogee is above Purdue University
time_vec = time_vec + 28*60*60;

for i=1:length(eci_coords)

    sidereal_angle = utc2siderealangle(2025,1,25,11,32,time_vec(i));

    JD = utc2jd(2025,1,25,11,32,time_vec(i));

    utc_datetime = datetime(2025,1,25,11,32,time_vec(i), 'TimeZone', 'UTC');
    % JD = juliandate(utc_datetime);

    pos_ecef(i,:) = (R3(deg2rad(sidereal_angle)) * eci_coords(i,1:3)')';

    % This function is buggy, need to fix
    % station_eci = llasa2eci(arms_lat, arms_lon, 0, sidereal_angle);
    station_eci = lla2eci([40.43094 -86.915798 0], datevec(utc_datetime));

    % visible_indices = check_station_visibility(eci_coords, station_eci/1000);

    % Matlab tool verification
    % utcTime = datetime(2025,1,25,11,32,time_vec(i));
    % pos_ecef(i,:) = eci2ecef(utcTime,eci_coords(i,1:3),eci_coords(i,4:6),[0 0 0]);

    % Compute Latitude
    latitude(i) = asind( pos_ecef(i,3)/norm(pos_ecef(i,:)) );
    % Compute Longitude
    longitude(i) = atan2d(pos_ecef(i,2),pos_ecef(i,1));

    % Try computing geodetic lat/long
    [lat(i), lon(i), h(i)] = ecef2geodetic(pos_ecef(i,1)*1000,...
        pos_ecef(i,2)*1000, pos_ecef(i,3)*1000);

    hours_timestamp = (time_vec - time_vec(1)) / 3600;

    [az(i),h(i),~,~,~,~] = getAzElRaDec(JD,station_eci,eci_coords(i,1:3)*1000);
    azel_idx(i) = i;

    if h(i) < 15
        h(i) = NaN;
        az(i) = NaN;
        azel_idx(i) = NaN;
    end

end

visible_mask = ~isnan(azel_idx);  % True where satellite is visible
start_idx = find(diff([0, visible_mask]) == 1); % Start of visibility
end_idx = find(diff([visible_mask, 0]) == -1);  % End of visibility

total_visibility_hours = 0;
for i = 1:length(start_idx)
    t_start = hours_timestamp(start_idx(i));
    t_end = hours_timestamp(end_idx(i));

    % Add the duration of this visibility period
    total_visibility_hours = total_visibility_hours + (t_end - t_start);
end
total_visibility_hours


% =========================================================================
% Ground Track Plot 
[ground_track_fig] = GroundTrackPlotSetup();
ground_track_fig.Position = [300 0 1200 600];
plot(lon,lat,'.r','MarkerSize',15);
plot(lon(visible_mask),lat(visible_mask),'.c','MarkerSize',15);
% plot(longitude,latitude,'.b','MarkerSize',15);
scatter(arms_lon,arms_lat,500,'p','MarkerEdgeColor','k','MarkerFaceColor','y')

legend('Ground Track','Satellite is visible','Neil Armstrong Hall')
% =========================================================================


% =========================================================================
% AzEl Plot
az_el_figure = figure();

subplot(2,1,1)
hold on
grid minor
plot(hours_timestamp,az,'Linewidth',4)
title('Azimuth vs. Time')
xlabel('Time [hours]')
ylabel('Azimuth [°]')
ax = gca;
ax.FontSize = 20;

subplot(2,1,2)
hold on
grid minor
plot(hours_timestamp,h,'Linewidth',4)
title('Elevation vs. Time')
xlabel('Time [hours]')
ylabel('Elevation [°]')
ax = gca;
ax.FontSize = 20;
az_el_figure.Position = [500 0 1000 800];








