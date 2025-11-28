% =========================================================================
% 
% Filename:       HW3.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE590 - Satellite Constellations and Formation
% Professor:      Dr. David Arnas
% Assignment:     HW 3
% Semester:       Spring 2025
% 
% Description: Homework 3
%
% =========================================================================

%% Problem 1 - 1 Orbit
clear;clc;close all

a = 7000.0015; % [km]
ex = 4.6628E-4;
ey = 0;
i = 97.8739; % [deg]
RAAN = 0;
argLat_Sat1 = 0; % [deg]
argLat_Sat2 = 90; % [deg]

T = 1.61784; % [hours]
t = 1:T*3600;

[initPos_Sat1,initVel_Sat1] = keplerian2eci(a,ex,ey,i,RAAN,argLat_Sat1,'alt');
[initPos_Sat2,initVel_Sat2] = keplerian2eci(a,ex,ey,i,RAAN,argLat_Sat2,'alt');

options = odeset('RelTol',1E-12,'AbsTol',1E-12);

[t_out_1,eci_out_1] = ode45(@(t,x) two_body_J2_ECI(t,x),...
                t, [initPos_Sat1 initVel_Sat1], options);

[t_out_2,eci_out_2] = ode45(@(t,x) two_body_J2_ECI(t,x),...
                t, [initPos_Sat2 initVel_Sat2], options);

orbFig = OrbitPlotSetup;
orbFig.Position = [0 0 1000 1000];
plot3(eci_out_1(:,1),eci_out_1(:,2),eci_out_1(:,3),'b','LineWidth',4);
plot3(eci_out_1(:,1),eci_out_1(:,2),eci_out_1(:,3),'r','LineWidth',4);

% Compute radial distance for plotting vs. latitude
radial_distance_sat1 = vecnorm(eci_out_1(:,1:3),2,2);
radial_distance_sat2 = vecnorm(eci_out_2(:,1:3),2,2);

% Compute geocentric latitude
for i = 1:length(eci_out_1)
    lat_gc_sat1(i) = atan2d(eci_out_1((i),3), sqrt(eci_out_1((i),1)^2 + eci_out_1((i),2)^2));
end

for i = 1:length(eci_out_2)
    lat_gc_sat2(i) = atan2d(eci_out_2((i),3), sqrt(eci_out_2((i),1)^2 + eci_out_2((i),2)^2));
end

dist_vs_lat_fig = figure('Color','white','Position',[500 0 1000 1000]);
subplot(2,1,1)
hold on; grid minor
plot(lat_gc_sat1,radial_distance_sat1,'r','LineWidth',2);
xlabel('Geocentric Latitude [°]')
ylabel('Radial Distance [km]')
title('Radial Distance vs. Latitude for 1 Orbital Period: \theta = 0°')
ax = gca;
ax.FontSize = 20;

subplot(2,1,2)
hold on; grid minor
plot(lat_gc_sat2,radial_distance_sat2,'b','LineWidth',2);
xlabel('Geocentric Latitude [°]')
ylabel('Radial Distance [km]')
title('Radial Distance vs. Latitude for 1 Orbital Period: \theta = 90°')

ax = gca;
ax.FontSize = 20;

%% Problem 1 - 100 Orbits
clear;clc;close all

a = 7000.0015; % [km]
ex = 4.6628E-4;
ey = 0;
i = 97.8739; % [deg]
RAAN = 0;
argLat_Sat1 = 0; % [deg]
argLat_Sat2 = 90; % [deg]

T = 1.61784; % [hours]
t = 1:T*3600*100;

[initPos_Sat1,initVel_Sat1] = keplerian2eci(a,ex,ey,i,RAAN,argLat_Sat1,'alt');
[initPos_Sat2,initVel_Sat2] = keplerian2eci(a,ex,ey,i,RAAN,argLat_Sat2,'alt');

options = odeset('RelTol',1E-12,'AbsTol',1E-12);

[t_out_1,eci_out_1] = ode45(@(t,x) two_body_J2_ECI(t,x),...
                t, [initPos_Sat1 initVel_Sat1], options);

[t_out_2,eci_out_2] = ode45(@(t,x) two_body_J2_ECI(t,x),...
                t, [initPos_Sat2 initVel_Sat2], options);

orbFig = OrbitPlotSetup;
orbFig.Position = [0 0 1000 1000];
plot3(eci_out_1(:,1),eci_out_1(:,2),eci_out_1(:,3),'b','LineWidth',4);
plot3(eci_out_1(:,1),eci_out_1(:,2),eci_out_1(:,3),'r','LineWidth',4);

% Compute radial distance for plotting vs. latitude
radial_distance_sat1 = vecnorm(eci_out_1(:,1:3),2,2);
radial_distance_sat2 = vecnorm(eci_out_2(:,1:3),2,2);

% Compute geocentric latitude
for i = 1:length(eci_out_1)
    lat_gc_sat1(i) = atan2d(eci_out_1((i),3), sqrt(eci_out_1((i),1)^2 + eci_out_1((i),2)^2));
end

for i = 1:length(eci_out_2)
    lat_gc_sat2(i) = atan2d(eci_out_2((i),3), sqrt(eci_out_2((i),1)^2 + eci_out_2((i),2)^2));
end

dist_vs_lat_fig = figure('Color','white','Position',[500 0 1000 1000]);
subplot(2,1,1)
hold on; grid minor
plot(lat_gc_sat1,radial_distance_sat1,'r','LineWidth',2);
xlabel('Geocentric Latitude [°]')
ylabel('Radial Distance [km]')
title('Radial Distance vs. Latitude for 100 Orbital Periods: \theta = 0°')
ax = gca;
ax.FontSize = 20;

subplot(2,1,2)
hold on; grid minor
plot(lat_gc_sat2,radial_distance_sat2,'b','LineWidth',2);
xlabel('Geocentric Latitude [°]')
ylabel('Radial Distance [km]')
title('Radial Distance vs. Latitude for 100 Orbital Periods: \theta = 90°')

ax = gca;
ax.FontSize = 20;


%% Problem 2
clear;clc;close all

a = 2.65618E4; % [km]
ecc = 0.74;
incl = 63.4271; % [°]
argp = 270; % [°]
RAAN = 0; % [°]
nu = 180; % [°]

num_orbits = 1;
T = num_orbits * (2*pi)*(sqrt(a^3/MU)); % Orbit period
t = 1:T;

% Create initial states for cartesian propagation
[r_ijk,v_ijk] = keplerian2eci(a,ecc,incl,RAAN,argp,nu);

% Form vector for orbital element propagation
element_vec = [a,ecc,incl,RAAN,argp,nu];

options = odeset('RelTol',1E-12,'AbsTol',1E-12);

% Cartesian propagation
[~, eci_vec] = ode45(@(t,x) two_body_J2_ECI(t,x),...
                t, [r_ijk v_ijk], options);

% Propagate orbital elements
[~, orb_el_out] = ode45(@(t,x) orbital_elements_eom_J2_gauss(t,x),...
                t, element_vec, options);


% Plot orbital elements over a period

% Compute orbital elements from cartesian
for i = 1:length(eci_vec)
    [a(i),ecc(i),incl(i),raan(i),argp(i),nu(i)] = ...
    eci2keplerian(eci_vec(i,1:3), eci_vec(i,4:6));

    raan(i) = mod(raan(i) + 180, 360) - 180;

    M(i) = rad2deg(  nu2m( deg2rad(nu(i)), ecc(i))  );
end

for i = 1:length(orb_el_out)
    orb_el_out(i,6) = rad2deg(nu2m(deg2rad(orb_el_out(i,6)),orb_el_out(i,2)));
end

t_hours = t/3600;

% Create cell arrays to index entire variables
M = mod(M + 180, 360) - 180;
cartesian = {a, ecc, incl, raan, argp, M};


orb_el = {orb_el_out(:,1), orb_el_out(:,2), orb_el_out(:,3),...
         mod(orb_el_out(:,4) + 180, 360) - 180, orb_el_out(:,5), mod(orb_el_out(:,6) + 180, 360) - 180};

for i = 1:6
    cartesian_mat(:,i) = cartesian{i};
    orb_el_mat(:,i) = orb_el{i};
end

% Compute the difference
diff_mat = cartesian_mat - orb_el_mat;


titles = {'a [km]', 'e', 'i [°]', 'RAAN [°]', '\omega [°]', 'M [°]'};
ylabels = {'', '', '', '', '', ''}; 

orbital_elements = figure('color','white'); 
sgtitle('Orbital Elements over 1 Period','FontSize',30)
orbital_elements.Position = [1000 50 1800 1000];

% For legend only
plot(NaN,NaN,'r','LineWidth',2)
plot(NaN,NaN,'b--','LineWidth',1)

for i = 1:6
    subplot(3,2,i);       % pick the i-th subplot
    hold on; grid minor;  % set up axes
    plot(t_hours, cartesian{i},'r', 'LineWidth', 4);
    plot(t_hours, orb_el{i},'b--', 'LineWidth', 4);
    % plot(t_hours(3:end), diff_mat(3:end,i),'r', 'LineWidth', 4);

    ax = gca;
    ax.FontSize = 15;

    title(titles{i}, 'FontSize', 20,'FontWeight', 'bold');
    ylabel(ylabels{i}, 'FontSize', 15);
    xlabel('Time [Hours]')
    
end

% legend('Cartesian Propagation','Orbital Element Propagation');


%% Problem 3
clear;clc;close all

% Given
Cd = 2.2; % Rule of thumb for a lot of space missions
Sm = 0.01; % m^2/kg

% Assume stationary Earth atmosphere
v_rel_hat = 1; 

% Evaluate altitudes from 400 to 800 km
altArray = 400:10:800;

for i = 1:length(altArray)

    alt = altArray(i);

    % Compute semi-major axis [km], circular orbit
    a = EARTH_RADIUS + alt;   

    % Compute velocity for circular orbit
    v_km_s = sqrt(MU / a);   % [km/s]

    % Table lookup for air density
    rho = getAirDensity(alt); % kg/m^3

    % Convert density from kg/m^3 to kg/km^3 -- (1 m^3 = 1e-9 km^3)
    rho_kg_km3 = rho * 1e9;  

    % Convert S/m to (km^2)/kg for consistent units with v in km/s
    S_over_m_km2 = Sm / (1e6);  % 1 m^2 = 1e-6 km^2

    % Drag acceleration gamma_f [km/s^2]
    gamma_f = 0.5 * rho_kg_km3 * Cd * S_over_m_km2 * (v_km_s^2);

    % da/dt
    da_dt_km_s = -2 * sqrt(a^3 / MU) * gamma_f; % [km/s]

    % Convert from km/s to m/day
    km_per_day = da_dt_km_s * 86400;  % multiply by seconds in a day
    m_per_day  = km_per_day * 1000;

    % Store the magnitude of the decay
    decay_m_per_day(i) = -m_per_day;

end


figure('Color','w','Position',[0 0 1400 1000]);
hold on;grid minor;
plot(altArray, decay_m_per_day, 'b-o','LineWidth',2);
xlabel('Altitude (km)');
ylabel('Altitude Decay (meters/day)');
title('Daily Altitude Decay vs. Altitude (400-800 km)');
ax = gca;
ax.FontSize = 20;


%% Problem 4

clear; clc; close all;

Cd = 2.2;       % Drag coefficient
Sm = 0.01;      % [m^2/kg], ballistic parameter

% Initial altitudes from 100 up to 700 km
altList = 100:10:700;

% For each initial altitude, simulate the orbit decay down to 100 km
for i = 1:length(altList)
    alt0 = altList(i);
    timeTo100_yrs(i) = computeDecayTime(alt0, Cd, Sm);
end

figure('Color','w','Position',[0 0 1400 1000]);
plot(altList, timeTo100_yrs, 'b','LineWidth',2);
xlabel('Initial Altitude (km)');
ylabel('Time to Reach 100 km (years)');
title('Time to Reenter vs. Initial Altitude');
grid minor;
set(gca,'FontSize',14);
ax = gca;
ax.FontSize = 20;


function Tyrs = computeDecayTime(alt0, Cd, Sm)

    alt = alt0;   % current altitude [km]
    a   = EARTH_RADIUS + alt;   % semi-major axis

    dt_s = 600;   % time step [s]
    time_s = 0;   % accumulate total time [s]

    while alt > 100

        a_km = EARTH_RADIUS + alt;  % semi-major axis
        v_km_s = sqrt(MU / a_km);

        rho = getAirDensity(alt);  % kg/m^3

        % Convert to kg/km^3
        rho_kg_km3 = rho * 1e9;  

        % Convert S/m from m^2/kg to km^2/kg
        Sm_km2 = Sm / 1e6;

        gamma_f = 0.5 * rho_kg_km3 * Cd * Sm_km2 * (v_km_s^2);

        % da/dt in km/s
        da_dt_km_s = -2 * sqrt(a_km^3 / MU) * gamma_f;

        dalt_dt_km_s = da_dt_km_s;

        % Update altitude
        alt_new = alt + dalt_dt_km_s * (dt_s);

        if alt_new < 100
            alt_new = 100;
        end

        % Accumulate time
        time_s = time_s + dt_s;

        % Overwrite altitude
        alt = alt_new;
    end

    % Convert total time to years
    Tyrs = time_s / (86400*365);

end

%% Problem 5

% Landsat 7 
% https://celestrak.org/satcat/table-satcat.php?NAME=Landsat%207&PAYLOAD=1&MAX=500

% LANDSAT 7               
% 1 25682U 99020A   25057.58092020  .00001356  00000+0  28287-3 0  9992
% 2 25682  97.8708  78.2056 0001107  89.2635   9.7584 14.61396596376283

i = 97.8708;
raan = 78.2056;
e = 0.0001107;
omega = 89.2635; 
M = 9.7584;     % mean anomaly [deg]
n = 14.61396596; % mean motion [revs/day]

% T = 2*pi*sqrt(a^3/MU);
T = 1/(n/86400);

a = (MU*(T/(2*pi))^2)^(1/3)

nu = rad2deg(m2nu(deg2rad(M),e))


