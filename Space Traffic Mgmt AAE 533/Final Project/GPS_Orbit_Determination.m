clear;clc;close all

mu = 3.986e14;
RE = 6.371E6; % Earth Radius [m]
c = 299792458; % [m/s]

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

plot(0,0,'g','Linewidth',5)
plot(0,0,'r','Linewidth',5)




% Get SWOT receiver observation data ======================================
% rinexFile = 'SWOT_L1_GPSP_RINEX_1280_20241031T155320_20241031T172120_PIC2_01.rnx';
% rinexFile = 'SWOT_L1_GPSP_RINEX_1280_20241031T002550_20241031T015920_PIC2_01.rnx';
rinexFile = 'SWOT_L1_GPSP_RINEX_1280_20241031T015930_20241031T032050_PIC2_01.rnx';
SWOT_Observation = rinexread(rinexFile);

% Time (GPS Time)
gps_time = SWOT_Observation.GPS.Time;

GPS_Leap_Seconds = 27;

utc_time = gps_time - seconds(GPS_Leap_Seconds);

% C1C: Represents the pseudorange measurement on L1 frequency using the 
% C/A (Coarse/Acquisition) code
L1_PseudoRange_CA = SWOT_Observation.GPS.C1C;

% S1C: Signal Strength
L1_CA_SSI = SWOT_Observation.GPS.S1C;

% C1W: Pseudorange measurement on the L1 frequency using the encrypted P(Y) code
L1_PseudoRange_PY = SWOT_Observation.GPS.C1W;

% S1W: Signal Strength
L1_PY_SSI = SWOT_Observation.GPS.S1W;

% C2W: Pseudorange measurement on the L2 frequency using the encrypted P(Y) code
L2_PseudoRange_PY = SWOT_Observation.GPS.C2W;

% S2W: Signal Strength
L2_PY_SSI = SWOT_Observation.GPS.S2W;


% L1C: Carrier phase measurement on L1 frequency using the C/A code
L1_CarrierPhase_CA = SWOT_Observation.GPS.L1C;

% L2W: Carrier phase measurement for the L2 frequency using the P(Y) code
L2_CarrierPhase_CA = SWOT_Observation.GPS.L2W;

% Get SWOT Precise Orbital Ephemeris ======================================
% Computed with through data collected from the satellite's onboard Doppler
% Orbitography and Radiopositioning Integrated by Satellite (DORIS) 
% instrument, which receives signals from globally distributed ground beacons
filePath = 'SWOT_VOR_AXVCNE20241121_171715_20241030_225923_20241101_005923.nc';

% info = ncinfo(filePath);  % Inspect the file metadata
% disp(info.Attributes);    % Display global attributes

% Read position data (ECEF position vectors)
position_POE = h5read(filePath, '/position');

% POSITION WAS FOUND TO BE ECEF!! WOW...
% Get metadata for a variable (e.g., 'x')
% xInfo = info.Variables(strcmp({info.Variables.Name}, 'position')); % Replace 'x' with the variable name
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
% POE_utcTime(1)
% POE_utcTime(end)

% Read orbit quality flags
orbit_qual = h5read(filePath, '/orbit_qual');


% Parse GPS ephemeris Oct 30 - Nov 1 ======================================
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
    
    t = 1:86400;
    
    options = odeset('RelTol', 1e-8,'AbsTol',1e-10);

    [r_ijk, v_ijk] = keplerian2ijk(a,e,i0,RAAN,omega,0);
    
    [T, Z] = ode45(@two_body_ode,t,[r_ijk, v_ijk],options);
    
    % plot3(Z(:,1), Z(:,2), Z(:,3), 'LineWidth', 2);
    % scatter3(Z(end,1), Z(end,2), Z(end,3), 300,'.');

    count = count + 1;

end

UserClockBias = 0;

% Estimate pseudo-range
rho_CA = L1_PseudoRange_CA + c*(SVClockBias - UserClockBias);

n = sqrt(mu / a^3);

% Compute satellite position for each time step
timeVec = 0:600:86400; % 1 day in 10-minute intervals
positions = [];
for t = timeVec
    % Mean anomaly
    M = mod(n * t, 2 * pi);

    % Solve Kepler's Equation for Eccentric Anomaly (E)
    E = M;
    for iter = 1:10
        E = M + e * sin(E);
    end

    % True anomaly
    v = atan2(sqrt(1 - e^2) * sin(E), cos(E) - e);

    % Argument of latitude
    u = v + omega;

    % Apply corrections
    dr = Crc + Crs; %TODO double check this
    di = Cic + Cic;
    c_radius = a * (1 - e * cos(E)) + dr; % Corrected radius
    i = i0 + di * cos(2 * u); % Corrected inclination

    % Compute satellite position in orbital plane
    x_orb = c_radius * cos(v);
    y_orb = c_radius * sin(v);

    % Convert to ECEF coordinates
    x = x_orb * cos(omega) - y_orb * sin(omega) * cos(i);
    y = x_orb * sin(omega) + y_orb * cos(omega) * cos(i);
    z = y_orb * sin(i);

    positions = [positions; x, y, z];

end

% SWOT_Observation
% GPS_ephemeris

% Compute satellite positions for 4 satellites
count = 1;
for i = [63,65,68,69,72]

    a = GPS_ephemeris.sqrtA(i)^2;
    e = GPS_ephemeris.Eccentricity(i);
    omega = GPS_ephemeris.omega(i) * 180/pi; % Radians->deg
    RAAN = GPS_ephemeris.OMEGA0(i) * 180/pi; % Radians->deg

    M0 = GPS_ephemeris.M0(i);

    omega = mod(omega, 360);
    if omega < 0
        omega = omega + 360;
    end

    RAAN = mod(RAAN, 360);
    if RAAN < 0
        RAAN = RAAN + 360;
    end

    mean_anomaly = M0*180/pi; % [deg]
    
    % Initial guess for E (E ~ M for small eccentricities)
    E = mean_anomaly;
    % Newton's method 
    for j = 1:20
        E_new = E - (E - e * sind(E) - mean_anomaly)...
                    / (1 - e * cosd(E));
        E = E_new;
    end

    % Compute true anomaly
    x = cosd(E) - e;
    y = sqrt(1-e^2) * sind(E);
    nu = (atan2d(y,x));
    if nu < 0
        nu = nu + 360;
    end
    
    i0 = GPS_ephemeris.i0(i) * 180/pi; % Radians->deg

    [satellite_positions(count,:) satellite_vel(count,:)] = keplerian2ijk(a,e,i0,RAAN,omega,nu);

    clk_bias_corr(count) = c*GPS_ephemeris.SVClockBias(i);

    count = count + 1;

end


options = odeset('RelTol', 1e-8,'AbsTol',1e-10);
max_sats = 33;
pseudoranges = [];
blockDoneFlag = 0;
svObsList = [];
ObsListIdx = [];
final_sv_list = [];
final_pseudoranges = [];
final_GPS_positions = [];
count = 1;

% Create pseudorange and GPS satellite position vectors
for i = 1:size(SWOT_Observation.GPS)-1

    % Check if adjacent times match and create list if so
    if (SWOT_Observation.GPS.Time(i) == SWOT_Observation.GPS.Time(i+1))
        % Append obs sv list
        svObsList = [svObsList SWOT_Observation.GPS.SatelliteID(i)];

        % Append obs idx list
        ObsListIdx = [ObsListIdx i];

        % Append pseudorange list
        pseudoranges = [pseudoranges SWOT_Observation.GPS.C1C(i)];

    else

        % Store final matching row and proceed
        svObsList = [svObsList SWOT_Observation.GPS.SatelliteID(i)];
        ObsListIdx = [ObsListIdx i];
        pseudoranges = [pseudoranges SWOT_Observation.GPS.C1C(i)];

        % Get index of closest ephemeris time
        [~, minIdx] = min(abs(SWOT_Observation.GPS.Time(i) - datetime(GPS_ephemeris.Time)));

        ephSVList = [];
        ephSVIdx = [];
        % Find and store all possible ephemerides for propogating
        for j = minIdx:minIdx+max_sats

            if abs( datetime(GPS_ephemeris.Time(j)) - datetime(GPS_ephemeris.Time(j+1))) < seconds(30*60)
                ephSVList = [ephSVList GPS_ephemeris.SatelliteID(j)];
                ephSVIdx = [ephSVIdx j];
            else
                break;
            end

        end

        % Store last SV
        ephSVList = [ephSVList GPS_ephemeris.SatelliteID(j)];
        ephSVIdx = [ephSVIdx j];

        % SV ID matching between obs and ephemeris list
        for l = 1:length(svObsList)

            for m = 1:length(ephSVList)

                if svObsList(l) == ephSVList(m)
                    % If there is a match, compute and store position data

                    % Compute and store 3D position from orbital elements
                    a = GPS_ephemeris.sqrtA( ephSVIdx(m) )^2;
                    e = GPS_ephemeris.Eccentricity( ephSVIdx(m) );
                    omega = GPS_ephemeris.omega( ephSVIdx(m) ) * 180/pi; % Radians->deg
                    RAAN = GPS_ephemeris.OMEGA0( ephSVIdx(m) ) * 180/pi; % Radians->deg

                    M0 = GPS_ephemeris.M0( ephSVIdx(m) );

                    omega = mod(omega, 360);
                    if omega < 0
                        omega = omega + 360;
                    end

                    RAAN = mod(RAAN, 360);
                    if RAAN < 0
                        RAAN = RAAN + 360;
                    end

                    mean_anomaly = M0*180/pi; % [deg]

                    % Initial guess for E (E ~ M for small eccentricities)
                    E = mean_anomaly;
                    % Newton's method 
                    for k = 1:20
                        E_new = E - (E - e * sind(E) - mean_anomaly)...
                                    / (1 - e * cosd(E));
                        E = E_new;
                    end

                    % Compute true anomaly
                    x = cosd(E) - e;
                    y = sqrt(1-e^2) * sind(E);
                    nu = (atan2d(y,x));
                    if nu < 0
                        nu = nu + 360;
                    end

                    i0 = GPS_ephemeris.i0( ephSVIdx(m) ) * 180/pi; % Radians->deg

                    % Compute 3D sat position at ephemeris time
                    [sat_pos, sat_vel] = keplerian2ijk(a,e,i0,RAAN,omega,nu);

                    % Compute time to propogate
                    tfinal = seconds( GPS_ephemeris.Time(ephSVIdx(m)) - SWOT_Observation.GPS.Time(ObsListIdx(l)) );

                    t = 1:abs(tfinal);

                    % Propagate sat position to pseudorange measurement time
                    if tfinal < 0 % Propogate backwards w/ negative velocity
                        [T, Z] = ode45(@two_body_ode,t,[sat_pos', -sat_vel'],options);
                    elseif tfinal > 0 % Propogate forward
                        [T, Z] = ode45(@two_body_ode,t,[sat_pos', sat_vel'],options);
                    elseif tfinal == 0 % Don't propagate
                        Z = sat_pos';
                    end

                    % Store propogated states and pseudorange
                    final_sv_list = [final_sv_list; svObsList(l)];
                    final_pseudoranges = [final_pseudoranges; pseudoranges(l)];
                    final_GPS_positions = [final_GPS_positions; Z(end,1:3)];

                end


            end


        end

        % final_pseudoranges
        [rcvr_pos(count,:), clk_bias(count)] = gps_receiver_position(final_pseudoranges, final_GPS_positions,clk_bias_corr);
        % disp(rcvr_pos(count,:))
        % disp(final_pseudoranges)
        % disp(final_GPS_positions)

        % Clear when done
        pseudoranges = [];
        blockDoneFlag = 0;
        svObsList = [];
        ObsListIdx = [];

        final_sv_list = [];
        final_pseudoranges = [];
        final_GPS_positions = [];

        count = count + 1;

    end

end

%


%

% test_pseudos = 1.0e+07*[2.5306;
%                         2.3995;
%                         2.0483;
%                         2.0181;
%                         2.2575];
% 
% test_sat_pos = 1.0e+07 *[1.6130    1.0999    1.8085;
%                         0.0984    1.5491    2.1821;
%                        -1.1124   -1.9387    1.4295;
%                         0.3916   -1.5111    2.1541;
%                         0.2535   -2.4878    0.9079];
% 
% [rcvr_pos_fake, clk_bias_fake] = gps_receiver_position(test_pseudos, test_sat_pos,clk_bias_corr);
% scatter3(rcvr_pos_fake(1,1),rcvr_pos_fake(1,2),rcvr_pos_fake(1,3),2000,'pc','filled','MarkerEdgeColor','k')
    

scatter3(rcvr_pos(1,1),rcvr_pos(1,2),rcvr_pos(1,3),200,'pg','filled','MarkerEdgeColor','k')
plot3(rcvr_pos(:,1), rcvr_pos(:,2), rcvr_pos(:,3),'g', 'LineWidth', 2);

pseudoranges = [SWOT_Observation.GPS.C1C(38);... %SVID8
                SWOT_Observation.GPS.C1C(39);... %SVID10
                SWOT_Observation.GPS.C1C(40);... %SVID13
                SWOT_Observation.GPS.C1C(41);... %SVID14
                SWOT_Observation.GPS.C1C(43)];   %SVID17

% pseudoranges = [SWOT_Observation.GPS.C1C(2);... %SVID8
%                 SWOT_Observation.GPS.C1C(3);... %SVID10
%                 SWOT_Observation.GPS.C1C(4);... %SVID13
%                 SWOT_Observation.GPS.C1C(5);... %SVID14
%                 SWOT_Observation.GPS.C1C(7)];   %SVID17

[rcvr_pos, clk_bias] = gps_receiver_position(pseudoranges, satellite_positions,clk_bias_corr);

% for i = 1:length(satellite_positions)
%     [sat_pos_matlab(i,:), sat_vel_matlab(i,:), ~] = eci2ecef([2024 10 31 2 00 00],satellite_positions(i,:),satellite_vel(i,:),[0 0 0]);
% end
% 
% [lla,~] = receiverposition(pseudoranges,sat_pos_matlab,zeros(5,1),sat_vel_matlab);
% matlab_recvr_position = lla2eci([lla(1) lla(2) lla(3)],[2024 10 31 2 00 00]);
% scatter3(matlab_recvr_position(1),matlab_recvr_position(2),matlab_recvr_position(3),800,'pm','filled','MarkerEdgeColor','k')

% Plot POE (Truth Orbit)

% Non-matlab ecef to eci implementation
% position_POE_ECI_test = zeros(3,length(position_POE));
% for i = 1:length(position_POE)
%     position_POE_ECI_test(:,i) = ecef_to_eci(position_POE(:,i), POE_utcTime(i));
% end
% plot3(position_POE_ECI_test(1,1085:1300), position_POE_ECI_test(2,1085:1300), position_POE_ECI_test(3,1085:1300),'cyan', 'LineWidth', 4);

% POE_utcTime_array = [year(POE_utcTime) month(POE_utcTime)...
%                     day(POE_utcTime) hour(POE_utcTime)...
%                     minute(POE_utcTime) second(POE_utcTime)];
% position_POE_ECI = zeros(3,length(position_POE));
% for i = 1:length(position_POE)
%     position_POE_ECI(:,i) = ecef2eci( POE_utcTime_array(i,:), position_POE(:,i));
% end
load('position_POE_ECI.mat');

scatter3(position_POE_ECI(1,1085),position_POE_ECI(2,1085),position_POE_ECI(3,1085),200,'pr','filled','MarkerEdgeColor','k')

% position_POE_ECI = ecef2eci( [2024 10 31 2 00 00], position_POE(:,1085));
% scatter3(position_POE_ECI(1),position_POE_ECI(2),position_POE_ECI(3),500,'pg','filled','MarkerEdgeColor','k')

% Truth orbit
% plot3(position_POE(1,1085:1300), position_POE(2,1085:1300), position_POE(3,1085:1300),'b', 'LineWidth', 2);
plot3(position_POE_ECI(1,1085:1300), position_POE_ECI(2,1085:1300), position_POE_ECI(3,1085:1300),'r', 'LineWidth', 2);

% Receiver Position
% scatter3(rcvr_pos(1),rcvr_pos(2),rcvr_pos(3),200,'py','filled','MarkerEdgeColor','k')

% GPS satellites
scatter3(satellite_positions(:,1),satellite_positions(:,2),satellite_positions(:,3),200,'.r')

plotA.Position = [0 100 1000 1000];

ax = gca;
ax.FontSize = 20;
legend('','Measured GPS Receiver Position','Truth Orbit')


%% SGP4

period = 48*60; % 48 hours * 60 minutes = 2 periods at GEO

t = 0:60:period*60; % increments of 60 seconds

cc=0; % set line counter

    fid = fopen('GEO_SAT.txt'); % load the TLE
    
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
        
            % how far shall the TLE be propagated [minutes]
            % tsince = period; 


            % extract position and velocity
            for i = 1:period
                tsince = i;

                [satrec, R_SGP4(i,:), V_SGP4(i,:)] = sgp4(satrec, tsince); 
            end
        

            % Numeric integrator ------------------------------------------

            orb_elements = str2double(strsplit(tline2, ' '));
            incl = orb_elements(3); % [deg]
            RAAN = orb_elements(4); % [deg]
            ecc = orb_elements(5)*10^-7; 
            arg_per = orb_elements(6); % [deg]
            mean_anomaly = orb_elements(7); % [deg]
            mean_motion  = orb_elements(8); % [rev/day]
            
            % Initial guess for E (E ~ M for small eccentricities)
            E = mean_anomaly;
            % Newton's method 
            for j = 1:20
                E_new = E - (E - ecc * sind(E) - mean_anomaly)...
                            / (1 - ecc * cosd(E));
                E = E_new;
            end
   
        
            % Compute true anomaly
            x = cosd(E) - ecc;
            y = sqrt(1-ecc^2) * sind(E);
            nu = (atan2d(y,x));
            if nu < 0
                nu = nu + 360;
            end
    
            % Compute semi-major axis 
            a = (mu / (4*pi^2*(mean_motion/86400)^2) )^(1/3);
        
            alt = (norm([r_ijk(1), r_ijk(2), r_ijk(3)]) - RE)/1000;
   
 
        end

    end


