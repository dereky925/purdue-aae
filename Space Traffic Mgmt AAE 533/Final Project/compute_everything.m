function [rcvr_pos, rcvr_vel, clock_bias, satellite_positions] = compute_everything(obsFile, gpsfile, obsVec, gpsVec)

    c = 299792458;
    pseudoranges = [];
    SWOT_Observation = rinexread(obsFile);
    ps = SWOT_Observation.GPS;
    
    for i = obsVec
        pseudoranges = [pseudoranges ps.C1C(i)];
    end
    
    % Compute satellite positions for 4 satellites
    combinedTimetable = vertcat(gpsfile);

    % Time might be unique'd automatically, so only check for SatelliteID after
    columnsToCheck = combinedTimetable(:, {'SatelliteID'});
    
    % Find unique rows based on the specified columns
    [~, uniqueIdx] = unique(columnsToCheck, 'rows', 'stable');
    
    % Filter the original timetable
    uniqueTimetable = combinedTimetable(uniqueIdx, :);

    GPS_ephemeris = sortrows(uniqueTimetable, {'Time', 'SatelliteID'});

    count = 1;
    for i = gpsVec
    
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

    utcTimeRef = ps.Time(obsVec(1));
    utcTime_array = [year(utcTimeRef) month(utcTimeRef)...
                    day(utcTimeRef) hour(utcTimeRef)...
                    minute(utcTimeRef) second(utcTimeRef)];

    for i = 1:length(satellite_positions)
        [sat_pos_matlab(i,:), sat_vel_matlab(i,:), ~] = eci2ecef(utcTime_array,satellite_positions(i,:),satellite_vel(i,:),[0 0 0]);
    end
    
    % pseudoranges
    [lla,rcvr_vel] = receiverposition(pseudoranges',sat_pos_matlab,zeros(length(pseudoranges),1),sat_vel_matlab);
    rcvr_pos = lla2eci([lla(1) lla(2) lla(3)],utcTime_array);

    clock_bias = 0;

    % [rcvr_pos, rcvr_vel, clock_bias] = gps_receiver_position(pseudoranges', satellite_positions, clk_bias_corr);

    % gps_receiver_position(pseudoranges', satellite_positions, clk_bias_corr);

end