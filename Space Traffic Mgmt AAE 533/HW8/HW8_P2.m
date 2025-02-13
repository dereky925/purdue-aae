% =========================================================================
% 
% Filename:       HW8.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE533 - Space Traffic Management
% Professor:      Dr. Carolin Frueh
% Contact:        cfrueh@purdue.edu
% Assignment:     HW 8
% Semester:       Fall 2024
% 
% Description:
% Plot all objects in space-track.org catalog
%
%
% =========================================================================
clear;clc;close all

mu = 3.986e14;
RE = 6.371E6; % Earth Radius [m]

data = readtable('TLEs.xlsx');

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
    flipud(imData), 'EdgeColor', 'none');

% Adjust axes
axis equal;
view(3);  % 3D view
grid on;
xlabel('Inertial x (m)');
ylabel('Inertial y (m)');
zlabel('Inertial z (m)');

plotA.Position = [1600 100 1000 1000];

target_datetime = datetime(2024, 11, 8, 0 , 0, 0,...
                                          'TimeZone', 'America/New_York');
count = 0;
len = size(data);

intelsat33E_loc = 1E7*[-4.1407 -0.7635 0.0015];
% intelsat33E_loc = 1E7*[3.3319 -2.5707 0.0010];
%%
for i = 1:len(1)
% for i = 15100

    if contains(data.satName{i},'INTELSAT 33E')
    % if contains(data.satName{i},'NAVSTAR')
    % if contains(data.satName{i},'STARLINK')
    
        disp(string(i/len(1)*100) + '%')
        count = count + 1;
        disp(count)

        epoch = (strsplit(data.line1{i}, ' '));
        epoch = char(epoch(4));
        year = 2000 + str2double(epoch(1:2));  % Adds 2000 for 21st century (adjust for other cases)
        day_of_year = str2double(epoch(3:5));
        fractional_day = str2double(['0' epoch(6:end)]);
        date = datetime(year, 1, 1, 'TimeZone', 'UTC') + days(day_of_year - 1) + days(fractional_day);

        prop_seconds = seconds(target_datetime - date);

        orb_elements = str2double(strsplit(data.line2{i}, ' '));
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

        % nu = atand( (sind(E)*sqrt(1-ecc^2)) / (cosd(E) - ecc));
        % if nu < 0
        %     nu = nu + 360;
        % end
    
        % Compute true anomaly
        x = cosd(E) - ecc;
        y = sqrt(1-ecc^2) * sind(E);
        nu = (atan2d(y,x));
        if nu < 0
            nu = nu + 360;
        end

 
        % Compute semi-major axis 
        a = (mu / (4*pi^2*(mean_motion/86400)^2) )^(1/3);


        % num_orbits = 1;
        % P_coast = num_orbits * (2*pi)/(sqrt(mu/(a^3))); %Coast time
        t = 1:1:prop_seconds;
        % t = 1:1:2;
    
    
        [r_ijk, v_ijk] = keplerian2ijk(a,ecc,incl,RAAN,arg_per,nu);


        alt = (norm([r_ijk(1), r_ijk(2), r_ijk(3)]) - RE)/1000;

    % if alt > 34000

        options = odeset('RelTol', 1e-8,'AbsTol',1e-10);
        % options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
        [T, Z] = ode45(@two_body_ode,t,[r_ijk, v_ijk],options);
        % [T, Z] = ode45(@(t,x) four_body_ode(t,x, date),...
        %         t, [r_ijk v_ijk], options);
        
        altitude = (norm([Z(end,1), Z(end,2), Z(end,3)]) - RE)/1000;
    
        % Plot GEO satellites
        if altitude > 30000
            
            name = data.satName{i};
            % name = name(3:end-2);

            if contains(data.satName{i},'INTELSAT 33E')
                plot3(Z(end,1), Z(end,2), Z(end,3), 'pk', 'MarkerFaceColor','y', 'MarkerSize', 20);
                Z(end,1:3)
                text(Z(end,1), Z(end,2), Z(end,3), name(3:end), 'FontSize', 12, 'Color', 'k');
            else
                
                plot3(Z(end,1), Z(end,2), Z(end,3), 'r.', 'MarkerSize', 20);

                if (norm(Z(end,1:3) - intelsat33E_loc) < 5E6)
                    text(Z(end,1), Z(end,2), Z(end,3), name, 'FontSize', 12, 'Color', 'k');
                end
                
            end
        % Plot MEO satellites
        elseif altitude > 2000
            plot3(Z(end,1), Z(end,2), Z(end,3), 'b.', 'MarkerSize', 3);
        % Plot LEO satellites
        % else 
            plot3(Z(end,1), Z(end,2), Z(end,3), 'g.', 'MarkerSize', 3);  
        end

    end

end







