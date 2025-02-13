% read in a TLE object, retrive osculating state 
%
%
%
% Author: C. Frueh
% creation date  9/5/2017
% last modified: 8/25/2019
% 
% inputs: TLE file
% dependencies: vallado subroutines as uploaded on BB
% - twoline2rv and its dependencies
% - sgp4 and its dependencies
% 
% error corrections:
%   - adhere more closely to pep8, 09/01/2022
%
%

clear;clc; close all;

mu = 3.986e14;
RE = 6.371E6; % Earth Radius [m]

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
                [satrec, r(i,:), v(i,:)] = sgp4(satrec, tsince); 
            end
        

            plot3(r(end,1)*1000, r(end,2)*1000, r(end,3)*1000, 'r.', 'MarkerSize', 20); 



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
        
        
            [r_ijk, v_ijk] = keplerian2ijk(a,ecc,incl,RAAN,arg_per,nu);
    
    
            alt = (norm([r_ijk(1), r_ijk(2), r_ijk(3)]) - RE)/1000;
    
    
            options = odeset('RelTol', 1e-8,'AbsTol',1e-10);

            [T, Z] = ode45(@two_body_ode,t,[r_ijk, v_ijk],options);

            % [T, Z] = ode45(@(t,x) four_body_ode(t,x, date),...
            %         t, [r_ijk v_ijk], options);

            plot3(Z(end,1), Z(end,2), Z(end,3), 'b*', 'MarkerSize', 30);  

        end

    end

% Plot sgp4
plot3(r(:,1)*1000, r(:,2)*1000, r(:,3)*1000,'r', 'LineWidth', 2); 

% Plot two body
plot3(Z(:,1), Z(:,2), Z(:,3), 'b--', 'LineWidth', 2); 
legend('','','','Two-body ODE45','SGP4')
ax=gca;
ax.FontSize = 16;


r_m = r .* 1000;
v_m = v.*1000;
diff = r_m - Z(1:end-1,1:3);
diff2 = v_m - Z(1:end-1,4:6);

b = figure(2);
subplot(2,1,1)
hold on
grid minor
plot(diff,'LineWidth',2)
title("SGP4 vs. ODE45 Difference ECI")
xlabel('Seconds [s]')
ylabel('Error [m]')
legend('X-Position Diff','Y-Position Diff','Z-Position Diff','Location','best')
ax = gca;
ax.FontSize = 16;


subplot(2,1,2)
hold on
grid minor
plot(diff2,'LineWidth',2)
legend('X-Velocity Diff','Y-Velocity Diff','Z-Velocity Diff','Location','best')
ax = gca;
ax.FontSize = 16;
b.Position = [10 10 1400 800];








