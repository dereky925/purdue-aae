% =========================================================================
% 
% Filename:       P1_TLE_Wrapper.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE533 - Space Traffic Management
% Professor:      Dr. Carolin Frueh
% Contact:        cfrueh@purdue.edu
% Assignment:     HW 9
% Semester:       Fall 2024
% 
% Description:
% Problem 5 Astra 2E satellite observation
%
%
% =========================================================================

clear;clc; close all;


mu = 3.986e14;
RE = 6.371E6; % Earth Radius [m]

% Create a datetime object for November 2, 2024, at 4 PM UTC
desiredUTC = datetime(2024, 11, 2, 16, 0, 0, 'TimeZone', 'UTC');

% SGP4 ====================================================================

cc=0; % set line counter

    fid = fopen('ASTRA2E.txt'); % load the TLE
    
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

            epoch = (strsplit(tline1, ' '));
            epoch = char(epoch(4));
            year = 2000 + str2double(epoch(1:2));  % Adds 2000 for 21st century (adjust for other cases)
            day_of_year = str2double(epoch(3:5));
            fractional_day = str2double(['0' epoch(6:end)]);
            date = datetime(year, 1, 1, 'TimeZone', 'UTC') + days(day_of_year - 1) + days(fractional_day);

            % Compute time to desired time in minutes
            period = minutes(desiredUTC - date); % [minutes]
        
            % how far shall the TLE be propagated [minutes]
            tsince = period; 

            % extract position and velocity
            for i = 1:period
                tsince = i;
                [satrec, r(i,:), v(i,:)] = sgp4(satrec, tsince); 
            end
        end
    end

AU = 1.496E11; % [m]

% Convert End of Orbit ECI position to ECEF
r_eci = [r(end,1)*1000 r(end,2)*1000 r(end,3)*1000];
v_eci = [v(end,1)*1000 v(end,2)*1000 v(end,3)*1000];
a_eci = [0 0 0];
[r_ecef,~,~] = eci2ecef(desiredUTC,r_eci, v_eci, a_eci);

% Congo, Africa
lla_zenith = ecef2lla(r_ecef');
lla_zenith = [lla_zenith(1) lla_zenith(2) 0];

% Convert zenith station into ECI
station_zenith_position = lla2eci(lla_zenith,[2024 11 2 16 0 0]);

% Compute rho (range from station to satellite)
rho = norm(r_eci - station_zenith_position);
% vecnorm(r_eci - station_zenith_position, 2, 2)

% Compute flat surface satellite travel function, lambertian
Tau_plate = 1/rho^2;

% Compute distance from the sun
r_sun = sun(juliandate(desiredUTC))*AU;
r_sunobject = norm(r_sun);

% I0 - Mean solar constant
Ls = 3.828E26; %[W] Mean Luminosity 
I0 = Ls / (4*pi*AU^2); % [W/m^2]

% Scale solar constant directly to a given object location
Is = I0*AU^2/r_sunobject^2;
% Is = I0; % [m] Assume Is is about equal to I0

theta_s = 45; %[deg] 
Is_plate = Is*cosd(theta_s);

I_obj = Is_plate * Tau_plate * theta_s*pi/180;

% Magnitude of Sun
mag_sun = -26.832;
M_sun_abs = 4.74;

% Magnitude of object
mag = mag_sun - 2.5*log10(Is_plate/I0);

% Compute Signal
D = 1; % [m] Aperture diameter
d = 0.06; % [m] 6 cm, seconday mirror diameter
lambda_ = 600E-9; % [m] 600 nm, Given, average visible light wavelength
h = 6.62607015E-34; % Js
c = 3E8; % [m/s]
k = 0.1; % [mag/airmass] for 600nm wavelength at sealevel
R = 1; % Normally computed with 1/cos(zenith angle)
gamma = k*R; % Atmospheric exctinction coefficient
L = 0.8; % Loss function 

S_obj = (pi*(D/2)^2 - pi*(d/2)^2)*lambda_/(h*c)*exp(gamma*R)*I_obj*L


% Compute FWHM (Full Width of Half Maximum) ===============================
r_D = D/2;
FWHM_airy = 1.028*lambda_/(2*r_D)


% Orbit Plot  =============================================================
plotA = figure();
hold on
imData = imread('2_no_clouds_4k.jpg');  % Load the texture image

[xS, yS, zS] = sphere(50);  % Generate a sphere for the globe
earth_radius = 6378.1370 * 1000;  % Earth radius in meters

xSE = earth_radius * xS;
ySE = earth_radius * yS;
zSE = earth_radius * zS;

% Create the rotation matrices
Ry = R3(lla_zenith(2)*pi/180+0.8);

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

% Mark Mission Start
scatter3(r(1,1)*1000, r(1,2)*1000, r(1,3)*1000, 200, 'pr', 'filled','MarkerEdgeColor','k'); 

% Mark Mission End
scatter3(r(end,1)*1000, r(end,2)*1000, r(end,3)*1000, 700, 'py', 'filled','MarkerEdgeColor','k'); 

% Plot Zenith Station Location
scatter3(station_zenith_position(1), station_zenith_position(2), station_zenith_position(3), 700, 'r.');

% Plot Orbit
plot3(r(:,1)*1000, r(:,2)*1000, r(:,3)*1000,'r', 'LineWidth', 2); 


% Adjust axes
axis equal;
view(3);  % 3D view
grid on;
xlabel('Inertial x (m)');
ylabel('Inertial y (m)');
zlabel('Inertial z (m)');

plotA.Position = [1600 100 1000 1000];
ax = gca();
ax.FontSize = 16;




% Define constants and settings:

% Given Image scale of 1 arcsec per pixel
% sigma = 5; % [px] standard deviation of Gaussian
FWHM_given = 2; % [arcseconds] or [pixels]
sigma = FWHM_given/(2*sqrt(2*log(2)));

covariance = [sigma^2, 0; 0, sigma^2]; % [px^2] covariance of Gaussian
imres = [30, 20]; % [px] image resolution (sz_x, sz_y)
mu = [15, 10]; % [px] mean location for the Gaussian center (x, y)
A_Gauss = 50; % [] arbitrary amplitude of the Gaussian
npts_pdf = 1000; % Number of points for PDF grid
rng(222); % Set random seed for reproducibility

int_t = 5; % [s], given integration time
Q = 0.65; % given quantum efficiency
C_all = S_obj * Q * int_t;
A_Gauss = 0.838*C_all/ (2*pi*sigma);

% Calculate the integral over each square pixel area in the image:
[xc, yc] = meshgrid(1:imres(1), 1:imres(2)); % [px] (sz_y, sz_x) centers of pixels
xyL = [xc(:), yc(:)] - 0.5; % [px] (npix x 2) list of x-y pairs giving lower values
xyU = xyL + 1; % [px] (npix x 2) list of x-y pairs giving upper values
Gauss_integral = mvncdf(xyL, xyU, mu, covariance); % [] (npix x 2) integrate pixel squares
Gauss_integral = reshape(Gauss_integral, imres(2), imres(1)); % [] (sz_y, sz_x)

V_Gauss = A_Gauss * sqrt(det(2 * pi * covariance)); % [] total volume of the Gaussian
E_img = V_Gauss * Gauss_integral; % [] expected values, still (npix x 2)
img = poissrnd(E_img);


% Signal-to-Noise Ratio
SNR = 100; 
noise_amplitude100 = A_Gauss / SNR;

SNR = 10;
noise_amplitude10 = A_Gauss / SNR;

SNR = 2;
noise_amplitude2 = A_Gauss / SNR;

noisy_img100 = img + noise_amplitude100 * abs(randn(size(img))); % Add Gaussian noise
noisy_img10 = img + noise_amplitude10 * abs(randn(size(img))); % Add Gaussian noise
noisy_img2 = img + noise_amplitude2 * abs(randn(size(img))); % Add Gaussian noise

% Generate a surface for the pdf:
xvec = linspace(mu(1) - 3 * sigma, mu(1) + 3 * sigma, npts_pdf);
yvec = linspace(mu(2) - 3 * sigma, mu(2) + 3 * sigma, npts_pdf);
[xmat, ymat] = meshgrid(xvec, yvec);
pdf = reshape(mvnpdf([xmat(:), ymat(:)], mu, covariance), npts_pdf, npts_pdf);
pdf = V_Gauss * pdf;

% Plot the image on the left:
c = figure();
title('Simulated image of a Gaussian distribution');
ax(2) = subplot(2, 2, 1);
imagesc(img);
axis equal
xlim([0, imres(1)] + 0.5);
ylim([0, imres(2)] + 0.5);
colormap(ax(2), gray(256));
xlabel('x [px]');
ylabel('y [px]');
title('No Noise');
ax = gca();
ax.FontSize = 16;

ax(2) = subplot(2, 2, 2);
imagesc(noisy_img100);
axis equal
xlim([0, imres(1)] + 0.5);
ylim([0, imres(2)] + 0.5);
colormap(ax(2), gray(256));
xlabel('x [px]');
ylabel('y [px]');
title('SNR = 100');
ax = gca();
ax.FontSize = 16;

ax(2) = subplot(2, 2, 3);
imagesc(noisy_img10);
axis equal
xlim([0, imres(1)] + 0.5);
ylim([0, imres(2)] + 0.5);
colormap(ax(2), gray(256));
xlabel('x [px]');
ylabel('y [px]');
title('SNR = 10');
ax = gca();
ax.FontSize = 16;

ax(2) = subplot(2, 2, 4);
imagesc(noisy_img2);
axis equal
xlim([0, imres(1)] + 0.5);
ylim([0, imres(2)] + 0.5);
colormap(ax(2), gray(256));
xlabel('x [px]');
ylabel('y [px]');
title('SNR = 2');
ax = gca();
ax.FontSize = 16;


c.Position = [1000 0 1000 1000];


b = figure();
% No Noise
ax(2) = subplot(2, 2, 1);
bplot = bar3(img, 1);
hold on
grid minor
pl = surf(xvec, yvec, pdf, 'LineStyle', 'none', 'FaceColor', [0.2, 0.4, 1], ...
    'FaceAlpha', 0.4);
AR = daspect(gca);
daspect([AR(1), AR(1), AR(3)]);
for k = 1:numel(bplot)
    zdata = bplot(k).ZData;
    bplot(k).CData = zdata;
end
colormap(ax(2), gray(256));
xlim([0,imres(1)] + 0.5)
ylim([0,imres(2)] + 0.5)
xlabel('x [px]')
ylabel('y [px]')
zlabel('Pixel intensities [counts]')
title('No Noise')

ax = gca();
ax.FontSize = 16;

% SNR = 100
ax(2) = subplot(2, 2, 2);
bplot = bar3(noisy_img100, 1);
hold on
grid minor
pl = surf(xvec, yvec, pdf, 'LineStyle', 'none', 'FaceColor', [0.2, 0.4, 1], ...
    'FaceAlpha', 0.4);
AR = daspect(gca);
daspect([AR(1), AR(1), AR(3)]);
for k = 1:numel(bplot)
    zdata = bplot(k).ZData;
    bplot(k).CData = zdata;
end
colormap(ax(2), gray(256));
xlim([0,imres(1)] + 0.5)
ylim([0,imres(2)] + 0.5)
xlabel('x [px]')
ylabel('y [px]')
zlabel('Pixel intensities [counts]')
title('SNR = 100')

ax = gca();
ax.FontSize = 16;

% SNR = 10
ax(2) = subplot(2, 2, 3);
bplot = bar3(noisy_img10, 1);
hold on
grid minor
pl = surf(xvec, yvec, pdf, 'LineStyle', 'none', 'FaceColor', [0.2, 0.4, 1], ...
    'FaceAlpha', 0.4);
AR = daspect(gca);
daspect([AR(1), AR(1), AR(3)]);
for k = 1:numel(bplot)
    zdata = bplot(k).ZData;
    bplot(k).CData = zdata;
end
colormap(ax(2), gray(256));
xlim([0,imres(1)] + 0.5)
ylim([0,imres(2)] + 0.5)
xlabel('x [px]')
ylabel('y [px]')
zlabel('Pixel intensities [counts]')
title('SNR = 10')

ax = gca();
ax.FontSize = 16;

% SNR = 2
ax(2) = subplot(2, 2, 4);
bplot = bar3(noisy_img2, 1);
hold on
grid minor
pl = surf(xvec, yvec, pdf, 'LineStyle', 'none', 'FaceColor', [0.2, 0.4, 1], ...
    'FaceAlpha', 0.4);
AR = daspect(gca);
daspect([AR(1), AR(1), AR(3)]);
for k = 1:numel(bplot)
    zdata = bplot(k).ZData;
    bplot(k).CData = zdata;
end
colormap(ax(2), gray(256));
xlim([0,imres(1)] + 0.5)
ylim([0,imres(2)] + 0.5)
xlabel('x [px]')
ylabel('y [px]')
zlabel('Pixel intensities [counts]')
title('SNR = 2')

b.Position = [0 0 1000 1000];
ax = gca();
ax.FontSize = 16;


%% Noise level is usually computed with RMS, but it seems safe to assume 
% that is approximately 1 sigma. Tested below:

desiredRMS = 1; % Set desired RMS value
noise = abs( desiredRMS * randn(1000, 1) ); % Generate Gaussian noise with zero mean
actualRMS = sqrt(mean(noise.^2)); % Verify the RMS
disp(actualRMS); % Should be close to 0.5




