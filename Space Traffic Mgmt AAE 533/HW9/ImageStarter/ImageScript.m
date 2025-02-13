% Little Routine to read in a fits file, crop to the object, display%%
%
% Author: Nathan Houtz
% date: 2021-09-18
% dependencies: 
% 00095337_labeled_stars.png for display
% 00095337.fit fits image
%
%

format shortg
format compact
close all
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants

MU_E = 398600.4418; %[km3 s-2], gravitational parameter
R_E = 6378.137; %[km] mean earth radius
omega_E = [0;0;7.2921150e-5]; %[rad s-1] rotation rate of Earth

% Telescope information:
lat = 32.903342; %[deg] site latitude
lon = -105.529344; %[deg] site longitude
alt = 2.225; %[km] site altitude above sea level

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load a FITS image

fname = '00095337.fit';
fInfo = fitsinfo(fname);
img = fitsread(fname);

% Crop the image to show just the object:
img_cropped = img(1980:2030,1720:1780);

% Load the labeled image
img_labeled = imread('00095337_labeled_stars.png');
img_labeled = img_labeled(102:863,605:1363,:);

% Get rid of "hot" pixels (cosmic rays, disfunctional pixels)
max_acceptable_value = 1300;
img(img>max_acceptable_value) = max_acceptable_value;


% Plot the images
f1 = figure();
tgroup1 = uitabgroup('Parent',f1);
tab(1) = uitab('Parent', tgroup1, 'Title', 'Raw image');
ax(1) = axes('parent',tab(1));
imagesc(img)
axis equal
axis([0,size(img,2),0,size(img,1)]+0.5)
colormap(gray(256));
xlabel('x [px]')
ylabel('y [px]')
title('Raw image of 29486')

tab(2) = uitab('Parent', tgroup1, 'Title', 'cropped image');
ax(2) = axes('parent',tab(2));
imagesc(img_cropped)
axis equal
axis([0,size(img_cropped,2),0,size(img_cropped,1)]+0.5)
colormap(gray(256));
xlabel('x [px]')
ylabel('y [px]')
title('Cropped image of 29486')

tab(3) = uitab('Parent', tgroup1, 'Title', 'Matched Stars');
ax(3) = axes('parent',tab(3));
imagesc(img_labeled);
axis equal
axis([0,size(img_labeled,2),0,size(img_labeled,1)]+0.5)
xlabel('x [px]')
ylabel('y [px]')
title('Background stars matched and labeled')

f1.Position = [0 0 1000 1000];


% Image Processing ======================================================

img = img_cropped;

% Compute background using median
bg_median = median(img(:));

% Compute background using mean
bg_mean = mean(img(:));

% Subtract out median and make negative values 0
img_median = img - bg_median;
img_median(img_median<0) = 0;

% Subtract out mean and make negative values 0
img_mean = img - bg_mean;
img_mean(img_mean<0) = 0;

% Compute centroid (median) - Center of Light Method
[x_center_median, y_center_median] = center_of_light(img_median);

% Compute centroid (mean) - Center of Light Method
[x_center_mean, y_center_mean] = center_of_light(img_mean);



% Define constants and settings:
res = size(img_cropped);
sigma = 1; % [px] standard deviation of Gaussian
covariance = [sigma^2, 0; 0, sigma^2]; % [px^2] covariance of Gaussian
imres = [res(2), res(1)]; % [px] image resolution (sz_x, sz_y)
mu = [x_center_mean, y_center_mean]; % [px] mean location for the Gaussian center (x, y)
A_Gauss = 20; % [] arbitrary amplitude of the Gaussian
npts_pdf = 100; % Number of points for PDF grid
rng(222); % Set random seed for reproducibility

% Calculate the integral over each square pixel area in the image:
[xc, yc] = meshgrid(1:imres(1), 1:imres(2)); % [px] (sz_y, sz_x) centers of pixels
xyL = [xc(:), yc(:)] - 0.5; % [px] (npix x 2) list of x-y pairs giving lower values
xyU = xyL + 1; % [px] (npix x 2) list of x-y pairs giving upper values
Gauss_integral = mvncdf(xyL, xyU, mu, covariance); % [] (npix x 2) integrate pixel squares
Gauss_integral = reshape(Gauss_integral, imres(2), imres(1)); % [] (sz_y, sz_x)

V_Gauss = A_Gauss * sqrt(det(2 * pi * covariance)); % [] total volume of the Gaussian
E_img = V_Gauss * Gauss_integral; % [] expected values, still (npix x 2)
% img = poissrnd(E_img);

% Generate a surface for the pdf:
xvec = linspace(mu(1) - 3 * sigma, mu(1) + 3 * sigma, npts_pdf);
yvec = linspace(mu(2) - 3 * sigma, mu(2) + 3 * sigma, npts_pdf);
[xmat, ymat] = meshgrid(xvec, yvec);
pdf = reshape(mvnpdf([xmat(:), ymat(:)], mu, covariance), npts_pdf, npts_pdf);
pdf = V_Gauss * pdf;

% Plot 
a = figure();
ax(1) = subplot(2, 1, 1);
hold on
imagesc(img);
scatter(x_center_mean,y_center_mean,100,'+','MarkerFaceColor','b')
axis equal
axis([0,size(img,2),0,size(img,1)]+0.5)
colormap(ax(1), gray(256));
xlabel('x [px]');
ylabel('y [px]');
title('Cropped Image Centroid');
legend('Center of Light Approx.')
ax = gca;
ax.FontSize = 16;
set(gca, 'YDir', 'reverse'); % Reverse the direction of the y-axis


% Create plane for background mean for plotting

% Define the x and y ranges for the plane
x = linspace(-10, res(2)+10, 50); 
y = linspace(-10, res(1)+10, 50);
[X, Y] = meshgrid(x, y);

% Define the specific z-height for the plane
z_height = bg_mean; % Replace with your desired z-height
Z = z_height * ones(size(X)); % Create a constant z plane


% Display a bar plot of pixel intensities:
ax(2) = subplot(2, 1, 2);

bplot = bar3(img, 1);
% Exclude the bar plot from the legend
for k = 1:numel(bplot)
    bplot(k).HandleVisibility = 'off'; % Set HandleVisibility for each bar
end

hold on
surf(X, Y, Z, 'FaceColor', [0, 1, 0], ... % Green color
     'FaceAlpha', 0.5, ... % Transparency (0 = fully transparent, 1 = opaque)
     'EdgeColor', 'none', ... % Remove edge lines
     'DisplayName', 'Mean Background Level = ' + string(bg_mean)); % Legend entry for the plane

% pl = surf(xvec, yvec, pdf, 'LineStyle', 'none', 'FaceColor', [0.2, 0.4, 1], ...
%     'FaceAlpha', 0.4);

AR = daspect(gca);
daspect([AR(1), AR(1), AR(3)]);
for k = 1:numel(bplot)
    zdata = bplot(k).ZData;
    bplot(k).CData = zdata;
end
colormap(ax(2), gray(256));

xlim([0,imres(1)+10])
ylim([-10,imres(2)])
xlabel('x [px]')
ylabel('y [px]')
zlabel('Pixel intensities [counts]')
view(-45, 10);
legend show;

a.Position = [1000 0 1000 1000];
ax = gca;
ax.FontSize = 16;



%% Backup

% Define constants and settings:
sigma = 1; % [px] standard deviation of Gaussian
covariance = [sigma^2, 0; 0, sigma^2]; % [px^2] covariance of Gaussian
imres = [20, 12]; % [px] image resolution (sz_x, sz_y)
mu = [8, 7]; % [px] mean location for the Gaussian center (x, y)
A_Gauss = 20; % [] arbitrary amplitude of the Gaussian
npts_pdf = 100; % Number of points for PDF grid
rng(222); % Set random seed for reproducibility

% Calculate the integral over each square pixel area in the image:
[xc, yc] = meshgrid(1:imres(1), 1:imres(2)); % [px] (sz_y, sz_x) centers of pixels
xyL = [xc(:), yc(:)] - 0.5; % [px] (npix x 2) list of x-y pairs giving lower values
xyU = xyL + 1; % [px] (npix x 2) list of x-y pairs giving upper values
Gauss_integral = mvncdf(xyL, xyU, mu, covariance); % [] (npix x 2) integrate pixel squares
Gauss_integral = reshape(Gauss_integral, imres(2), imres(1)); % [] (sz_y, sz_x)

V_Gauss = A_Gauss * sqrt(det(2 * pi * covariance)); % [] total volume of the Gaussian
E_img = V_Gauss * Gauss_integral; % [] expected values, still (npix x 2)
img = poissrnd(E_img);

% Generate a surface for the pdf:
xvec = linspace(mu(1) - 3 * sigma, mu(1) + 3 * sigma, npts_pdf);
yvec = linspace(mu(2) - 3 * sigma, mu(2) + 3 * sigma, npts_pdf);
[xmat, ymat] = meshgrid(xvec, yvec);
pdf = reshape(mvnpdf([xmat(:), ymat(:)], mu, covariance), npts_pdf, npts_pdf);
pdf = V_Gauss * pdf;

% Plot the image on the left:
a = figure();
ax(1) = subplot(1, 2, 1);
imagesc(img);
axis equal
xlim([0, imres(1)] + 0.5);
ylim([0, imres(2)] + 0.5);
colormap(ax(1), gray(256));
xlabel('x [px]');
ylabel('y [px]');
title('Simulated image of a Gaussian distribution');

% Display a bar plot of pixel intensities:
ax(2) = subplot(1, 2, 2);
bplot = bar3(img, 1);
hold on
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

a.Position = [1000 0 1000 1000];


