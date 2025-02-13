% rectangular_dome_ray_refraction.m
% This script performs a 3D ray‐tracing analysis through a dome‐shaped window 
% whose base is a rectangle of fixed length and width, but whose curved (bulged)
% surface is fitted to the rectangle and whose maximum height (H) is swept parametrically.
% The outer dome surface is defined by the implicit equation
%    F_outer(x,y,z) = z - H * cos(pi*x/L) * cos(pi*y/W) = 0,
% and the inner dome surface (representing a window of thickness d_thick) is defined by
%    F_inner(x,y,z) = z - (H - d_thick) * cos(pi*x/L) * cos(pi*y/W) = 0.
% The rectangular base spans x in [-L/2, L/2] and y in [-W/2, W/2], so that at the edges
% (x = ±L/2 or y = ±W/2) cos(pi/2) = 0 and hence z = 0. At the center (x = 0, y = 0) the 
% dome attains its maximum height, z = H (outer surface) or H - d_thick (inner surface).
% An incident ray (in air, n = 1) is launched from a fixed origin and directed toward the dome.
% When the ray encounters the outer surface it is refracted into the dome material (sapphire, n = 1.77)
% using Snell’s law (via refractRay). The ray then propagates inside until it meets the inner surface,
% where it is refracted back into air. The angular deviation between the incident and emergent rays
% is computed. For each dome height H (from the array H_values), the script finds the ray–surface 
% intersections numerically (using fzero) and displays the resulting geometry and ray paths in a 3D subplot.

clear; close all; clc;

% Physical parameters
n_air = 1.0;
n_sapphire = 1.77;
d_thick = 5;         % Dome (window) thickness in mm
L = 200;             % Base rectangle length in mm (x-direction)
W = 150;             % Base rectangle width in mm (y-direction)

% Dome heights to sweep (H is the maximum height at the center of the dome)
H_values = [10, 20, 30, 40, 50];  % in mm

% Incident ray parameters in 3D
% The rectangular dome is centered at (0,0) in the x-y plane.
% Choose an incident ray that will hit the dome (preferably in the first quadrant).
rayOrigin = [0, 10, 100];  % (mm)
% Define a target point roughly toward the dome center; here we choose a point with a modest z value.
target = [0, 0, 0];

rayDir = target - rayOrigin;
rayDir = rayDir / norm(rayDir);

% Create figure for subplots.
Nplots = numel(H_values);
figure('Color','white','Name','Rectangular Dome: Varying Height','NumberTitle','off');
set(gcf, 'Position', [100, 100, 1500, 1000]);

for i = 1:Nplots
    H = H_values(i);
    % 1. Compute intersection of the incident ray with the outer dome surface.
    % The outer surface is defined by F_outer(P) = z - H*cos(pi*x/L)*cos(pi*y/W) = 0.
    F_outer = @(t) domeFunction(rayOrigin + t*rayDir, H, L, W, 'outer', d_thick);
    t_guess = 100; % initial guess for t
    t_outer = fzero(F_outer, t_guess);
    P_outer = rayOrigin + t_outer * rayDir;
    
    % Compute the outer surface normal at P_outer.
    n_outer = domeNormal(P_outer, H, L, W, 'outer', d_thick);
    % Ensure the normal points from air into the dome.
    if dot(rayDir, n_outer) > 0
        n_outer = -n_outer;
    end
    
    % 2. Refract the ray into the dome (air -> sapphire).
    d_inside = refractRay(rayDir, n_outer, n_air, n_sapphire);
    
    % 3. Compute intersection of the ray (inside the dome) with the inner dome surface.
    % The inner surface is defined by F_inner(P) = z - (H - d_thick)*cos(pi*x/L)*cos(pi*y/W) = 0.
    F_inner = @(t) domeFunction(P_outer + t*d_inside, H, L, W, 'inner', d_thick);
    t_inner = fzero(F_inner, 50);
    P_inner = P_outer + t_inner * d_inside;
    
    % Compute the inner surface normal at P_inner.
    n_inner = domeNormal(P_inner, H, L, W, 'inner', d_thick);
    % For refraction from sapphire to air, the normal should point out of the dome.
    if dot(d_inside, n_inner) > 0
        n_inner = -n_inner;
    end
    
    % 4. Refract the ray out of the dome (sapphire -> air).
    d_exit = refractRay(d_inside, n_inner, n_sapphire, n_air);
    
    % 5. Compute the angular deviation (in degrees) between the incident and emergent rays.
    deviation_deg = acosd(dot(rayDir, d_exit));
    
    % 6. Plot the geometry and ray paths.
    subplot(2, ceil(Nplots/2), i);
    hold on; grid on; view(3); axis equal;
    
    % Plot the outer dome surface over the rectangular base.
    ngrid = 50;
    x_vals = linspace(-L/2, L/2, ngrid);
    y_vals = linspace(-W/2, W/2, ngrid);
    [X, Y] = meshgrid(x_vals, y_vals);
    Z_outer = H * cos(pi*X/L) .* cos(pi*Y/W);
    surf(X, Y, Z_outer, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'b');
    
    % Plot the inner dome surface.
    Z_inner = (H - d_thick) * cos(pi*X/L) .* cos(pi*Y/W);
    surf(X, Y, Z_inner, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'r');
    
    % Plot the incident ray (from rayOrigin to P_outer).
    t_vals = linspace(0, t_outer, 100);
    ray_in = rayOrigin + (t_vals.' * rayDir);
    plot3(ray_in(:,1), ray_in(:,2), ray_in(:,3), 'k', 'LineWidth', 2);
    
    % Plot the ray inside the dome (from P_outer to P_inner).
    t_vals2 = linspace(0, t_inner, 100);
    ray_mid = P_outer + (t_vals2.' * d_inside);
    plot3(ray_mid(:,1), ray_mid(:,2), ray_mid(:,3), 'g', 'LineWidth', 2);
    
    % Plot the emergent ray (from P_inner onward).
    t_vals3 = linspace(0, 50, 100);
    ray_out = P_inner + (t_vals3.' * d_exit);
    plot3(ray_out(:,1), ray_out(:,2), ray_out(:,3), 'm', 'LineWidth', 2);
    
    % Mark the intersection points.
    plot3(P_outer(1), P_outer(2), P_outer(3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    plot3(P_inner(1), P_inner(2), P_inner(3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    
    % Annotate the subplot.
    title(sprintf('H = %.0f mm, dev = %.2f°', H, deviation_deg), 'FontSize', 16);
    xlabel('x (mm)', 'FontSize', 14);
    ylabel('y (mm)', 'FontSize', 14);
    zlabel('z (mm)', 'FontSize', 14);
    xlim([-L, L]);
    ylim([-W, W]);
    zlim([0, max(H_values)+20]);
end

% Local function: domeFunction computes the implicit function value for the dome surface.
% For type 'outer': F(x,y,z) = z - H*cos(pi*x/L)*cos(pi*y/W) = 0.
% For type 'inner': F(x,y,z) = z - (H - d_thick)*cos(pi*x/L)*cos(pi*y/W) = 0.
function fval = domeFunction(P, H, L, W, type, d_thick)
    x = P(1); y = P(2); z = P(3);
    switch type
        case 'outer'
            fval = z - H * cos(pi*x/L) * cos(pi*y/W);
        case 'inner'
            fval = z - (H - d_thick) * cos(pi*x/L) * cos(pi*y/W);
        otherwise
            error('Unknown dome type.');
    end
end

% Local function: domeNormal computes the unit normal of the dome surface at point P.
% For the outer surface: F(x,y,z) = z - H*cos(pi*x/L)*cos(pi*y/W).
% Its partial derivatives (assuming the ray hits in a region where x and y are nonnegative) are:
%   dF/dx = H*(pi/L)*sin(pi*x/L)*cos(pi*y/W)
%   dF/dy = H*(pi/W)*cos(pi*x/L)*sin(pi*y/W)
%   dF/dz = 1.
% For the inner surface, H is replaced by (H - d_thick).
function n = domeNormal(P, H, L, W, type, d_thick)
    x = P(1); y = P(2); 
    switch type
        case 'outer'
            A = H;
        case 'inner'
            A = H - d_thick;
        otherwise
            error('Unknown dome type.');
    end
    dFdx = A * (pi/L) * sin(pi*x/L) * cos(pi*y/W);
    dFdy = A * (pi/W) * cos(pi*x/L) * sin(pi*y/W);
    dFdz = 1;
    gradF = [dFdx, dFdy, dFdz];
    n = gradF / norm(gradF);
end

% Local function: refractRay computes the refracted ray direction using Snell's law.
% d: incident unit vector.
% n: unit normal at the interface (pointing from medium 1 into medium 2).
% n1: refractive index of the incident medium.
% n2: refractive index of the refracted medium.
function d_refracted = refractRay(d, n, n1, n2)
    if dot(d, n) > 0
        n = -n;
    end
    eta = n1 / n2;
    cosThetaI = -dot(d, n);
    sin2ThetaT = eta^2 * (1 - cosThetaI^2);
    if sin2ThetaT > 1
        d_refracted = [NaN, NaN, NaN];
        return;
    end
    cosThetaT = sqrt(1 - sin2ThetaT);
    d_refracted = eta * d + (eta * cosThetaI - cosThetaT) * n;
    d_refracted = d_refracted / norm(d_refracted);
end