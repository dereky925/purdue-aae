% This MATLAB script performs a 3D ray-tracing analysis through a dome-shaped window 
% whose horizontal (base) radius is fixed while its maximum vertical height is varied. 
% The dome is defined by two surfaces that are not part of a perfect sphere. The outer 
% surface (the dome window) is given by the equation:
%       z = H * cos((pi/2) * (r/R0))
% and the inner (back) surface is defined as:
%       z = (H - d_thick) * cos((pi/2) * (r/R0))
% where r = sqrt(x^2 + y^2), H is the dome’s maximum height at the center, R0 is the 
% fixed base radius, and d_thick is the dome’s thickness. The dome is assumed to be fixed 
% on the ground (its rim is at z = 0). An incident ray is specified with a fixed origin 
% and direction; the ray is aimed toward the dome. For each dome height (each H value in a 
% swept array), the code uses a numerical method (via MATLAB’s fzero function) to find the 
% intersection of the ray with the outer dome surface by solving an implicit function 
% (domeFunction). At the intersection point, the surface normal is calculated using the 
% gradient of the dome function (via the helper function domeNormal). Snell’s law 
% (implemented in the function refractRay) is then applied at the outer interface to 
% refract the ray from air (n = 1) into the dome material (sapphire, n = 1.77). The 
% refracted ray is then propagated inside the dome until it intersects the inner surface, 
% again using fzero and domeFunction. The normal at the inner surface is computed, and 
% Snell’s law is applied a second time to refract the ray from sapphire back into air. 
% The script computes the angular deviation between the incident and emergent rays. 
% Finally, for each dome height, the code generates a 3D subplot showing the outer dome 
% surface (in blue, semitransparent), the inner dome surface (in red, semitransparent), 
% the incident ray (black), the ray path inside the dome (green), and the emergent ray 
% (magenta). This visualization helps illustrate how changing the dome’s maximum height 
% while keeping its base radius fixed affects the ray’s deviation.

% This script performs a 3D ray‐tracing analysis through a dome (window)
% whose base (horizontal) radius is fixed while its maximum height is varied.
% The outer dome surface is defined by:
%     z = H * cos( (pi/2) * (sqrt(x^2+y^2)/R0) )
% and the inner (back) surface is defined by:
%     z = (H - d_thick) * cos( (pi/2) * (sqrt(x^2+y^2)/R0) )
% where H is the dome’s maximum height (at x = y = 0) and R0 is the fixed base radius.
% The dome is fixed on the ground (z = 0 along its rim).
%
% An incident ray (in air, n = 1) is launched from a fixed point and aimed at the dome.
% At the outer interface the ray refracts into sapphire (n = 1.77), then it propagates
% until it hits the inner surface, where it refracts back into air. The angular deviation
% between the incident and emergent rays is computed.
%
% Adjust the parameters as needed.

clear; close all; clc;

% Physical parameters
n_air = 1.0;
n_sapphire = 1.77;
d_thick = 5;         % window (dome) thickness in mm
R0 = 100;            % fixed base (horizontal) radius in mm

% Dome heights to sweep (H values; note: for a hemisphere H = R0)
H_values = [10, 20, 30, 40, 50];  % in mm

% Incident ray parameters
% (Choose a ray origin outside the dome; the dome sits with rim at z=0.)
rayOrigin = [20, 0, 100];  
% Choose a fixed ray direction (for example, roughly toward the dome).
rayDir = [0.1, 0, -1];
rayDir = rayDir / norm(rayDir);

% Create figure for subplots.
Nplots = numel(H_values);
figure('Color','white','Name','Fixed Base Dome (R0 fixed); Varying Height','NumberTitle','off');
set(gcf, 'Position', [100, 100, 1500, 1000]);

for i = 1:Nplots
    H = H_values(i);
    
    % 1. --- Find intersection of the incident ray with the outer dome ---
    % Outer dome function: f_outer(x,y,z) = z - H*cos((pi/2)*(r/R0)) = 0.
    F_outer = @(t) domeFunction(rayOrigin + t*rayDir, H, R0, 'outer', d_thick);
    % Find a positive root (we assume the ray does hit the dome).
    t0 = 100;  % initial guess (adjust if needed)
    t_outer = fzero(F_outer, t0);
    P_outer = rayOrigin + t_outer*rayDir;
    
    % Compute the outer dome surface normal at P_outer.
    n_outer = domeNormal(P_outer, H, R0, 'outer', d_thick);
    % Ensure the normal points from air into the dome (i.e. against the incident ray).
    if dot(rayDir, n_outer) > 0
        n_outer = -n_outer;
    end
    
    % 2. --- Refract into the dome (air -> sapphire) ---
    d_inside = refractRay(rayDir, n_outer, n_air, n_sapphire);
    
    % 3. --- Find intersection with the inner dome ---
    % Inner dome function: f_inner(x,y,z) = z - (H - d_thick)*cos((pi/2)*(r/R0)) = 0.
    F_inner = @(t) domeFunction(P_outer + t*d_inside, H, R0, 'inner', d_thick);
    t_inner = fzero(F_inner, 50);  % initial guess; adjust if necessary
    P_inner = P_outer + t_inner*d_inside;
    
    % Compute the inner dome surface normal at P_inner.
    n_inner = domeNormal(P_inner, H, R0, 'inner', d_thick);
    % For refraction from sapphire to air, the normal should point out of the dome.
    if dot(d_inside, n_inner) > 0
        n_inner = -n_inner;
    end
    
    % 4. --- Refract out of the dome (sapphire -> air) ---
    d_exit = refractRay(d_inside, n_inner, n_sapphire, n_air);
    
    % 5. --- Compute angular deviation (in degrees) between incident and emergent rays ---
    deviation_deg = acosd(dot(rayDir, d_exit));
    
    % 6. --- Plot the geometry and ray paths ---
    subplot(2, ceil(Nplots/2), i);
    hold on; grid on; view(3); axis equal;
    
    % Generate a grid for plotting the dome surfaces.
    ngrid = 50;
    u = linspace(0, R0, ngrid);
    v = linspace(0, 2*pi, ngrid);
    [U,V] = meshgrid(u,v);
    % Outer dome surface: z = H*cos((pi/2)*(U/R0))
    X_outer = U .* cos(V);
    Y_outer = U .* sin(V);
    Z_outer = H * cos((pi/2)*(U/R0));
    surf(X_outer, Y_outer, Z_outer, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'b');
    
    % Inner dome surface: z = (H - d_thick)*cos((pi/2)*(U/R0))
    X_inner = X_outer; Y_inner = Y_outer;
    Z_inner = (H - d_thick) * cos((pi/2)*(U/R0));
    surf(X_inner, Y_inner, Z_inner, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'r');
    
    % Plot the incident ray (from rayOrigin to P_outer).
    t_vals = linspace(0, t_outer, 100);
    ray_in = rayOrigin + (t_vals.' * rayDir);
    plot3(ray_in(:,1), ray_in(:,2), ray_in(:,3), 'k', 'LineWidth', 2);
    
    % Plot the ray inside the dome (from P_outer to P_inner).
    t_vals2 = linspace(0, t_inner, 100);
    ray_mid = P_outer + (t_vals2.' * d_inside);
    plot3(ray_mid(:,1), ray_mid(:,2), ray_mid(:,3), 'g', 'LineWidth', 2);
    
    % Plot the emerging ray (from P_inner onward).
    t_vals3 = linspace(0, 50, 100);
    ray_out = P_inner + (t_vals3.' * d_exit);
    plot3(ray_out(:,1), ray_out(:,2), ray_out(:,3), 'm', 'LineWidth', 2);
    
    % Mark intersection points.
    plot3(P_outer(1), P_outer(2), P_outer(3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    plot3(P_inner(1), P_inner(2), P_inner(3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    
    % Annotate subplot.
    title(sprintf('H = %.0f mm, dev = %.2f°', H, deviation_deg), 'FontSize', 16);
    xlabel('x (mm)', 'FontSize', 14); ylabel('y (mm)', 'FontSize', 14); zlabel('z (mm)', 'FontSize', 14);
    xlim([-150, 250]); ylim([-150, 150]); zlim([0, max(H_values)+20]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local function definitions

function fval = domeFunction(P, H, R0, type, d_thick)
    % domeFunction evaluates the implicit function for the dome surface.
    % For type 'outer', the surface is defined by:
    %      f(x,y,z) = z - H*cos((pi/2)*(r/R0)) = 0.
    % For type 'inner', it is defined by:
    %      f(x,y,z) = z - (H - d_thick)*cos((pi/2)*(r/R0)) = 0.
    % P is a 3-element vector [x, y, z].
    x = P(1); y = P(2); z = P(3);
    r = sqrt(x^2 + y^2);
    switch type
        case 'outer'
            fval = z - H * cos((pi/2)*(r/R0));
        case 'inner'
            fval = z - (H - d_thick) * cos((pi/2)*(r/R0));
        otherwise
            error('Unknown dome type.');
    end
end

function n = domeNormal(P, H, R0, type, d_thick)
    % domeNormal computes the unit surface normal for the dome at point P.
    % For the dome function f(x,y,z) given in domeFunction, its gradient is:
    %    ∇f = [ A*x/r, A*y/r, 1 ], where A = {H or (H-d_thick)}*(pi/(2*R0)) * sin((pi/2)*(r/R0))
    % (with r = sqrt(x^2+y^2); if r==0, use [0,0,1]).
    x = P(1); y = P(2); z = P(3);
    r = sqrt(x^2+y^2);
    switch type
        case 'outer'
            A = H * (pi/(2*R0)) * sin((pi/2)*(r/R0));
        case 'inner'
            A = (H - d_thick) * (pi/(2*R0)) * sin((pi/2)*(r/R0));
        otherwise
            error('Unknown dome type.');
    end
    if r < 1e-6
        grad_f = [0, 0, 1];
    else
        grad_f = [A*x/r, A*y/r, 1];
    end
    n = grad_f / norm(grad_f);
end

function d_refracted = refractRay(d, n, n1, n2)
    % refractRay computes the refracted ray direction using Snell's law.
    % d: incident unit vector.
    % n: unit normal at the interface (pointing from medium 1 into medium 2).
    % n1: refractive index of the incident medium.
    % n2: refractive index of the refracted medium.
    % Returns the unit refracted ray vector.
    if dot(d, n) > 0
        n = -n;
    end
    eta = n1 / n2;
    cosThetaI = -dot(d, n);
    sin2ThetaT = eta^2 * (1 - cosThetaI^2);
    if sin2ThetaT > 1
        % Total internal reflection: return NaN vector.
        d_refracted = [NaN, NaN, NaN];
        return;
    end
    cosThetaT = sqrt(1 - sin2ThetaT);
    d_refracted = eta * d + (eta*cosThetaI - cosThetaT)*n;
    d_refracted = d_refracted / norm(d_refracted);
end

function t = findIntersection(F, t_guess)
    % A simple wrapper that uses fzero to find a root of F(t)=0 starting from t_guess.
    % (Not used separately in the main code.)
    options = optimset('TolX',1e-6);
    t = fzero(F, t_guess, options);
end