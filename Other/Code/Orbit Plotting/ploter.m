R_E = 6378.1363; % [km]

earthy(R_E, "Earth", 0.25); hold on;
axis equal

%plotOrbit(x_cartesian, 'r', 0.3, 1); hold on;
%plotOrbit(x_keplerian_cartesian, 'g', 0.3, 1); hold on;
%plotOrbit(x_modified_equinoctial_cartesian, 'b', 0.3, 1); hold on;
%plotOrbit(x_milankovitch_cartesian, 'c', 0.3, 1); hold off;

%legend("", "", "Cartesian", "", "Keplerian", "", "Modified Equinoctial", "", "Milankovitch")

%%
% Sat altitude
plot(vecnorm(x_cartesian(:, 1:3), 2, 2) - R_E)

function [] = plotOrbit(x, color, scale, grade)
    plot3(x(:, 1),x(:, 2),x(:, 3), color); hold on;
    plotOrbitWithArrows(x(:, 1),x(:, 2),x(:, 3), size(x, 1) / 50, color, scale, grade);
end

function [] = plotOrbitWithArrows(x, y, z, n, color, scale, grade)
    
    % Calculate velocity components (derivatives of position)
    vx = grade*gradient(x);  % Velocity in x direction (approximate derivative)
    vy = grade*gradient(y);  % Velocity in y direction (approximate derivative)
    vz = grade*gradient(z);
    
    hold on;
    
    % Add arrowheads every n points
    idx = round(1:n:length(x));  % Select points every n points for arrow placement
    
    % Plot arrowheads only (no body)
    quiver3(x(idx), y(idx), z(idx), vx(idx), vy(idx), vz(idx), scale, color, 'MaxHeadSize', 1, 'AutoScale', 'off');  % Set 0 for arrow body size
    
end