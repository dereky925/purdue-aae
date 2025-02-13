function [orbitPlot] = OrbitPlotSetup(lon)

    if nargin < 1
        lon = 0;  % Default value for 'lon'
    end

    % Inputs:
    % lon - Longitude to rotate Earth image [°]

    % Outputs:
    % orbitPlot - Figure()

    orbitPlot = figure('color','white');
    hold on
    imData = imread('no_clouds_4k.jpg');  % Load the texture image
    
    [xS, yS, zS] = sphere(50);  % Generate a sphere for the globe
    earth_radius = 6378.1370;  % Earth radius [km]
    
    xSE = earth_radius * xS;
    ySE = earth_radius * yS;
    zSE = earth_radius * zS;
    
    % Create the rotation matrices
    lon = deg2rad(lon);
    Ry = R3(-lon+1.9);
    
    Rx = R1(0);
    
    % Combine the rotations (Ry first for longitude, then Rx for latitude)
    rotationMatrix = Rx * Ry;
    
    % Apply the rotation to the coordinates
    rotated_coords = rotationMatrix * [xSE(:)'; ySE(:)'; zSE(:)'];
    
    % Reshape the rotated coords back to match the surface grid dimensions
    xSE_rot = reshape(rotated_coords(1, :), size(xSE));
    ySE_rot = reshape(rotated_coords(2, :), size(ySE));
    zSE_rot = reshape(rotated_coords(3, :), size(zSE));
    
    % Plot the Earth with texture mapping, using the rotated coordinates
    surface(xSE_rot, ySE_rot, zSE_rot, 'FaceColor', 'texturemap',...
        'CData', flipud(imData), 'EdgeColor', 'none');
    
    % Adjust axes
    axis equal;
    view(3);  % 3D view
    grid on;
    xlabel('Inertial X [km]');
    ylabel('Inertial Y [km]');
    zlabel('Inertial Z [km]');

    ax = gca;
    ax.FontSize = 20;

end