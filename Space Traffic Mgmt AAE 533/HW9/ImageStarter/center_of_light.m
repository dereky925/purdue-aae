function [x_center, y_center] = center_of_light(image)
    % Input: 
    % - image: 2D array representing the pixel intensities
    %
    % Output:
    % - x_center: The x-coordinate of the center of light
    % - y_center: The y-coordinate of the center of light

    % Get the dimensions of the image
    [rows, cols] = size(image);
    
    % Create coordinate grids for x and y
    [x, y] = meshgrid(1:cols, 1:rows);
    
    % Compute the total intensity (denominator of the equations)
    total_intensity = sum(image(:));
    
    % Compute the weighted sums for x and y (numerators of the equations)
    x_weighted_sum = sum(sum(image .* x));
    y_weighted_sum = sum(sum(image .* y));
    
    % Calculate the center of light
    x_center = x_weighted_sum / total_intensity;
    y_center = y_weighted_sum / total_intensity;
end