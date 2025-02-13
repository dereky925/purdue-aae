function [receiver_position, clock_bias] = gps_receiver_position(pseudoranges, satellite_positions, sat_clk_bias)
    % Solve GPS receiver position using pseudoranges and satellite positions
    % Inputs:
    % pseudoranges: Nx1 vector of measured pseudoranges (in meters)
    % satellite_positions: Nx3 matrix of satellite positions [x, y, z] (in meters)
    % Outputs:
    % receiver_position: 1x3 vector of receiver position [x, y, z] (in meters)
    % clock_bias: scalar representing the receiver clock bias (in meters)

    % Speed of light (m/s)
    c = 299792458;

    % pseudoranges = pseudoranges + sat_clk_bias;
    
    % Number of satellites
    num_sats = size(satellite_positions, 1);
    
    % Ensure at least 4 satellites for position determination
    if num_sats < 4
        error('At least 4 satellites are required to solve for position.');
    end
    
    % Initial estimates for receiver position and clock bias
    x0 = 0;
    y0 = 0;
    z0 = 0;
    b0 = 0; % Clock bias in meters
    
    % Convergence threshold
    threshold = 1e-6;
    
    % Maximum iterations
    max_iterations = 100;
    
    % Iterative solution
    for iter = 1:max_iterations
        % Compute the geometric distances to satellites
        distances = sqrt((satellite_positions(:, 1) - x0).^2 + ...
                         (satellite_positions(:, 2) - y0).^2 + ...
                         (satellite_positions(:, 3) - z0).^2);
        
        % Predicted pseudoranges
        predicted_pseudoranges = distances + b0;
        
        % Residuals
        residuals = pseudoranges - predicted_pseudoranges;
        
        % Geometry matrix H
        % H = [-(satellite_positions(:, 1) - x0) ./ distances, ...
        %      -(satellite_positions(:, 2) - y0) ./ distances, ...
        %      -(satellite_positions(:, 3) - z0) ./ distances, ...
        %      c*ones(num_sats, 1)];

        H = [-(satellite_positions(:, 1) - x0) ./ distances, ...
         -(satellite_positions(:, 2) - y0) ./ distances, ...
         -(satellite_positions(:, 3) - z0) ./ distances, ...
         ones(num_sats, 1)];
        
        % Least squares solution
        delta = (H' * H) \ (H' * residuals);
        
        % Update estimates
        x0 = x0 + delta(1);
        y0 = y0 + delta(2);
        z0 = z0 + delta(3);
        b0 = b0 + delta(4);
        
        % Check for convergence
        if norm(delta) < threshold
            break;
        end
    end
    
    % Output results
    receiver_position = [x0, y0, z0];
    clock_bias = b0 / c; % Convert clock bias to seconds
    
    % Warn if maximum iterations reached
    if iter == max_iterations
        warning('Solution did not converge within the maximum number of iterations.');
    end

end
