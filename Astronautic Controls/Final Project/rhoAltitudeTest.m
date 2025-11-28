h = sdpvar(1,1);      % Altitude (decision variable)
rho = sdpvar(1,1);    % Density (decision variable)
lambda = sdpvar(2,1); % Convex weights for the interval

% Suppose we know h lies between 0 and 10 km for a given interval
constraints = [];
constraints = [constraints, lambda(1) + lambda(2) == 1, lambda >= 0];
constraints = [constraints, h == lambda(1)*0 + lambda(2)*10];  % linear interpolation for altitude
constraints = [constraints, rho == lambda(1)*1.225 + lambda(2)*0.4135];
% Add additional constraints for other intervals or use binary variables to switch intervals.
