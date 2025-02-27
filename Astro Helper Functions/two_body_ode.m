function xdot = two_body_ode(~,x)

% Unpack the state
r_vector = x(1:3,1); % 3x1 vector
vel_vector = x(4:6,1); % 3X1 vector

r = sqrt(r_vector(1)^2 + r_vector(2)^2 + r_vector(3)^2);

%2BEOM
rddot = -MU/r^3 * r_vector; % 3x1 vector

xdot = [vel_vector;rddot]; % 6x1 vector

end