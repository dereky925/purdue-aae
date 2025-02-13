function xdot = two_body_ode(t,x)

%t is not used but is needed for ode45
%x = [position ; velocity]

mu = 398600.4418; % [km^3/s^2]

% Unpack the state
r_vector = x(1:3,1); % 3x1 vector
vel_vector = x(4:6,1); % 3X1 vector

r = sqrt(r_vector(1)^2 + r_vector(2)^2 + r_vector(3)^2);

%2BEOM
rddot = -mu/r^3 * r_vector; % 3x1 vector

xdot = [vel_vector;rddot]; % 6x1 vector

end