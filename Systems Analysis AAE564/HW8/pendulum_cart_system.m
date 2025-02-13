function dx = pendulum_cart_system(t,x)

% Parameters
m0 = 2;    % Mass of the cart
m = 1;    % Mass of the pendulum
l = 1;    % Length of the pendulum
g = 1; % Gravitational acceleration

% State variables
y         = x(1);    % Position of the cart
y_dot     = x(2);    % Velocity of the cart
theta     = x(3);    % Angle of the pendulum
theta_dot = x(4);    % Angular velocity of the pendulum

% Input
k = 0.5;
c = 5;
a = 0.1;
omega = 1;
w = a*sin(omega*t);
u = -k*y - c*y_dot + w;
% u = 0.000;

% y_ddot = l/cos(theta) * (-(m0+m)*g*tan(theta)-m*l*sin(theta)*theta_dot^2+u)/((m0+m*l/cos(theta))-m*l*cos(theta)) + g*tan(theta);
% theta_ddot = (-(m0+m)*g*tan(theta-m*l*sin(theta)*theta_dot^2+u)) / ((m0+m*l/cos(theta))-m*l*cos(theta));

theta_ddot = ( (m * l * cos(theta) * u) / (m0 + m) - (m^2 * l^2 * cos(theta) * sin(theta) * theta_dot^2) / (m0 + m) - m * l * g * sin(theta) ) / ( m * l^2 - (m^2 * l^2 * cos(theta)^2) / (m0 + m) );
y_ddot = (u + m * l * cos(theta) * theta_ddot - m * l * sin(theta) * theta_dot^2) / (m0 + m);


% y_ddot = (u - m * g * sin(theta) * cos(theta) - m * l * sin(theta) * theta_dot^2) / (m0 + m - m * cos(theta)^2);
% theta_ddot = (cos(theta) * y_ddot - g * sin(theta)) / l;
    

dx = [y_dot; y_ddot; theta_dot; theta_ddot];

end