function dx = linearizedSinglePendulum(t,x)

% m0 = 2; % Mass of cart
% m = 1; % Mass of pendulum
% l = 1; % Length of rod 
% g  = 1; % Gravity

y = x(1);
y_dot = x(2);
th = x(3);
th_dot = x(4);

% A = [0 -0.5; 0 -1.5];
% B = [0.5; 1.5];

% A = [zeros(2,2), eye(2,2); [0 -0.5; 0 -1.5] zeros(2,2)];
% B = [0; 0; 0.5; 1.5];
% 
% dx = A * x + B * u;

% Input
k = 1;
c = 2;
omega = 1000;
w = sin(omega*t);
u = -k*y - c*y_dot + w;
% u = 0.000;

y_ddot = -0.5*x(3) + u/2;
th_ddot = -1.5*x(3) + 1.5*u;


dx = [y_dot; y_ddot; th_dot; th_ddot];


end