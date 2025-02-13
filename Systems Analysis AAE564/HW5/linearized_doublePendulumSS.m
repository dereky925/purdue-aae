function dx = doublePendulumSS(t,x)

m0 = 2; % Mass of cart
m1 = 1; % Mass of pendulum 1
m2 = 1; % Mass of pendulum 2
l1 = 1; % Length of rod 1
l2 = 0.5; % Length of rod 2
g  = 1; % Gravity

y = x(1);
y_dot = x(2);
theta1_e = x(3);
th1_dot = x(4);
theta2_e = x(5);
th2_dot = x(6);

% A = [   m0+m1+m2         m1*l1*cos(th1)     m2*l2*cos(th2);
%      -m1*l1*cos(th1)        m1*l1^2               0;
%      -m2*l2*cos(th2)           0                 m2*l2^2];
% 
% B = -[m1*l1*sin(th1)*th1_dot^2 + m2*l2*sin(th2)*th2_dot^2;
%                     m1*l1*g*sin(th1);
%                     m2*l2*g*sin(th2)];

A = [0, -(g*l1*m1*cos(theta1_e)^2)/(- l1*m1*cos(theta1_e)^2 - l1*m2*cos(theta2_e)^2 + l1*m0 + l1*m1 + l1*m2),  -(g*l2*m2*cos(theta2_e)^2)/(- l2*m1*cos(theta1_e)^2 - l2*m2*cos(theta2_e)^2 + l2*m0 + l2*m1 + l2*m2);
     0, -(g*l1*m1*cos(theta1_e)*(- m2*cos(theta2_e)^2 + m0 + m1 + m2))/(l1^2*m1^2 + l1^2*m0*m1 + l1^2*m1*m2 - l1^2*m1^2*cos(theta1_e)^2 - l1^2*m1*m2*cos(theta2_e)^2),  -(g*l2*m2*cos(theta1_e)*cos(theta2_e)^2)/(- l1*l2*m1*cos(theta1_e)^2 - l1*l2*m2*cos(theta2_e)^2 + l1*l2*m0 + l1*l2*m1 + l1*l2*m2);
     0,  -(g*l1*m1*cos(theta1_e)^2*cos(theta2_e))/(- l1*l2*m1*cos(theta1_e)^2 - l1*l2*m2*cos(theta2_e)^2 + l1*l2*m0 + l1*l2*m1 + l1*l2*m2), -(g*l2*m2*cos(theta2_e)*(- m1*cos(theta1_e)^2 + m0 + m1 + m2))/(l2^2*m2^2 + l2^2*m0*m2 + l2^2*m1*m2 - l2^2*m2^2*cos(theta2_e)^2 - l2^2*m1*m2*cos(theta1_e)^2)];

B = [1/(- m1*cos(theta1_e)^2 - m2*cos(theta2_e)^2 + m0 + m1 + m2);
cos(theta1_e)/(- l1*m1*cos(theta1_e)^2 - l1*m2*cos(theta2_e)^2 + l1*m0 + l1*m1 + l1*m2);
cos(theta2_e)/(- l2*m1*cos(theta1_e)^2 - l2*m2*cos(theta2_e)^2 + l2*m0 + l2*m1 + l2*m2)];

q = A\B

dx = [y_dot; q(1); th1_dot; q(2); th2_dot; q(3)];

end