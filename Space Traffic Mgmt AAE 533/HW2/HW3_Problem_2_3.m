clear;clc

% Define symbolic variables
syms m0 m1 m2 l1 l2 theta1_e theta2_e g real

% Define the matrix M (from the provided image)
M = [m0 + m1 + m2, -m1*l1*cos(theta1_e), -m2*l2*cos(theta2_e);
     -m1*l1*cos(theta1_e), m1*l1^2, 0;
     -m2*l2*cos(theta2_e), 0, m2*l2^2];

% Invert the matrix
M_inv = inv(M);

M_ = [0 0 0;
      0 m1*l1*g*cos(theta1_e) 0;
      0 0 m2*l2*g*cos(theta2_e)];

M_inv = M_inv * M_;

disp(M_inv);

% Substitute specific values 
M_sub = subs(M_inv, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 1, 1, 0, 0});
% M_sub = subs(M_inv, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 1, 1, pi, pi});
% M_sub = subs(M_inv, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 0.5, 1, 0, 0});
% M_sub = subs(M_inv, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 0.5, 1, pi, pi});

% Display the numeric inverse matrix
disp(M_sub);

%

A = [0, 1/2, 1/2;
     0, 3/2, 1/2;
     0, 1/2, 3/2];

A = [0, 1/2, 1/2;
     0, -3/2, -1/2;
     0, -1/2, -3/2];

A = [0, 1/2, 1/2;
     0, 3/2, 1/2;
     0, 1, 3];

A = [0, 1/2, 1/2;
     0, -3/2, -1/2;
     0, -1, -3];

B = [1;0;0];

C = eye(3);

D = [0; 0; 0];

sys = ss(A,B,C,D);
poles = pole(sys)
zeros = tzero(sys)



