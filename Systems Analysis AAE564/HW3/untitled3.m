clear;clc;close all
% Define symbolic variables for the states, inputs, and parameters
syms y theta1 theta2 dy dtheta1 dtheta2 u
syms m0 m1 m2 l1 l2 g

% Define the equations of motion (substitute the equations given in the problem)
eq1 = (m0 + m1 + m2)*diff(y,2) - m1*l1*cos(theta1)*diff(theta1,2) - m2*l2*cos(theta2)*diff(theta2,2) + ...
       m1*l1*sin(theta1)*diff(theta1,2)^2 + m2*l2*sin(theta2)*diff(theta2,2)^2 == u;
   
eq2 = -m1*l1*cos(theta1)*diff(y,2) + m1*l1^2*diff(theta1,2) + m1*l1*g*sin(theta1) == 0;

eq3 = -m2*l2*cos(theta2)*diff(y,2) + m2*l2^2*diff(theta2,2) + m2*l2*g*sin(theta2) == 0;

% Define the state vector
x = [y; dy; theta1; dtheta1; theta2; dtheta2];

c = [m0; m1; m2; l1; l2; g];
c_1 = [2; 1; 1; 1; 1; 1];
c_1 = [2; 1; 1; 1; 0.5; 1];

% Linearize the system around E1 = (y^e = 0, theta1^e = 0, theta2^e = 0)
x_eq = [0; 0; 0; 0; 0; 0];
u_eq = 0;

% Compute the Jacobians (A = df/dx, B = df/du, C = dg/dx, D = dg/du)
A = jacobian([eq1; eq2; eq3], x);
B = jacobian([eq1; eq2; eq3], u);
C = jacobian([y; theta1; theta2], x);
D = jacobian([y; theta1; theta2], u);

% Substitute the equilibrium points and parameters
A = subs(A, [x; u], [x_eq; u_eq]);
B = subs(B, [x; u], [x_eq; u_eq]);
C = subs(C, [x; u], [x_eq; u_eq]);
D = subs(D, [x; u], [x_eq; u_eq]);

% Evaluate the state-space matrices for a given set of parameters
A_num = (subs(A, c, c_1));
B_num = (subs(B, c, c_1));
C_num = (subs(C, c, c_1));
D_num = (subs(D, c, c_1));

% Display the results
disp('A matrix:'); disp(A_num);
disp('B matrix:'); disp(B_num);
disp('C matrix:'); disp(C_num);
disp('D matrix:'); disp(D_num);


% % A = [0,   0,   0,   0,   0,   0;
% %      0,   0,   1,   0,   0,   0;
% %      0,   0,   0,   0,  1/2 , 0];
% 
% A = [0,   0,   0,   0,   0,   0;
%      0,   0,   1,   0,   0,   0;
%      0,   0,   0,   0,   1,   0];
% 
% B = [1; 0; 0];
% 
% C = [1, 0, 0, 0, 0, 0;
%      0, 0, 1, 0, 0, 0;
%      0, 0, 0, 0, 1, 0];
% 
% D = [0; 0; 0];
% 
% sys = ss(A,B,C,D);
% poles = pole(sys);
% zeros = tzero(sys);


