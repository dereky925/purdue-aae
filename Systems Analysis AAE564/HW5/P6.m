clear;clc;close all

% Define symbolic variables
syms m0 m1 m2 l1 l2 theta1_e theta2_e g real
% Define the matrix M (from the provided image)
M = [m0 + m1 + m2, -m1*l1*cos(theta1_e), -m2*l2*cos(theta2_e);
     -m1*l1*cos(theta1_e), m1*l1^2, 0;
     -m2*l2*cos(theta2_e), 0, m2*l2^2];

% Invert the matrix
M_inv = inv(M);
M_ = [0 0 0;
      0 -m1*l1*g*cos(theta1_e) 0;
      0 0 -m2*l2*g*cos(theta2_e)];

M_final = M_inv * M_;

b = inv(M) * [1; 0; 0];

A_top = [zeros(3), eye(3)];
A_bottom = [M_final, zeros(3)];
A = [A_top; A_bottom];

M_final = A;

% b_sub = subs(b, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 1, 1, 0, 0})

% disp(M_final);

%8.b.
% L7 = subs(M_final, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 0.5, 1, 0, 0})

% L8 = subs(M_final, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 0.5, 1, pi, pi});



% 8.a. 
% Substitute specific values
L1 = subs(M_final, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 1, 1, 0, 0});
L2 = subs(M_final, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 1, 1, pi, pi});
L3 = subs(M_final, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 0.99, 1, 0, 0});
L4 = subs(M_final, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 0.99, 1, pi, pi});
L5 = subs(M_final, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 0.5, 1, 1, 1, 0, 0});
L6 = subs(M_final, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 0.5, 1, 1, 1, pi, pi});
L7 = subs(M_final, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 0.5, 1, 0, 0});
L8 = subs(M_final, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 0.5, 1, pi, pi});

L1_eig = double(eig(L1))
L2_eig = double(eig(L2))
L3_eig = double(eig(L3))
L4_eig = double(eig(L4))
L5_eig = double(eig(L5))
L6_eig = double(eig(L6))
L7_eig = double(eig(L7))
L8_eig = double(eig(L8))

% [~,L1_vec]= eig(L1)
% [~,L2_vec]= eig(L2)
% [~,L3_vec]= eig(L3)
% [~,L4_vec]= eig(L4)
% [~,L5_vec]= eig(L5)
% [~,L6_vec]= eig(L6)
[~,L7_vec]= eig(double(L7))
[~,L8_vec]= eig(double(L8))


