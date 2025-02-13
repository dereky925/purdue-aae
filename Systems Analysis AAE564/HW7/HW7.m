clear;clc;close all

A = [ 0 1 0 0; 0 0 1 0; 0 0 0 1; 564 1 0 0];

A*A*A*A



%4
A = [0 1; 1 0];
t = log(2);
tA = t*A;


A = [3 -5;-5 3];
% Jordan Method
[T,V] = eig(A);
V2 = inv(T)*A*T;
syms t
jordan = T*expm(V*t)*inv(T);


% Matlab method
B = expm(A*t);


% STM method
STM = eye(2) + tA + 1/2*tA^2 + 1/factorial(3)*tA^3 + 1/factorial(4)*tA^4 + 1/factorial(5)*tA^5




