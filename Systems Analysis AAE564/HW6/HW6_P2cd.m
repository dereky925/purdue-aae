clear;clc;close all

syms Omega w
A = [0 0 1 0; 0 0 0 1; Omega^2-w^2/2 w^2/2 0 0; w^2/2 Omega^2-w^2/2 0 0];

[a,b] = eig(A)






