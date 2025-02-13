clear;clc;close all

mu = 3.9860044e14;
r = [ -2949.07357  4820.806427 -3779.9543]*1000;

F = [0 0 0 1 0 0; 
     0 0 0 0 1 0; 
     0 0 0 0 0 1; 
     -mu/norm(r)^3 0 0 0 0 0; 
     0 -mu/norm(r)^3 0 0 0 0; 
     0 0 -mu/norm(r)^3 0 0 0];

% F = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];

F = [zeros(3) eye(3);
    3*mu*r(1)^2/norm(r)^5-mu/norm(r)^3 3*mu*r(1)*r(2)/norm(r)^5 3*mu*r(1)*r(3)/norm(r)^5 0 0 0;
    3*mu*r(1)*r(2)/norm(r)^5 3*mu*r(2)^2/norm(r)^5-mu/norm(r)^3 3*mu*r(2)*r(3)/norm(r)^5 0 0 0;
    3*mu*r(1)*r(3)/norm(r)^5 3*mu*r(2)*r(3)/norm(r)^5 3*mu*r(3)^2/norm(r)^5-mu/norm(r)^3 0 0 0]

expm(F*1)














