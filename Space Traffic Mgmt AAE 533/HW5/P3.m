clear;clc;close all

syms x y z r R rho rho_norm RA dec X Y Z real

r = [x;y;z]; % satellite position
R = [X;Y;Z]; % station position

rho= r-R;
rho_norm = sqrt(rho(1)^2 + rho(2)^2 + rho(3)^2);

RA = atan2(rho(2), rho(1));
dec = asin(rho(3) / rho_norm);

dRA_dx = diff(RA, x)
dRA_dy = diff(RA, y)
dRA_dz = diff(RA, z)

dRA_X = diff(RA, X)
dRA_Y = diff(RA, Y)
dRA_Z = diff(RA, Z)


ddec_dx = diff(dec, x);
ddec_dy = diff(dec, y);
ddec_dz = diff(dec, z);

ddec_dX = diff(dec, X);
ddec_dY = diff(dec, Y);
ddec_dZ = diff(dec, Z);

pretty(ddec_dx)
pretty(ddec_dy)
pretty(ddec_dz)

pretty(ddec_dX)
pretty(ddec_dY)
pretty(ddec_dZ)

