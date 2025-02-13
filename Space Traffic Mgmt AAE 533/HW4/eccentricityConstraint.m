function [rho_dot, rho] = eccentricityConstraint ( RA, dec, RA_prop, dec_prop, R_ECI, R_dot)



R = 6371;           %[km] %Radius of Earth
mu = 3.986004418e5; %[km^3/s^2] Gravitational parameter

%convert to radians
RA = RA * pi/180;
RA_prop = RA_prop * pi/180;
dec = dec * pi/180;
dec_prop = dec_prop * pi/180;


%Position Unit Vectors
up = [cos(RA)*cos(dec);...
      sin(RA)*cos(dec);...
      sin(dec)];
 
ua = [-sin(RA)*cos(dec);...
      cos(RA)*cos(dec);...
      0];
  
ud = [-cos(RA)*sin(dec);...
      -sin(RA)*sin(dec);...
      cos(dec)];

%Range
n = 1e5;
rho = linspace(1,2*R,n);
RA_dot = sum(diff(RA_prop))/length(RA_prop);
RA_dot = 3.1893e-04;
dec_dot = sum(diff(dec_prop))/length(dec_prop);
dec_dot = -7.9012e-05;
  
%
h1 = cross(R_ECI, up);
h2 = cross(up, (RA_dot*ua + dec_dot*ud).');
h3 = cross(up, R_dot) + cross(R_ECI, (RA_dot*ua + dec_dot*ud).');
h4 = cross(R_ECI, R_dot);

c0 = norm(h1)^2;
c1 = dot(2*h1, h2);
c2 = dot(2*h1, h3);
c3 = dot(2*h1, h4);
c4 = norm(h2)^2;
c5 = dot(2*h2, h3);
c6 = dot(2*h2, h4) + norm(h3)^2;
c7 = dot(2*h3, h4);
c8 = norm(h4)^2;

w0 = norm(R_ECI)^2;
w1 = 2*dot(R_dot,up);
w2 = (RA_dot^2*cos(dec)^2) + dec_dot^2;
w3 = (2*RA_dot*dot(R_dot, ua))+(2*dec_dot*dot(R_dot,ud));
w4 = norm(R_dot)^2;
w5 = 2*dot(R_ECI,up);

F_rho = w2.*rho.^2+w3.*rho+w4-(2.*mu./sqrt(rho.^2+w5.*rho+w0));
P_rho = c1*rho.^2 + c2*rho + c3;
U_rho = c4*rho.^4 + c5*rho.^3 + c6*rho.^2 + c7*rho + c8;

a4 = c0;
a3 = P_rho + c0*w1;
a2 = U_rho + c0*F_rho + w1*P_rho;
a1 = F_rho.*P_rho + w1*U_rho;
a0 = F_rho.*U_rho + ((mu^2)*(1-.4^2));

rho_dot = zeros(4,100000);
for i = 1:length(rho)
    rho_dot(:,i) = real(roots([a4, a3(i), a2(i), a1(i), a0(i)]));
end
% rho_dot(abs(imag(rho_dot)) > 0) = nan;

end