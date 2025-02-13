function [ecc_1,ecc_2,ecc_3, ecc_4, a_upper, a_lower, zero_upper, zero_lower, rho, orbits, x_random, y_random] = get_a_ecc_constraints(R,R_dot,RA,RA_dot,dec,dec_dot,numRE,t,e)

% R [m]
% R [m/s]
% RA [rad]
% RA_dot [rad/s]
% dec [rad]
% dec_dot [rad/s]

mu = 3.9860044e14; %m^3/s^2

% Rho guess vector
RE = 6371000; %[m]
rho = 0:10000:RE*numRE;

u_p = [cos(RA)*cos(dec); sin(RA)*cos(dec); sin(dec)];
u_RA = [-sin(RA)*cos(dec); cos(RA)*cos(dec); 0];
u_dec = [-cos(RA)*sin(dec); -sin(RA)*sin(dec); cos(dec)];

w_0 = norm(R)^2;
w_1 = 2*(R_dot*u_p);
w_2 = RA_dot^2*cos(dec)^2 + dec_dot^2;
w_3 = 2*RA_dot*(R_dot*u_RA) + 2*dec_dot*(R_dot*u_dec);
w_4 = norm(R_dot)^2;
w_5 = 2*(R'*u_p);

F = w_2.*rho.^2 + w_3.*rho + w_4 - ( 2.*mu ./ (sqrt(rho.^2 + w_5.*rho + w_0))  );

% Zero energy constraint
E = 0;
rho_dot_upper = -w_1/2 + (sqrt( (w_1/2)^2 - F + 2*E));
rho_dot_lower = -w_1/2 - (sqrt( (w_1/2)^2 - F + 2*E));

rho_dot_upper(abs(imag(rho_dot_upper)) > 0) = nan;
rho_dot_lower(abs(imag(rho_dot_lower)) > 0) = nan;

zero_upper = rho_dot_upper;
zero_lower = rho_dot_lower;

% Semi-major axis constraint
E = -mu/(2*RE+600000);
rho_dot_upper_a_constraint = -w_1/2 + (sqrt( (w_1/2)^2 - F + 2*E));
rho_dot_lower_a_constraint = -w_1/2 - (sqrt( (w_1/2)^2 - F + 2*E));

rho_dot_upper_a_constraint(abs(imag(rho_dot_upper_a_constraint)) > 0) = nan;
rho_dot_lower_a_constraint(abs(imag(rho_dot_lower_a_constraint)) > 0) = nan;

a_upper = rho_dot_upper_a_constraint;
a_lower = rho_dot_lower_a_constraint;

h_1 = cross( R , u_p)';
h_2 = cross( u_p , (RA_dot*u_RA + dec_dot*u_dec) );
h_3 = cross( u_p , R_dot )' + cross( R , (RA_dot*u_RA + dec_dot*u_dec));
h_4 = cross( R , R_dot )';

c_0 = norm(h_1)^2;
c_1 = 2*h_1*h_2;
c_2 = 2*h_1*h_3;
c_3 = 2*h_1*h_4;
c_4 = norm(h_2)^2;
c_5 = 2*h_2'*h_3;
c_6 = 2*h_2'*h_4 + norm(h_3)^2;
c_7 = 2*h_3'*h_4;
c_8 = norm(h_4)^2;

P = c_1*rho.^2 + c_2*rho + c_3;
U = c_4*rho.^4 + c_5*rho.^3 + c_6*rho.^2 + c_7*rho + c_8;

% Eccentricity Constraint
a_0 = F.*U + mu^2.*(1-e^2);
a_1 = F.*P + w_1*U;
a_2 = U + c_0*F + w_1*P;
a_3 = P + c_0*w_1;
a_4 = c_0;

%
% R_ECI = R;
% up = u_p;
% ua = u_RA;
% ud = u_dec;
% 
% h1 = cross(R_ECI, up);
% h2 = cross(up, (RA_dot*ua + dec_dot*ud).');
% h3 = cross(up, R_dot) + cross(R_ECI, (RA_dot*ua + dec_dot*ud).');
% h4 = cross(R_ECI, R_dot);
% 
% c0 = norm(h1)^2;
% c1 = dot(2*h1, h2);
% c2 = dot(2*h1, h3);
% c3 = dot(2*h1, h4);
% c4 = norm(h2)^2;
% c5 = dot(2*h2, h3);
% c6 = dot(2*h2, h4) + norm(h3)^2;
% c7 = dot(2*h3, h4);
% c8 = norm(h4)^2;
% 
% w0 = norm(R_ECI)^2;
% w1 = 2*dot(R_dot,up);
% w2 = (RA_dot^2*cos(dec)^2) + dec_dot^2;
% w3 = (2*RA_dot*dot(R_dot, ua))+(2*dec_dot*dot(R_dot,ud));
% w4 = norm(R_dot)^2;
% w5 = 2*dot(R_ECI,up);
% 
% F_rho = w2.*rho.^2+w3.*rho+w4-(2.*mu./sqrt(rho.^2+w5.*rho+w0));
% P_rho = c1*rho.^2 + c2*rho + c3;
% U_rho = c4*rho.^4 + c5*rho.^3 + c6*rho.^2 + c7*rho + c8;
% 
% a_4 = c0;
% a_3 = P_rho + c0*w1;
% a_2 = U_rho + c0*F_rho + w1*P_rho;
% a_1 = F_rho.*P_rho + w1*U_rho;
% a_0 = F_rho.*U_rho + ((mu^2)*(1-.4^2));

%

% syms rho_dot
sol = zeros(4,length(a_1));
for i = 1:length(a_1)
%     eqn = a_4*rho_dot^4 + a_3(i)*rho_dot^3 + a_2(i)*rho_dot + a_1(i)*rho_dot + a_0 == 0;
%     sol(:,i) = real(  double(  solve(eqn, rho_dot)  )  );
    sol(:,i) = (  (  roots([a_4, a_3(i), a_2(i), a_1(i), a_0(i)])  )  );
end

sol(abs(imag(sol)) > 0) = nan;

ecc_1 = sol(1,:); 
ecc_2 = sol(2,:);  
ecc_3 = sol(3,:); 
ecc_4 = sol(4,:); 

% Sample uniform points =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
x_min = 0;
x_max = 0.86*RE;

% Generate random points
n_points = 200;  % Number of random points to generate
x_random = x_min + (x_max - x_min) * rand(n_points, 1);  % Random x-values between x_min and x_max

E = 0;
F_rand = w_2*x_random.^2 + w_3*x_random + w_4 - (2*mu)./( sqrt(x_random.^2 + w_5*x_random + w_0) );
y_upper =  -w_1/2 + real( sqrt( (w_1/2)^2 - F_rand + 2*E));
y_lower =  -w_1/2 - real( sqrt( (w_1/2)^2 - F_rand + 2*E));

% Generate random y-values between the two functions
y_random = y_lower + (y_upper - y_lower) .* rand(n_points, 1);


r_rand = repmat(R',length(x_random),1) + x_random*u_p';
r_dot_rand = repmat(R_dot,length(x_random),1) + y_random*u_p' + x_random*RA_dot*u_RA' + x_random*dec_dot*u_dec';

options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
j = 1;
for i=1:length(r_rand)
    [~, orbits(:,j:j+5)] = ode45(@two_body_ode,t,[r_rand(i,:)' r_dot_rand(i,:)'],options);
    j = j+6;
end


end