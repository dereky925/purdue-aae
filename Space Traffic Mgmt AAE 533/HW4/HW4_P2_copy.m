clear;clc;close all

mu = 3.9860044e14; %m^3/s^2

% Speed of light [m/s]
c = 3E8;

lat = 33.275825;   %[°]
lon = -111.740806; %[°]
alt = 377;         %[m]

lla = [lat, lon, alt];

offset = 300;

%2024-09-25T12:00:00.000 
utc_datetime = datetime(2024, 9, 25, 0, 00, 000, 'TimeZone', 'UTC');

pos_ISS_ECI = [5117.808565310730 -2877.911425706390 3417.677000888100]' .* 1000; %[m] 
vel_ISS_ECI = [5.04434922610630 3.47434703163623 -4.60296327886324]' .* 1000; %[m/s] 

[a,ecc,incl,RAAN,argp,nu,truelon,arglat,lonper] = ijk2keplerian(pos_ISS_ECI, vel_ISS_ECI);


num_orbits = 1.96;
P_coast = num_orbits * (2*pi)/(sqrt(mu/(a^3))); %Coast time
t = 1:1:P_coast;

options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
[Tout, Z] = ode45(@two_body_ode,t,[pos_ISS_ECI vel_ISS_ECI],options);


earth_rot = 7.272E-5; % Earth rotation [rad/s]
RE = 6371*1000; % Earth Radius [m]
max_rho_guess = RE;
rho_guess = 0:100000:max_rho_guess;


% OBSERVATION 1 =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
obs1time = utc_datetime + seconds(Tout(end-offset));
JD1 = juliandate(obs1time);
pos_eci1 = lla2eci(lla, datevec(obs1time)); %[m]

% Nutation Transform
N = nutation(datevec(obs1time));

% Precession Transform
P = precession(JD1);

% Rotate to J2000
station_vec1 = P' * N' * pos_eci1';

[az1,h1,ra1,dec1] = getAzElRaDec(JD1,station_vec1', Z(end-offset,1:3));

% Compute finite difference
[~,~,ra1_FD_1,dec1_FD_1] = getAzElRaDec(JD1,station_vec1', Z(end-offset-1,1:3));
[~,~,ra1_FD_2,dec1_FD_2] = getAzElRaDec(JD1,station_vec1', Z(end-offset+1,1:3));

% Right Ascension finite difference
ra1_dot_1 = (ra1 - ra1_FD_1);
ra1_dot_2 = (ra1_FD_2 - ra1);
ra1_dot = mean([ra1_dot_1,ra1_dot_2]) * pi/180; % Convert to rad/s

% Declination finite difference
dec1_dot_1 = (dec1 - dec1_FD_1);
dec1_dot_2 = (dec1_FD_2 - dec1);
dec1_dot = mean([dec1_dot_1,dec1_dot_2]) * pi/180; % Convert to rad/s

% Compute velocity of station in ECI
R_dot_ECI = cross([0 0 earth_rot],station_vec1);

u_p_1 = [cosd(ra1)*cosd(dec1) sind(ra1)*cosd(dec1) sind(dec1)]';
u_ra_1 = [-sind(ra1)*cosd(dec1) cosd(ra1)*cosd(dec1) 0]';
u_dec_1 = [-cosd(ra1)*sind(dec1) -sind(ra1)*sind(dec1) cosd(dec1)]';

w0_1 = norm(station_vec1)^2;
w1_1 = 2*(R_dot_ECI*u_p_1);
w2_1 = ra1_dot^2*cosd(dec1)^2 + dec1_dot^2;
w3_1 = 2*ra1_dot*(R_dot_ECI*u_ra_1) + 2*dec1_dot*(R_dot_ECI*u_dec_1);
w4_1 = norm(R_dot_ECI)^2;
w5_1 = 2*(station_vec1'*u_p_1);

F_1 = w2_1*rho_guess.^2 + w3_1*rho_guess + w4_1 - (2*mu)./( sqrt(rho_guess.^2 + w5_1*rho_guess + w0_1) );

% % Another way to compute drho
% radical = sqrt( (w1_1/2)^2 - F_1 );
% valid_rho = real(radical)>0;
% % rho = rho(valid_rho);
% % rho = [rho,fliplr(rho)];
% radical = radical(valid_rho);
% drho = -w1_1/2 + radical;

[ecc_1,ecc_2,ecc_3, ecc_4, a_upper, a_lower, zero_upper, zero_lower, rho] = get_a_ecc_constraints(station_vec1,R_dot_ECI,ra1*pi/180,ra1_dot,dec1*pi/180,dec1_dot,0.5);

E = 0;
rho_dot_1_pos =  -w1_1/2 + real( sqrt( (w1_1/2)^2 - F_1 + 2*E));
rho_dot_1_neg =  -w1_1/2 - real( sqrt( (w1_1/2)^2 - F_1 + 2*E));

% Assume object is in LEO < 600 km
E = -mu/(2*RE+600000);
rho_dot_1_pos_a_cons =  -w1_1/2 + real( sqrt( (w1_1/2)^2 - F_1 + 2*E));
rho_dot_1_neg_a_cons =  -w1_1/2 - real( sqrt( (w1_1/2)^2 - F_1 + 2*E));

% Eccentricity constraint
h_1_1 = cross(station_vec1, u_p_1)';
h_2_1 = cross(u_p_1, (ra1_dot*u_ra_1 + dec1_dot*u_dec_1) );
h_3_1 = cross(u_p_1,R_dot_ECI)' + cross(station_vec1,(ra1_dot*u_ra_1 + dec1_dot*u_dec_1));
h_4_1 = cross(station_vec1,R_dot_ECI);

c_0_1 = norm(h_1_1)^2;
c_1_1 = 2*h_1_1*h_2_1;
c_2_1 = 2*h_1_1*h_3_1;
c_3_1 = 2*h_1_1*h_4_1';
c_4_1 = norm(h_2_1)^2;
c_5_1 = 2*h_2_1'*h_3_1;
c_6_1 = 2*h_2_1'*h_4_1' + norm(h_3_1)^2;
c_7_1 = 2*h_3_1'*h_4_1';
c_8_1 = norm(h_4_1)^2;

P_1 = c_1_1*rho_guess.^2 + c_2_1*rho_guess + c_3_1;
U_1 = c_4_1*rho_guess.^4 + c_5_1*rho_guess.^3 + c_6_1*rho_guess.^2 + c_7_1*rho_guess + c_8_1;

% Eccentricity constraint of <0.4
e_constraint = 0.4;

a_4 = c_0_1;
a_3 = P_1 + c_0_1*w1_1;
a_2 = U_1 + c_0_1*F_1 + w1_1*P_1;
a_1 = F_1*P_1' + w1_1*U_1;
a_0 = F_1*U_1' + mu^2*(1-e_constraint^2);

% syms rho_dot
% sol = zeros(4,length(a_1));
% for i = 1:length(a_1)
%     eqn = a_4*rho_dot^4 + a_3(i)*rho_dot^3 + a_2(i)*rho_dot + a_1(i)*rho_dot + a_0 == 0;
%     sol(:,i) = real(  double(  solve(eqn, rho_dot)  )  )
% end

% Sample uniform points =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
x_min = 0;
x_max = 0.438*RE;

% Generate random points
n_points = 200;  % Number of random points to generate
x_random = x_min + (x_max - x_min) * rand(n_points, 1);  % Random x-values between x_min and x_max

E = 0;
F_1_rand = w2_1*x_random.^2 + w3_1*x_random + w4_1 - (2*mu)./( sqrt(x_random.^2 + w5_1*x_random + w0_1) );
y_upper =  -w1_1/2 + real( sqrt( (w1_1/2)^2 - F_1_rand + 2*E));
y_lower =  -w1_1/2 - real( sqrt( (w1_1/2)^2 - F_1_rand + 2*E));

% Generate random y-values between the two functions
y_random = y_lower + (y_upper - y_lower) .* rand(n_points, 1);


num_orbits = 1;
P_coast = num_orbits * (2*pi)/(sqrt(mu/(a^3))); %Coast time
t = 1:1:P_coast;

r_rand = repmat(station_vec1',length(x_random),1) + x_random*u_p_1';
r_dot_rand = repmat(R_dot_ECI,length(x_random),1) + y_random*u_p_1' + x_random*ra1_dot*u_ra_1' + x_random*dec1_dot*u_dec_1';

options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
j = 1;
for i=1:length(r_rand)
    [Tout, Z_rand(:,j:j+5)] = ode45(@two_body_ode,t,[r_rand(i,:)' r_dot_rand(i,:)'],options);
    j = j+6;
end


% OBSERVATION 2 =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
obs2time = utc_datetime + seconds(Tout(end));
JD2 = juliandate(obs2time);

% Compute station ECI at end of orbit propogation
pos_eci2 = lla2eci(lla, datevec(obs2time));

% Nutation Transform
N = nutation(datevec(obs2time));

% Precession Transform
P = precession(JD2);

% Rotate to J2000
station_vec2 = P' * N' * pos_eci2';

[az2,h2,ra2,dec2] = getAzElRaDec(JD2,station_vec2',Z(end,1:3));

u_p_2 = [cosd(ra2)*cosd(dec2) sind(ra2)*cosd(dec2) sind(dec2)]';
u_ra_2 = [-sind(ra2)*cosd(dec2) cosd(ra2)*cosd(dec2) 0]';
u_dec_2 = [-cosd(ra2)*sind(dec2) -sind(ra2)*sind(dec2) cosd(dec2)]';


% Observation 1 Admissible Region  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
close all
a = figure(1);
hold on
grid minor
plot(rho/RE,zero_upper./1000,'r')
plot(rho/RE,zero_lower./1000,'r')

plot(rho/RE,a_upper./1000,'b')
plot(rho/RE,a_lower./1000,'b')

plot(rho/RE,ecc_1./1000,'k')
plot(rho/RE,ecc_2./1000,'k')
plot(rho/RE,ecc_3./1000,'k')
plot(rho/RE,ecc_4./1000,'k')

scatter(x_random/RE,y_random./1000)

% plot(drho)
xlabel('Range [# of Earth Radii]','FontSize',25)
ylabel('Range rate [km/s]','FontSize',25)

a.Position = [100 100 1000 1000];

% Orbit Plot  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
b = figure(3);
hold on
imData = imread('2_no_clouds_4k.jpg');  % Load the texture image

[xS, yS, zS] = sphere(50);  % Generate a sphere for the globe
earth_radius = 6378.1370 * 1000;  % Earth radius in meters

xSE = earth_radius * xS;
ySE = earth_radius * yS;
zSE = earth_radius * zS;


% Create the rotation matrices
Ry = R3(-lon+.6);

Rx = R1(0);

% Combine the rotations (Ry first for longitude, then Rx for latitude)
rotationMatrix = Rx * Ry;

% Apply the rotation to the coordinates
rotated_coords = rotationMatrix * [xSE(:)'; ySE(:)'; zSE(:)'];

% Reshape the rotated coordinates back to match the surface grid dimensions
xSE_rot = reshape(rotated_coords(1, :), size(xSE));
ySE_rot = reshape(rotated_coords(2, :), size(ySE));
zSE_rot = reshape(rotated_coords(3, :), size(zSE));

% Plot the Earth with texture mapping, using the rotated coordinates
surface(xSE_rot, ySE_rot, zSE_rot, 'FaceColor', 'texturemap', 'CData', flipud(imData), 'EdgeColor', 'none');

% Adjust axes
axis equal;
view(3);  % 3D view
grid on;
xlabel('Inertial x (m)');
ylabel('Inertial y (m)');
zlabel('Inertial z (m)');

plot3(Z(:,1), Z(:,2), Z(:,3), 'r', 'LineWidth', 4);

j = 1;
for i=1:length(r_rand)
    plot3(Z_rand(:,j), Z_rand(:,j+1), Z_rand(:,j+2));
    j = j+6;
end

% Plot observer location in ECI
scatter3(pos_eci1(1), pos_eci1(2), pos_eci1(3), 1000, 'c.');
scatter3(pos_eci2(1), pos_eci2(2), pos_eci2(3), 1000, 'm.');

% Mark the end point of the orbit
plot3(Z(end-offset,1), Z(end-offset,2), Z(end-offset,3), 'pm', 'MarkerSize', 20,'MarkerFaceColor','r','MarkerEdgeColor','k');
plot3(Z(end,1), Z(end,2), Z(end,3), 'pm', 'MarkerSize', 20,'MarkerFaceColor','r','MarkerEdgeColor','k');

% Get axis children (just as in the original code)
ch = get(gca, 'children');

b.Position = [1600 100 1000 1000];



















