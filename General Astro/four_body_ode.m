function xdot = four_body_ode(t,x,utcTime)

r_vector = x(1:3,1); % 3x1 vector
vel_vector = x(4:6,1); % 3X1 vector
r = norm(r_vector);

utcTime = utcTime + seconds(t);
JD = juliandate(utcTime);

omega = 7.2921159e-5; % Earth's angular velocity in rad/s, sidereal day
RE = 6371*1000; % Earth Radius [m]
mu_Earth = 3.9860044e14; %m^3/s^2
mu_Sun = 1.327e20;
mu_Moon = 4.9048e12;
au = 1.496e11; % Astronomical Unit
E = 1361; % W/m^2 Solar flux at Earth, 1 AU
c = 3e8; % m/s Speed of light

ISS_DragArea = 1500; % m^2 
% https://www.esa.int/Education/Space_In_Bytes/ATV_a_very_special_delivery_-_Lesson_notes

ISS_Mass = 450000; % kg
% https://en.wikipedia.org/wiki/International_Space_Station

C_d = 0.4; % ISS is mostly made of alluminum
C = (1/4 + 1/9*C_d);

r_sun = sun(JD)'*au; % Distance to sun
r_moon = moon(JD)'*au; % Distance to moon

S = (x(1:3)-r_sun)/norm(x(1:3)-r_sun); % Sun unit vector

earth_rate = [0 0 15*3.1415/180/3600]; %rad/s

v_ = vel_vector - cross(earth_rate,r_vector)'; % m/s relative air molecule velocity to satellite

air_density_ISS = 3.8e-12 ; % kg/m³ air density at 400 km 
% https://www.esa.int/Education/Space_In_Bytes/ATV_a_very_special_delivery_-_Lesson_notes

% Compute J2 Effects ======================================================

%https://www.vcalc.com/wiki/eng/J2
J2 = 0.00108262668;

% Rotate satellite position to ECEF frame
t = JD * 24 * 60 * 60; % define time in seconds since reference epoch;

% Calculate the Earth rotation angle
theta = omega * t;

% % Nutation Transform
% N = nutation(datevec(utcTime));
% 
% % Precession Transform
% P = precession(JD);
% 
% % Rotate to J2000
% r_vector = P' * N' * r_vector;

r_vec_ECEF = R3(theta) * r_vector;

r_vec_ECEF_norm = norm(r_vec_ECEF);

s = r_vec_ECEF(1)/r_vec_ECEF_norm;
t = r_vec_ECEF(2)/r_vec_ECEF_norm;
u = r_vec_ECEF(3)/r_vec_ECEF_norm;

u_ECEF = [s; t; u];

u_J2 = 3/2 * [5*s*u^2 - s; 5*t*u^2 - t; 5*u^3 - 3*u];

% 2-Body dynamics in ECEF
two_body_ECEF = -mu_Earth/r_vec_ECEF_norm^2*u_ECEF; 
J = mu_Earth*J2/r_vec_ECEF_norm^2*(RE/r_vec_ECEF_norm)*u_J2; % J2

a_J2_ECEF = two_body_ECEF + J;

% Rotate back to ECI
a_J2_ECI = R3(-theta) * a_J2_ECEF;

% Define Accelerations ====================================================

% Earth Two-Body
% -mu_Earth/r^3 * r_vector 

% Sun Perturbations
a_sun = -mu_Sun*((r_vector - r_sun)...
    /norm(r_vector-r_sun)^3 + r_sun/norm(r_sun)^3);

% Moon Perturbations
a_moon = - mu_Moon*((r_vector - r_moon)...
    /norm(r_vector-r_moon)^3 + r_moon/norm(r_moon)^3);

% Solar Radiation Pressure
SRP = -ISS_DragArea/ISS_Mass * E / c * au^2/(norm(x(1:3)-r_sun)^2)*C*S;

% Drag
drag = -C/2*air_density_ISS*ISS_DragArea/ISS_Mass*v_.^2.*v_./norm(v_);

% Sum accelerations =======================================================

rddot =  a_sun + a_moon + SRP + drag + a_J2_ECI;
xdot = [vel_vector;rddot];

end