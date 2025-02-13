function xdot = four_body_ode(~,x,JD)

mu_Earth = 3.9860044e14;
mu_Sun = 1.327e20;
mu_Moon = 4.9048e12;

au = 1.496e11; % Astronomical Unit

E = 1361; % W/m^2 Solar flux at Earth, 1 AU

c = 3e8; % m/s Speed of light

r_vector = x(1:3,1); % 3x1 vector
vel_vector = x(4:6,1); % 3X1 vector

ISS_DragArea = 1500; % m^2 
% https://www.esa.int/Education/Space_In_Bytes/ATV_a_very_special_delivery_-_Lesson_notes

ISS_Mass = 450000; % kg
% https://en.wikipedia.org/wiki/International_Space_Station

C_d = 0.4; % ISS is mostly made of alluminum
C = (1/4 + 1/9*C_d);

r = norm(r_vector);

r_sun = sun(JD)'*au; % Distance to sun
r_moon = moon(JD)'*au; % Distance to moon

S = (x(1:3)-r_sun)/norm(x(1:3)-r_sun); % Sun unit vector

earth_rate = [0 0 15*3.1415/180/3600]; %rad/s

v_ = vel_vector - cross(earth_rate,r_vector)'; % m/s relative air molecule velocity to satellite

air_density_ISS = 3.8e-12 ; % kg/m³ air density at 400 km 
% https://www.esa.int/Education/Space_In_Bytes/ATV_a_very_special_delivery_-_Lesson_notes

rddot = -mu_Earth/r^3 * r_vector ... % Earth
- mu_Sun*((r_vector - r_sun)/norm(r_vector-r_sun)^3 + r_sun/norm(r_sun)^3) ... % Sun
- mu_Moon*((r_vector - r_moon)/norm(r_vector-r_moon)^3 + r_moon/norm(r_moon)^3)... % Moon 
- ISS_DragArea/ISS_Mass * E / c * au^2/(norm(x(1:3)-r_sun)^2)*C*S... % Solar Radiation Pressure
- C/2*air_density_ISS*ISS_DragArea/ISS_Mass*v_.^2.*v_./norm(v_); % Drag

xdot = [vel_vector;rddot];

end