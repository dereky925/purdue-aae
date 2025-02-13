function [sidereal_angle] = UTC_TO_JD()

% Compute Julian Date
JD = data.TimeMJD + 2400000.5; 

% Compute Julian Centuries, T1
T1 = (JD - 2451545) ./ 36525; 

% Compute JD0, keep only integer portion
JD0 = floor(JD); 

% Compute T0
T0 = (JD0 - 2451545) ./ 36525; 

% Mean sidereal time of Greenwich 0h Earth time (UT)
UT = 24110.54841 + 8640184.812866*T0 + 0.093104*T1.^2 - 0.0000062*T1.^3 + (1.0027279093*(mod(JD,1)-.5).*24).*3600;
UT = UT ./ 3600;

% Compute sidereal time
sidereal_time = mod( UT + data.StationGeocentricLongitudedeg./15 , 24);

% % Compute sidereal angle [deg]
sidereal_angle = sidereal_time * 15;

end