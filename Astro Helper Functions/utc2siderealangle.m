function [sidereal_angle] = utc2siderealangle(year,month,day,hour,minute,second)

% Outputs
% sidereal_angle: [degrees]

% If the month is January or February, adjust year and month
if month <= 2
    year = year - 1;
    month = month + 12;
end

% Calculate A and B
A = floor(year / 100);
B = 2 - A + floor(A / 4);

% Calculate the Julian Day
JD = floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + day + B - 1524.5;

% Compute fractional day
frac_minutes = minute + floor(second/60) + mod(second,60)/60;
frac_hours = hour + floor(frac_minutes/60) + mod(frac_minutes,60)/60;
fractional_day = floor(frac_hours/24) + mod(frac_hours,24)/24;
JD = JD + fractional_day;

% Compute Julian Centuries, T1
T1 = (JD - 2451545) ./ 36525; 

% Compute JD0, keep only integer portion
JD0 = floor(JD); 

% Compute T0
T0 = (JD0 - 2451545) ./ 36525; 

% Compute GMST
GMST = 24110.54841 + 8640184.812866*T0 + 0.093104*T1.^2 - 0.0000062*T1.^3 + (1.0027279093*(mod(JD,1)-.5).*24).*3600;
GMST = GMST ./ 3600;

% Compute sidereal time
EARTH_RATE = 15; % degrees/hr
sidereal_time = mod( GMST, 24);

% Compute sidereal angle [deg]
sidereal_angle = sidereal_time * EARTH_RATE;

end