function [LMSTd , GMSTd ] = getSidereal (JD , GD_Long )
%   Takes : JD - Julian Date (days)
%           GDLong - Geodedic Longitude (deg)
%   Returns : LMSTd - Local Mean Sidereal Time (deg)
%             GMSTd - Greenwich Mean Sidereal Time (deg)
% JD0 = floor(JD + 0.5) - 0.5;
% T0 = ( JD0 - 2451545) /36525;
% Get our intermediate values
JD0 = floor(JD + 0.5) - 0.5;
T0 = ( JD0 - 2451545) /36525;
T1 = (JD - 2451545) /36525;
% Get the value of UT , remember the 12 hour offset
UT = (JD - floor(JD + 0.5)) *60*60*24 + 86400/2;
GMSTs = 24110.54841 ...
    + 8640184.812866.* T0 ...
    + 0.093104.* T1 .^2 ...
    - 0.0000062.* T1 .^3 + 1.0027279093.* UT;
% Convert GMST to degrees and apply modulo 360 to wrap around
GMSTd = mod(GMSTs * 360 / 60 / 60 / 24, 360);

% Calculate LMST by adding Longitude (convert to seconds) and apply modulo 360
LMSTs = GMSTs + GD_Long * 86400 / 360;
LMSTd = mod(LMSTs * 360 / 60 / 60 / 24, 360);
end