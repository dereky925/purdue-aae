function [elevation , azimuth ] = getEleAzi (lat , dec , hourAngle )
%   Takes : GDLat - Geodedic Latitude (deg)
%           TCDec - Topocentric Declination (deg)
%           hourAngle - Hour angle tau (deg)
%   Returns : elevation - Elevation angle (deg)
%             azimuth - Azimuth angle (deg)

% Solve with some trig from the third equation :
elevation = asind (sind(lat).* sind(dec) ...
    + cosd(lat).* cosd(dec).* cosd( hourAngle ));
% Set some intermediate variables to make solving easier , from 1st and 2nd
% equations

x = (sind(lat).* cosd(dec).* cosd( hourAngle ) - cosd(lat).* sind(dec));
y = cosd(dec).* sind( hourAngle );
% Solve for azimuth with atan2 , and make sure azimuth is greater than 0.
azimuth = atan2d (y, x);
azimuth ( azimuth < 0) = azimuth ( azimuth < 0) + 360;
end