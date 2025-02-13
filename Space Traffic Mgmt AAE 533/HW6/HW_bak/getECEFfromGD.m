function [ R_ECEF ] = getECEFfromGD (GDLat , GDLong , Alt)
% Takes : GDLat - Geodedic Latitude (deg)
%         GDLong - Geodedic Longitude (deg)
%         Alt - Alitude (km)
%         Returns : R_ECEF - Position of the topocenter in ECEF coordinates (km)
%         Relevant constants , Earth radius a, flattening f
a = 6378.137; % km
f = 1./298.257223563;
% Find the eccentricity of the Earth
e = sqrt (2 .* f - f.^2);
% Find the GD radius N
N = a./ sqrt (1 - e .^2.* sind( GDLat ).^2);
% Plug in the values to get our ECEF position
R_ECEF = [(N + Alt).* cosd( GDLat ).* cosd( GDLong )...
    (N + Alt).* cosd( GDLat ).* sind( GDLong )...
    (N.*(1 - e.^2) + Alt).* sind( GDLat )];
end