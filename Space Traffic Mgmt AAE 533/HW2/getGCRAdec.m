function [RA , dec] = getGCRAdec ( r_ECI )
%       Takes : r_ECI - Spacecraft position vector in ECI (km)
%     Returns : RA - Geocentric Right Ascension (deg)
%               dec - Geocentric Declination (deg)
%               Split the SC vector into components
x = r_ECI(:, 1);
y = r_ECI(:, 2);
z = r_ECI(:, 3);
% Find the magnitude of the projection on the xy plane .
rbar = sqrt(x.^2 + y.^2);
% Use atan to get the declination from the third equation .
dec = atand (z./ rbar); % deg

dec(rbar == 0) = sign(z(rbar == 0)); % deg
% Use atan2 to get the right ascension from the other two equations .
RA = atan2d (y, x); % deg
% Just add some checks to make sure we get a non -NAN value and within our
% bounds
RA ((x == 0) & (y == 0)) = 0; % deg
RA = wrapTo180 (RA); % deg
end