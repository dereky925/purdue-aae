function [GCLat , GCLong ] = getGCFromECEF ( R_ECEF )
% Takes : r_ECEF - Spacecraft position vector in ECEF (km)

    % Returns : GCLat - Geocentric Latitude (deg)
    %           GCLong - Geocentric Longitude (deg)
    %           Grab the components of the ECEF vectors
    
xr = R_ECEF (: ,1);
yr = R_ECEF (: ,2);
zr = R_ECEF (: ,3);
% Find the latitude and longitude using trigonometry :
GCLat = 90 - acosd (zr ./( sqrt(xr .^2 + yr .^2 + zr .^2)));
GCLong = atan2d (yr , xr);
end