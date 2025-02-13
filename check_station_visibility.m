function visible_indices = check_station_visibility(orbit_eci, station_eci)
    % INPUTS:
    % orbit_eci   - (Nx3) array of satellite ECI coordinates at each timestep [km]
    % station_eci - (Mx3) array of ground station ECI coordinates [km]
    %
    % OUTPUT:
    % visible_indices - Array of indices where the station is within view
    
    % Constants
    Re = 6378.137; % Earth's mean equatorial radius [km]
    FOV = 10; % Field of view (half-cone angle) in degrees
    
    % Initialize output
    visible_indices = [];
    
    % Iterate through orbit vector
    for i = 1:size(orbit_eci, 1)
        % Compute satellite altitude dynamically
        sat_radius = norm(orbit_eci(i, :)); % Distance from Earth's center [km]
        altitude = sat_radius - Re; % Altitude above Earth's surface [km]

        % Compute the ground footprint radius dynamically
        rho = atan((Re + altitude) * tan(deg2rad(FOV)) / Re); % Angle from nadir
        footprint_radius = Re * rho; % Approximate projected footprint [km]

        % Convert satellite ECI to geodetic latitude & longitude
        [sat_lat, sat_lon] = eci2latlon(orbit_eci(i, :));

        % Check each station
        for j = 1:size(station_eci, 1)
            % Convert station ECI to geodetic latitude & longitude
            [stat_lat, stat_lon] = eci2latlon(station_eci(j, :));

            % Compute great-circle distance between satellite and station
            distance = great_circle_distance(sat_lat, sat_lon, stat_lat, stat_lon);

            % Check if the station is within the satellite's viewing footprint
            if distance <= footprint_radius
                visible_indices = [visible_indices; i];
                break; % Only need to store index once per timestep
            end
        end
    end
end

%% Helper Function: Convert ECI to Geodetic Latitude & Longitude
function [lat, lon] = eci2latlon(eci)
    % Convert ECI coordinates to geodetic latitude and longitude
    % Inputs: eci = [x, y, z] (km)
    
    Re = 6378.137; % Earth’s mean equatorial radius (km)
    
    % Compute longitude
    lon = atan2d(eci(2), eci(1)); % atan2(y, x) gives correct quadrant
    
    % Compute geodetic latitude (approximate)
    r_xy = sqrt(eci(1)^2 + eci(2)^2);
    lat = atan2d(eci(3), r_xy); % Approximate latitude

    % Ensure longitude is in range [-180, 180]
    if lon > 180
        lon = lon - 360;
    end
end

%% Helper Function: Compute Great Circle Distance (Haversine Formula)
function d = great_circle_distance(lat1, lon1, lat2, lon2)
    % Compute great-circle distance between two latitude/longitude points
    % Inputs in degrees, output in km
    
    Re = 6378.137; % Earth’s mean equatorial radius (km)
    
    % Convert degrees to radians
    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    lat2 = deg2rad(lat2);
    lon2 = deg2rad(lon2);
    
    % Haversine formula
    delta_lat = lat2 - lat1;
    delta_lon = lon2 - lon1;
    a = sin(delta_lat / 2)^2 + cos(lat1) * cos(lat2) * sin(delta_lon / 2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));
    
    % Distance in km
    d = Re * c;
end