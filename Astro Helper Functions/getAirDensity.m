function rho = getAirDensity(h_km)

% Inputs:
%   h_km  - Geometric altitude above Earth's surface [km]
%
% Output:
%   rho   - Atmospheric density [kg/m^3]
%
% Reference:
%   Vallado, D.A., "Fundamentals of Astrodynamics and Applications,"
%   (Professor's extended table up to 1000 km).

% Define the piecewise table:
% [ alt_min, alt_max, base_alt(h0), nominal_density(rho0), scale_height(H) ]
% Units: altitude in km, density in kg/m^3, scale height in km
    densityTable = [
       0,   25,    0,       1.225e+00,    7.249
      25,   30,   25,       3.899e-02,    6.349
      30,   40,   30,       1.774e-02,    6.682
      40,   50,   40,       3.972e-03,    7.554
      50,   60,   50,       1.057e-03,    8.382
      60,   70,   60,       3.206e-04,    7.714
      70,   80,   70,       8.770e-05,    6.549
      80,   90,   80,       1.905e-05,    5.799
      90,  100,   90,       3.396e-06,    5.382
     100,  110,  100,       5.297e-07,    5.877
     110,  120,  110,       9.661e-08,    7.263
     120,  130,  120,       2.438e-08,    9.473
     130,  140,  130,       8.484e-09,   12.636
     140,  150,  140,       3.845e-09,   16.149
     150,  180,  150,       2.070e-09,   22.523
     180,  200,  180,       5.464e-10,   29.740
     200,  250,  200,       2.789e-10,   37.105
     250,  300,  250,       7.248e-11,   45.546
     300,  350,  300,       2.418e-11,   53.628
     350,  400,  350,       9.518e-12,   53.298
     400,  450,  400,       3.725e-12,   58.515
     450,  500,  450,       1.585e-12,   60.828
     500,  600,  500,       6.967e-13,   63.822
     600,  700,  600,       1.454e-13,   71.835
     700,  800,  700,       3.614e-14,   88.667
     800,  900,  800,       1.170e-14,  124.640
     900, 1000,  900,       5.245e-15,  181.050
     1000, 999999999, 1000,  3.019e-15,  268.000
    ];

    % Clamp altitude if < 0 or > 1000
    if h_km < 0
        error('Altitude cannot be negative for this model.');
    elseif h_km > 1000
        h_km = 1001;
    end

    % Find the row for the given altitude
    rowIndex = -1;
    for i = 1:size(densityTable,1)
        if h_km >= densityTable(i,1) && h_km < densityTable(i,2)
            rowIndex = i;
            break;
        end
    end

    % If not found (e.g. h_km == 1000 exactly), handle boundary
    if rowIndex == -1
        rowIndex = size(densityTable,1);  % last row
    end

    % Extract table parameters
    h0   = densityTable(rowIndex, 3);  % base altitude
    rho0 = densityTable(rowIndex, 4);  % nominal density
    H    = densityTable(rowIndex, 5);  % scale height

    % Exponential formula
    rho = rho0 * exp(-(h_km - h0)/H);

end