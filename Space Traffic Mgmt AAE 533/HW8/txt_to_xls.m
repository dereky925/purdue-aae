% Load the TLE file
filename = 'GEO.txt'; % Replace with your TLE file path
fileID = fopen(filename, 'r');
data = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);
lines = data{1};

% Preallocate storage for TLEs
numTLEs = length(lines) / 3;
satName = cell(numTLEs, 1);
line1 = cell(numTLEs, 1);
line2 = cell(numTLEs, 1);

% Extract each TLE set (satellite name, line 1, line 2)
for i = 1:numTLEs
    satName{i} = lines{3*i - 2}; % Satellite name
    line1{i} = lines{3*i - 1};    % Line 1 of TLE
    line2{i} = lines{3*i};        % Line 2 of TLE
end

% Create a table
TLETable = table(satName, line1, line2, 'VariableNames', {'SatelliteName', 'Line1', 'Line2'});

% Write the table to Excel
writetable(TLETable, 'GEO_TLEs.xlsx');

