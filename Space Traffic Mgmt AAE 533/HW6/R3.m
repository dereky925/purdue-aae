function [matrix] = R3(angle)

% Computes R3 matrix from [rad] input angle

matrix = [cos(angle)    sin(angle)   0;
         -sin(angle)    cos(angle)   0;
               0             0       1 ];

end