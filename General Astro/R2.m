function [matrix] = R2(angle)

% Computes R2 matrix from [rad] input angle

matrix = [cos(angle)    0  -sind(angle);
              0         1       0;
          sin(angle)   0   cos(angle) ];

end