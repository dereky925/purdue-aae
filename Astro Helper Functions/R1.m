function [matrix] = R1(angle)

% Computes R1 matrix from [rad] input angle

matrix = [1          0           0;
          0    cos(angle)   sin(angle);
          0    -sin(angle)   cos(angle)];

end