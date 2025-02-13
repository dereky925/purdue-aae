clear;clc
% Define symbolic variables and constants
syms y theta1 theta2 dy dtheta1 dtheta2 u real
syms m0 m1 m2 l1 l2 g real

% Define the state vector
x = [y; theta1; theta2; dy; dtheta1; dtheta2];

% Now we assume you have your A and B matrices properly defined in terms of x and u
% For the sake of illustration, I'm using identity matrices as placeholders
% You will replace this with your system-specific equations

% Define the A matrix (6x6 to represent 6 states)
A_L1 = [0 1 0 0 0 0;
        0 0 1 0 0 0;
        0 0 0 1 0 0;
        -1 0 0 0 0 0;
        0 -1 0 0 0 0;
        0 0 -1 0 0 0];

% Define the B matrix (6x1 to match 6 states and 1 input)
B_L1 = [0; 0; 0; 1; 0; 0];

% Define the output matrix C (3x6 to match 3 outputs and 6 states)
C = [1 0 0 0 0 0;    % Output y
     0 1 0 0 0 0;    % Output theta1
     0 0 1 0 0 0];   % Output theta2

% Define the feedthrough matrix D (3x1 to match 3 outputs and 1 input)
D = [0; 0; 0];

% Check matrix sizes to make sure they are correct
disp('Size of A_L1:');
disp(size(A_L1)); % Should be 6x6

disp('Size of B_L1:');
disp(size(B_L1)); % Should be 6x1

disp('Size of C:');
disp(size(C));    % Should be 3x6

disp('Size of D:');
disp(size(D));    % Should be 3x1

% Create state-space model
sys_L1 = ss(A_L1, B_L1, C, D);

% Compute poles and zeros
poles_L1 = pole(sys_L1);
zeros_L1 = tzero(sys_L1);

% Display poles and zeros
disp('Poles for L1:');
disp(poles_L1);
disp('Zeros for L1:');
disp(zeros_L1);