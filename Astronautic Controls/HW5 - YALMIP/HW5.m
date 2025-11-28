% =========================================================================
% 
% Filename:       HW5.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE590 - Applied Control in Astronautics
% Professor:      Dr. Kenshiro Oguri
% Assignment:     HW 5
% Semester:       Spring 2025
% 
% Description: Homework 5
%
% =========================================================================

%% Problem a.1
clear; clc; close all

% Define matrix A (5 x 9) and vector y (5 x 1)
A = [-5   -11  5  -6   1   0  -6   6  -6   -3
      8    2  -5   0   1  -2   8   -3   7   4
      1   -6   3   4  -3  -4  -1   11  -3  -8
      0   -4   9   0  -7  -3  -10  -5   3   0
      8   2   -5   -5  3  -5  -3   -5  -10  3];

y = [1; 2; 3; 4; 5];

B = A(:,1:5);

yalmip('clear')

% Dimension of x must match number of columns of A
n = length(A);

x = sdpvar(n,1); % Decision variable
objective   = norm(x,2); 
constraints = norm(A*x+B*y,2) <= 0.1;

% Solve
options = sdpsettings('verbose',1,'solver','mosek'); 
sol = optimize(constraints, objective, options);

% Check solution status
if sol.problem == 0
    disp('Solution is feasible and optimal:')
    disp(['  Minimum ||x||_2 = ', num2str(value(objective), '%.6f')])
    disp('  Optimal x:')
    disp(value(x))
else
    warning(['Something went wrong (Code = ' num2str(sol.problem) ')']);
    disp(sol.info);
end

%% Problem a.2
clear; clc; close all

% Define matrix A (5 x 9) and vector y (5 x 1)
A = [-5   -11  5  -6   1   0  -6   6  -6   -3
      8    2  -5   0   1  -2   8   -3   7   4
      1   -6   3   4  -3  -4  -1   11  -3  -8
      0   -4   9   0  -7  -3  -10  -5   3   0
      8   2   -5   -5  3  -5  -3   -5  -10  3];

y = [1; 2; 3; 4; 5];

B = A(:,1:5);

yalmip('clear')

% Dimension of x must match number of columns of A
n = length(A);

x = sdpvar(n,1); % Decision variable
objective   = norm(x,2); 
constraints = norm(A*x-y,2) >= 0.1;

% Solve
options = sdpsettings('verbose',1,'solver','mosek'); 
sol = optimize(constraints, objective, options);

% Check solution status
if sol.problem == 0
    disp('Solution #1 is feasible and optimal:')
    disp(['  Minimum ||x||_2 = ', num2str(value(objective), '%.6f')])
    % disp('  Optimal x:')
    % disp(value(x))
else
    warning(['Something went wrong (Code = ' num2str(sol.problem) ')']);
    disp(sol.info);
end



%% Problem a.3
clear; clc; close all

% Define matrix A (5 x 9) and vector y (5 x 1)
A = [-5   -11  5  -6   1   0  -6   6  -6   -3
      8    2  -5   0   1  -2   8   -3   7   4
      1   -6   3   4  -3  -4  -1   11  -3  -8
      0   -4   9   0  -7  -3  -10  -5   3   0
      8   2   -5   -5  3  -5  -3   -5  -10  3];

y = [1; 2; 3; 4; 5];

B = A(:,1:5);

yalmip('clear')

% Dimension of x must match number of columns of A
n = length(A);

x = sdpvar(n,1); % Decision variable
objective   = norm(x,2)^2; 
constraints = norm(A*x+B*y,2) <= 0.1;

% Solve
options = sdpsettings('verbose',1,'solver','mosek'); 
sol = optimize(constraints, objective, options);

% Check solution status
if sol.problem == 0
    disp('Solution is feasible and optimal:')
    disp(['  Minimum ||x||_2 = ', num2str(value(objective), '%.6f')])
    disp('  Optimal x:')
    disp(value(x))
else
    warning(['Something went wrong (Code = ' num2str(sol.problem) ')']);
    disp(sol.info);
end

%% Problem a.4
clear; clc; close all

% Define matrix A (5 x 9) and vector y (5 x 1)
A = [-5   -11  5  -6   1   0  -6   6  -6   -3
      8    2  -5   0   1  -2   8   -3   7   4
      1   -6   3   4  -3  -4  -1   11  -3  -8
      0   -4   9   0  -7  -3  -10  -5   3   0
      8   2   -5   -5  3  -5  -3   -5  -10  3];

y = [1; 2; 3; 4; 5];

B = A(:,1:5);

yalmip('clear')

% Dimension of x must match number of columns of A
n = length(A);

x = sdpvar(n,1); % Decision variable
objective   = -norm(x,1); 
constraints = norm(A*x-y,2) <= 0.1;

% Solve
options = sdpsettings('verbose',1,'solver','mosek'); 
sol = optimize(constraints, objective, options);

% Check solution status
if sol.problem == 0
    disp('Solution is feasible and optimal:')
    disp(['  Minimum ||x||_1 = ', num2str(value(objective), '%.6f')])
    disp('  Optimal x:')
    disp(value(x))
else
    warning(['Something went wrong (Code = ' num2str(sol.problem) ')']);
    disp(sol.info);
end

%% Problem a.5
clear; clc; close all

% Define matrix A (5 x 9) and vector y (5 x 1)
A = [-5   -11  5  -6   1   0  -6   6  -6   -3
      8    2  -5   0   1  -2   8   -3   7   4
      1   -6   3   4  -3  -4  -1   11  -3  -8
      0   -4   9   0  -7  -3  -10  -5   3   0
      8   2   -5   -5  3  -5  -3   -5  -10  3];

y = [1; 2; 3; 4; 5];

B = A(:,1:5);

yalmip('clear')

% Dimension of x must match number of columns of A
n = length(A);

x = sdpvar(n,1); % Decision variable
z = sdpvar(n,1);
objective   = dot(x,z); 
constraints = (norm(A*x-y,2) + norm(z,1) ) <= 0.1;

% Solve
options = sdpsettings('verbose',1,'solver','mosek'); 
sol = optimize(constraints, objective, options);

% Check solution status
if sol.problem == 0
    disp('Solution is feasible and optimal:')
    disp(['  Minimum  = ', num2str(value(objective), '%.6f')])

    disp('  Optimal x:')
    disp(value(x))

    disp('  Optimal z:')
    disp(value(z))
else
    warning(['Something went wrong (Code = ' num2str(sol.problem) ')']);
    disp(sol.info);
end

%% Problem a.6 - Spectral norm
clear; clc; close all

% Define matrix A (5 x 9) and vector y (5 x 1)
A = [-5   -11  5  -6   1   0  -6   6  -6   -3
      8    2  -5   0   1  -2   8   -3   7   4
      1   -6   3   4  -3  -4  -1   11  -3  -8
      0   -4   9   0  -7  -3  -10  -5   3   0
      8   2   -5   -5  3  -5  -3   -5  -10  3];

y = [1; 2; 3; 4; 5];

B = A(:,1:5);

yalmip('clear')

% Dimension of x must match number of columns of A
n = length(A);
m = height(A);

x = sdpvar(m,1); % Decision variable
objective   = norm(A - (x+y)*ones(1,n),2); 
constraints = [];

% Solve
options = sdpsettings('verbose',1,'solver','mosek'); 
sol = optimize(constraints, objective, options);

% Check solution status
if sol.problem == 0
    disp('Solution is feasible and optimal:')
    disp(['  Minimum  = ', num2str(value(objective), '%.6f')])

    disp('  Optimal x:')
    disp(value(x))

else
    warning(['Something went wrong (Code = ' num2str(sol.problem) ')']);
    disp(sol.info);
end

%% Problem a.7 - Frobenius norm PSD
clear; clc; close all

% Define matrix A (5 x 9) and vector y (5 x 1)
A = [-5   -11  5  -6   1   0  -6   6  -6   -3
      8    2  -5   0   1  -2   8   -3   7   4
      1   -6   3   4  -3  -4  -1   11  -3  -8
      0   -4   9   0  -7  -3  -10  -5   3   0
      8   2   -5   -5  3  -5  -3   -5  -10  3];

y = [1; 2; 3; 4; 5];

B = A(:,1:5);

yalmip('clear')

% Dimension of x must match number of columns of A
n = length(A);
m = height(A);

X = sdpvar(m,m); % Decision variable
objective   = norm(X-B,'fro'); 
constraints = [X >= 0 , X == X']; % Ensure positive semi-definite and symmetric

% Solve
options = sdpsettings('verbose',1,'solver','mosek'); 
sol = optimize(constraints, objective, options);

% Check solution status
if sol.problem == 0
    disp('Solution is feasible and optimal:')
    disp(['  Minimum  = ', num2str(value(objective), '%.6f')])

    disp('  Optimal x:')
    disp(value(X))

else
    warning(['Something went wrong (Code = ' num2str(sol.problem) ')']);
    disp(sol.info);
end

%% Problem a.8 - Frobenius norm NSD
clear; clc; close all

% Define matrix A (5 x 9) and vector y (5 x 1)
A = [-5   -11  5  -6   1   0  -6   6  -6   -3
      8    2  -5   0   1  -2   8   -3   7   4
      1   -6   3   4  -3  -4  -1   11  -3  -8
      0   -4   9   0  -7  -3  -10  -5   3   0
      8   2   -5   -5  3  -5  -3   -5  -10  3];

y = [1; 2; 3; 4; 5];

B = A(:,1:5);

yalmip('clear')

% Dimension of x must match number of columns of A
n = length(A);
m = height(A);

X = sdpvar(m,m,'symmetric'); % Decision variable, force symmetry
objective   = norm(X-B,'fro'); 
constraints = X <= 0; % Ensure negative semi-definite 

% Solve
options = sdpsettings('verbose',1,'solver','mosek'); 
sol = optimize(constraints, objective, options);

% Check solution status
if sol.problem == 0
    disp('Solution is feasible and optimal:')
    disp(['  Minimum  = ', num2str(value(objective), '%.6f')])

    disp('  Optimal x:')
    disp(value(X))

else
    warning(['Something went wrong (Code = ' num2str(sol.problem) ')']);
    disp(sol.info);
end

%% Problem a.9 - 2 norm PSD nad lambda max <= 1
clear; clc; close all

% Define matrix A (5 x 9) and vector y (5 x 1)
A = [-5   -11  5  -6   1   0  -6   6  -6   -3
      8    2  -5   0   1  -2   8   -3   7   4
      1   -6   3   4  -3  -4  -1   11  -3  -8
      0   -4   9   0  -7  -3  -10  -5   3   0
      8   2   -5   -5  3  -5  -3   -5  -10  3];

y = [1; 2; 3; 4; 5];

B = A(:,1:5);

yalmip('clear')

% Dimension of x must match number of columns of A
n = length(A);
m = height(A);

X = sdpvar(m,m,'symmetric'); % Decision variable, force symmetry
objective   = norm(X-B,2); 
constraints = [X >= 0, X<=eye(m)]; % Ensure positive semi-definite & lambda max <= 1 

% Solve
options = sdpsettings('verbose',1,'solver','mosek'); 
sol = optimize(constraints, objective, options);

% Check solution status
if sol.problem == 0
    disp('Solution is feasible and optimal:')
    disp(['  Minimum  = ', num2str(value(objective), '%.6f')])

    disp('  Optimal x:')
    disp(value(X))

else
    warning(['Something went wrong (Code = ' num2str(sol.problem) ')']);
    disp(sol.info);
end


