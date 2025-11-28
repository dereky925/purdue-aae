%% Practice Problem 1
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
constraints = norm(A*x-y,2) <= 0.1;

% Solve
options = sdpsettings('verbose',1,'solver','mosek'); 
sol = optimize(constraints, objective, options);

% Check solution status
if sol.problem == 0
    disp('Solution is feasible and optimal:')
    disp(['  Minimum ||x||_2 = ', num2str(value(objective), '%.16f')])
    % disp('  Optimal x:')
    % disp(value(x))
else
    warning(['Something went wrong (Code = ' num2str(sol.problem) ')']);
    disp(sol.info);
end

%% Practice Problem 2
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
constraints = norm(A*x-y,2) <= 0.1;

% Solve
options = sdpsettings('verbose',1,'solver','mosek'); 
sol = optimize(constraints, objective, options);

% Check solution status
if sol.problem == 0
    disp('Solution #1 is feasible and optimal:')
    disp(['  Minimum ||x||_2 = ', num2str(value(objective), '%.16f')])
    % disp('  Optimal x:')
    % disp(value(x))
else
    warning(['Something went wrong (Code = ' num2str(sol.problem) ')']);
    disp(sol.info);
end

%% Practice Problem 3
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
objective   = norm(x,1); 
constraints = norm(A*x - y,2) <= 0.1;

% Solve
options = sdpsettings('verbose',1,'solver','mosek'); 
sol = optimize(constraints, objective, options);

% Check solution status
if sol.problem == 0
    disp('Solution #1 is feasible and optimal:')
    disp(['  Minimum ||x||_2 = ', num2str(value(objective), '%.16f')])
    % disp('  Optimal x:')
    % disp(value(x))
else
    warning(['Something went wrong (Code = ' num2str(sol.problem) ')']);
    disp(sol.info);
end

%% Practice Problem 4 - M spectral norm
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
objective   = norm(A-x*ones(1,n),2); 
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

%% Practice Problem 5 - Frobenius norm 
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
constraints = X >= 0; % Ensure positive semi-definite

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

%% Practice Problem 6 - Constraints
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
objective   = norm(X-0.5*(B+B'),2); 
constraints = [X >= 0, X<=eye(m)]; % Ensure positive semi-definite and lambda max <= 1

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


%% Practice Problem 7: Constrained Problem and Lagrangian Relaxation
clear; clc; close all

% Define matrix A (5 x 10) and vector y (5 x 1)
A = [-5   -11  5  -6   1   0  -6   6  -6   -3
      8    2  -5   0   1  -2   8   -3   7   4
      1   -6   3   4  -3  -4  -1   11  -3  -8
      0   -4   9   0  -7  -3  -10  -5   3   0
      8   2   -5   -5  3  -5  -3   -5  -10  3];

y = [1; 2; 3; 4; 5];

B = A(:,1:5);

yalmip('clear')

% Define the decision variable x.
% Note: n should match the number of columns of A. Here A is 5x10.
n = size(A,2);
x = sdpvar(n,1); 

% Define the objective and constraint for the original constrained problem.
objective   = norm(x,2); 
constraints = norm(A*x - y,2) <= 0.1;

% Solve the constrained problem
options = sdpsettings('verbose',1,'solver','mosek'); 
sol = optimize(constraints, objective, options);

if sol.problem == 0
    disp('Constrained solution is feasible and optimal:')
    disp(['  Minimum ||x||_2 = ', num2str(value(objective), '%.16f')])
else
    warning(['Something went wrong (Code = ' num2str(sol.problem) ')']);
    disp(sol.info);
end

% Retrieve the Lagrange multiplier (dual variable) for the constraint
% In YALMIP, each constraint has an associated dual variable.
lambda = dual(constraints(1));

% Build and solve the corresponding unconstrained Lagrangian problem
% The constraint can be written as: norm(A*x-y,2) - 0.1 <= 0.
% Thus, the Lagrangian is:
%    L(x,lambda) = norm(x,2) + lambda*(norm(A*x-y,2) - 0.1)
% Since we already obtained lambda from the constrained problem,
% we use that value in the Lagrangian objective.

Lagrangian_obj = norm(x,2) + lambda*(norm(A*x-y,2) - 0.1);

% Solve the unconstrained problem with the Lagrangian objective.
% (There are no constraints here)
sol2 = optimize([], Lagrangian_obj, options);

if sol2.problem == 0
    disp('Unconstrained Lagrangian problem solved:')
    disp(['  Minimum Lagrangian value = ', num2str(value(Lagrangian_obj), '%.16f')])
    % Optionally, display the optimal x:
    disp('  Optimal x for the unconstrained problem:')
    disp(value(x))
else
    warning(['Something went wrong in the Lagrangian problem (Code = ' num2str(sol2.problem) ')']);
    disp(sol2.info);
end

disp('The Lagrange multiplier (dual variable) for the constraint is:')
disp(lambda)

