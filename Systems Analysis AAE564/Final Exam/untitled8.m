% Fall 2021 Exam

%% P1 a
clear;clc

A = [1 -1 0; 0 1 -1; -1 0 1];
eig(A)
vecnorm(eig(A),2,2)

% Unstable because for discrete systems if norm(eig) > 1, unbounded

%% P1 b
clear;clc

A = [0 1 0 0;0 0 1 0; 0 0 0 1; 1 0 0 0];
eig(A)


% The presence of a positive real eigenvalue (lambda = +1) indicates 
% that solutions grow unbounded over time if perturbed from the equilibrium


% Unstable

%% P2 
clear;clc

A = [0 1 0 0; 0 0 1 0; 0 0 0 1; 2 3 -1 -3];
B = [ 0 0 0 1]';

C = [-1 0 1 0];
syms s

pretty(simplify(C*inv(s*eye(4)-A)*B))

%% P3
clear;clc

A = [-1 -1 -1; 0 0 -1; 0 0 1];

[v,e] = eig(A)


%% P5
clear;clc

A = [-1 0 0;-1 0 0;-1 -1 1];
B = [1 2 2]';

Qc = [B A*B A^2*B]
rank(Qc)

[v,e] = eig(A)

rank([-1 0 0 1; 1 0 0 2;1 -1 1 2])
rank([-2 0 0 1; 1 -1 0 2;1 -1 0 2])




%% Test 2 P4
clear;clc

Qc = [1 -1; -2 2];
rank(Qc)









