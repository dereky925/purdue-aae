%P1/P3
clear;clc;close all

Q = [0 1; 0 1];
rref(Q);
N = null(ans,"rational");
disp("System 2 basis:");disp(N);

Q = [1 1; 1 1];
rref(Q);
N = null(ans,"rational");
disp("System 3 basis:");disp(N);

Q = [-2 1;4 -2];
rref(Q);
N = null(ans,"rational");
disp("System 4 basis:");disp(N);



%% P2
clear;clc;close all

syms Omg k m
A = [0 1 0 0; Omg^2-k/m 0 k/m 0;0 0 0 1; k/m 0 Omg^2-k/m 0];
C = [1 0 0 0];

Q = [C;C*A;C*A^2;C*A^3]

rank(Q)

%% P4
clear;clc;close all

A = [-1 0; 0 1];
eig(A)

A = [1 0; 0 1];
eig(A)

A = [0 1; 4 0];
eig(A)

%% P5
clear;clc;close all

A = [5 -1 -2;1 3 -2; -1 -1 4];
C = [1 1 0; 0 1 1];

Q = [C;C*A;C*A^2]
rank(Q)

eig(A)

syms L

matrix = [5-L -1 -2; 1 3-L -2; -1 -1 4-L; 1 1 0; 0 1 1];

rank(double(subs(matrix,L,6)))
rank(double(subs(matrix,L,4)))
rank(double(subs(matrix,L,2)))

%% P6
clear;clc;close all

lam1 = 4;
lam2 = 5;
eigen1 = 2;
eigen2 = 3;
c1 = 0;
c2 = 0;
m = [lam1 - eigen1 0; 0 lam2 - eigen2;c1 c2];
obs = [c1 c2;lam1 - eigen1 0; 0 lam2 - eigen2]
rank(obs)
rref(m)

%% P7

clear;clc;close all

c = 0;
k = 0;

% Define symbolic variables
syms m0 m1 m2 l1 l2 theta1_e theta2_e g real
% Define the matrix M 
M = [m0 + m1 + m2, -m1*l1*cos(theta1_e), -m2*l2*cos(theta2_e);
     -m1*l1*cos(theta1_e), m1*l1^2, 0;
     -m2*l2*cos(theta2_e), 0, m2*l2^2];

% Invert the matrix
K = [-k 0 0;
      0 -m1*l1*g*cos(theta1_e) 0;
      0 0 -m2*l2*g*cos(theta2_e)];

C = [-c 0 0; 0 0 0; 0 0 0];

Phi = [1; 0; 0];

A_top = [zeros(3), eye(3)];
A_bottom = [M\K, M\C];
A = [A_top; A_bottom];

B = [0; 0; 0; M\Phi];

%L1
A_ans = double( subs(A, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 1, 1, 0, 0}));
B = double( subs(B, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 1, 1, 0, 0}));
C = [1 0 0 0 0 0];
Q0 = [C;C*A_ans;C*A_ans^2;C*A_ans^3;C*A_ans^4;C*A_ans^5];
eigen = eig(A_ans);
disp("L1:")
disp(Q0)
disp("Rank = " + string(rank(Q0)))
disp("  ")
disp("Eigenvalues: ")
disp(eigen)

rank([A_ans - 0*eye(6);C])
rank([A_ans - 1.4142i*eye(6);C])
rank([A_ans + 1.4142i*eye(6);C])
rank([A_ans - 1i*eye(6);C])
rank([A_ans + 1i*eye(6);C])

%L3
A_ans = double( subs(A, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 0.99, 1, 0, 0}));
B = double( subs(B, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 0.99, 1, 0, 0}));
C = [1 0 0 0 0 0];
Q0 = [C;C*A_ans;C*A_ans^2;C*A_ans^3;C*A_ans^4;C*A_ans^5];
eigen = eig(A_ans);
disp("L3:")
disp(Q0)
disp("Rank = " + string(rank(Q0)))
disp("  ")
disp("Eigenvalues: ")
disp(eigen)

%L7
A_ans = double( subs(A, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 0.5, 1, 0, 0}));
B = double( subs(B, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 0.5, 1, 0, 0}));
C = [1 0 0 0 0 0];
Q0 = [C;C*A_ans;C*A_ans^2;C*A_ans^3;C*A_ans^4;C*A_ans^5];
eigen = eig(A_ans);
disp("L7:")
disp(Q0)
disp("Rank = " + string(rank(Q0)))
disp("  ")
disp("Eigenvalues: ")
disp(eigen)



