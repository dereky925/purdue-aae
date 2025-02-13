%% P1.a
clear;clc;close all

A = [-1 0;0 1];
B = [1;1];

Qc = ctrb(A,B)

disp("Rank = " + string(rank(Qc)));

%% P1.b
clear;clc;close all

A = [-1 0;0 1];
B = [0;1];

Qc = ctrb(A,B)

disp("Rank = " + string(rank(Qc)));

%% P1.c
clear;clc;close all

A = [1 0;0 1];
B = [1;1];

Qc = ctrb(A,B)

disp("Rank = " + string(rank(Qc)));


%% P2
clear;clc;close all

A = [5 1 -1; -1 3 -1;-2 -2 4];
B = [1 0; 1 1; 0 1];

Qc = ctrb(A,B)

disp("Rank = " + string(rank(Qc)));

eig(A)

syms L

matrix = [A-L*eye(3) B]

rank(double(subs(matrix,L,2)))
rank(double(subs(matrix,L,6)))
rank(double(subs(matrix,L,4)))

%% P3

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
Qc = ctrb(A_ans,B);
eigen = eig(A_ans);
disp("L1 Qc:")
disp(Qc)
disp("Rank = " + string(rank(Qc)))
disp("  ")
disp("Eigenvalues: ")
disp(eigen)

rank([A_ans - 0*eye(6) B])
rank([A_ans - 1.4142i*eye(6) B])
rank([A_ans + 1.4142i*eye(6) B])
rank([A_ans - 1i*eye(6) B])
rank([A_ans + 1i*eye(6) B])

%L3
A_ans = double( subs(A, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 0.99, 1, 0, 0}));
B = double( subs(B, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 0.99, 1, 0, 0}));
Qc = ctrb(A_ans,B);
eigen = eig(A_ans);
disp("L3 Qc:")
disp(Qc)
disp("Rank = " + string(rank(Qc)))
disp("  ")
disp("Eigenvalues: ")
disp(eigen)

%L7
A_ans = double( subs(A, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 0.5, 1, 0, 0}));
B = double( subs(B, {m0, m1, m2, l1, l2, g, theta1_e, theta2_e}, {2, 1, 1, 1, 0.5, 1, 0, 0}));
C = [1 0 0 0 0 0];
Qc = ctrb(A_ans,B);
eigen = eig(A_ans);
disp("L7 Qc:")
disp(Qc)
disp("Rank = " + string(rank(Qc)))
disp("  ")
disp("Eigenvalues: ")
disp(eigen)


%% P4
clear;clc;close all

syms Omg k m

A = [0 1 0 0; Omg^2 - k/m 0 k/m 0; 0 0 0 1; k/m 0 Omg^2 - k/m 0];

B = [0; 0; 0; 1/m];

Qc = [B A*B A^2*B A^3*B]

A = double( subs(A, {Omg, k, m}, {1, 0.5, 3}));
B = double( subs(B, {Omg, k, m}, {0.5, 0.5, 1}));

Qc = ctrb(A,B)

rref(Qc)

disp("Rank = " + string(rank(Qc)));

%% P5
clear;clc
syms Omg k m

A = [0 1 0 0; Omg^2-k/m 0 k/m 0;0 0 0 1; k/m 0 Omg^2-k/m 0];

B = [0;-1/m;0;1/m];

% Qc = [B A*B A^2*B A^3*B]

% rank(Qc)
% 
% rref(Qc')
% Qc'

eig(A)

x = sqrt(-(-m*Omg^2+2*k)/m);

rank([A - Omg*eye(4) B])
rank([A - x*eye(4) B])
rank([A + x*eye(4) B])
rank([A + Omg*eye(4) B])

T = [0  -1/m 1  0;
    -1/m  0  1  0;
     0  1/m  0  1;
    1/m   0  1  1];

rref(T)
rankT = rank(T)

inv(T)*A*T
inv(T)*B

%% P6
clear;clc;

lam1 = 4;
eigen1 = 5;

lam2 = 4;
eigen2 = 5;

b1 = 1;
b2 = 1;
ctr = [lam1 - eigen1 0 b1 ; 0 lam2 - eigen2 b2] 
rank(ctr)
rref(ctr)
% rref(ctr)



%% P7
clear;clc;

A = [0 1;4 0];
B = [1;2];

[w,e] = eig(A);

w = [-2; 1];

w'*A
w'*B

-2*w'

A*w
-2*w

A'*w
-2*w

%% P8
clear;clc

A = [1 -2;2 1];
B = [1 ;1];

eig(A)

w = [1j;1];

w'*A
(1+2j)*w'

w'*B