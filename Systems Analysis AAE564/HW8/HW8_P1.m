clear;clc;close all

% Example 108

syms t
e = [1 t; 0 1];

C = [1 0];
B = [0;1];

C*e*B

%% 1.a.
clear;clc

A = [-5 2; -12 5];
% C = [1 1];
B = [1;1];
[T,V] = eig(A);
syms t
% jordan = T*expm(V*t)*inv(T)

expm(A*t)*B

%% 1.b.
clear;clc

A = [-5 2; -12 5];
B = [1;2];
[T,V] = eig(A);
syms t
expm(A*t)*B

%% 1.c.
clear;clc

A = [-5 2; -12 5];
B = [1;3];
[T,V] = eig(A);
syms t
expm(A*t)*B


%% 2
clear;clc

A = [-1 0 0 0; 1 -2 0 0; 1 0 -3 0; 1 1 1 1];
C = [0 0 1 0];
B = [0; 0; 1; 0];
[T,V] = eig(A);
syms t
C*expm(A*t)*B
V

%%
clc

A = [0 1; -2 -3];
B = [0; 1];
C = [3 -1];
D = 0;

pretty(C*inv(eye(2)*s-A)*B+D)

[num, den] = ss2tf(A, B, C, D)

%% 4
clc
A = [0 1; -2 -3];
B = [0; 1];
C = [-1 -3];
D = [1];
syms s
pretty(C*inv(eye(2)*s-A)*B+D)

