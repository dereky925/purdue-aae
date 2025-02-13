clear;clc;close all


%P2 
A = [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 0 0 0];

[V,D] = eig(A);


%P3
A = [0 1 0;0 0 1;-1 1 1];
[V,D] = eig(A);

%P4
A = [0 1 0;0 0 1;-6 -11 -6]; % Companion Matrix
[V,D] = eig(A);

%P5
A = [0 1;-2 2];

[V,D] = eig(A);





