clear;clc;


syms s

A = [-s^2+1 -s^2+1; -s^2-1 s^2+1];
B = [s^2+1 ; s^2-1];

simplify(A*B)

%% F21 P5
clear;clc

A = [2 -1 -1; 
     -1 3 -2;
     -1 0  1];
characteristic_poly = poly(A);
vpa_eigenvalues = vpa(eig(sym(A)));
eig(A)

%%
clear;clc
A = [2 -3 0; 2 -3 0; 3 -5 1];
eig(A)

%%
clear;clc;
A = [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 0 0 0];

%%
%Exam problem TF
numerator = [1 3 2];
denominator = [1 2 -3];

[a,b,c,d] = tf2ss(numerator,denominator)