clear;clc

A = [1 2 3 4; ...
    2 3 4 5; ...
    3 4 5 6];

rref(A);
% null(ans,"rational")
rank(A)
