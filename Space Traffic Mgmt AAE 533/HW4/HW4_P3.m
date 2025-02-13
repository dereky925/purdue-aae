
%% Adding Two Gaussian PDFs
clear;clc;close all

vals = 1:10000;

mean1 = 1000;
var1 = 100;
gauss1 = normpdf(vals, mean1, var1);

mean2 = 7000;
var2 = 1000;
gauss2 = normpdf(vals, mean2, var2);

gauss_sum = gauss1 + gauss2;

a = figure(1);
hold on
grid minor

plot(vals,gauss1,'LineWidth',4)
plot(vals,gauss2,'LineWidth',4)
plot(vals,gauss_sum,':','LineWidth',4)
xlabel("Value",'FontSize',25)
ylabel("Probability",'FontSize',25)
legend("Mean = 1000, Var = 100","Mean = 7000, Var = 1000","PDF Sum",'FontSize',25)

a.Position = [100 100 1000 1000];

%% Adding Two Guassian Random Variables
clear;clc;close all

a = 1:100;
mean1 = 10;
var1 = 10;
gauss1 = normrnd(mean1,var1,100);

mean2 = 70;
var2 = 3;
gauss2 = normrnd(mean2,var2,100);

gauss_var_sum = gauss1 + gauss2;

b = figure(2);
hold on
grid minor
histogram(gauss1)
histogram(gauss2)
histogram(gauss_var_sum)
xlabel("Value",'FontSize',25)
ylabel("Count",'FontSize',25)
legend("Mean = 10, Var = 10","Mean = 70, Var = 3","Random Variable Sum",'FontSize',25)


%% Uniform distributions
clear;clc;close all

% Uniform distribution, number between 0 and 1, 10,000 values
n1 = rand(1,10000);
n2 = rand(1,10000);
n3 = rand(1,10000);
n4 = rand(1,10000);
n5 = rand(1,10000);

n = (n1 + n2 + n3 + n4 + n5)/5;

c = figure(3);
hold on
grid minor
histogram(n1)
histogram(n2)
histogram(n3)
histogram(n4)
histogram(n5)
histogram(n)

xlabel("Random Draw Value",'FontSize',25)
ylabel("Count",'FontSize',25)
legend("n1","n2","n3","n4","n5","n",'FontSize',25)

c.Position = [100 100 1000 1000];
% histogram([n1 n2 n3 n4 n5])


