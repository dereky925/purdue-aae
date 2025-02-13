clear;clc;close all

% Example Poisson Distribution


lambda = 4;  % Rate parameter (mean number of events)

% Generate random numbers from the Poisson distribution
random_poisson = poissrnd(lambda, 1000, 1);

% Compute the mean of the generated data
computed_mean = mean(random_poisson);

a = figure(1);
hold on
grid minor
xlabel("Value",'FontSize',25)
ylabel("Count",'FontSize',25)
histogram(random_poisson)

a.Position = [100 100 1000 1000];


mean = 15/1000*0+76/1000*1+144/1000*2+204/1000*3+203/1000*4+134/1000*5+103/1000*6+67/1000*7+30/1000*8+15/1000*9+7/1000*10+2/1000*11;






