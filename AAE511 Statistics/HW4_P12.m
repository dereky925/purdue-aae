clear; clc; close all

%% Given
% Mean time between loads = 0.4 years
% This implies a rate of 1 load every 0.4 years => lambda = 1 / 0.4 = 2.5 loads per year
lambda = 1 / 0.4;

% Time interval of interest (in years)
T = 4;

% Poisson mean over 4 years
mu = lambda * T;  % mu = 2.5 * 4 = 10

%% Part (a): Expected number of loads in 4 years
expected_loads = mu;  % For Poisson, expected value = mu

fprintf('(a) Expected number of loads in 4 years: %.0f loads\n', expected_loads);

%% Part (b): Probability that more than 11 loads occur
% We want P(X > 11) = 1 - P(X <= 11)
% Compute P(X <= 11) using the Poisson PMF definition:

P_leq_11 = 0;  % initialize sum for cumulative probability

for k = 0:11
    % Poisson PMF: P(X = k) = (e^-mu * mu^k) / k!
    term = (exp(-mu) * mu^k) / factorial(k);
    P_leq_11 = P_leq_11 + term;
end

P_more_than_11 = 1 - P_leq_11;

fprintf('(b) P(X > 11) = 1 - P(X <= 11) = %.3f\n', P_more_than_11);

%% Part (c): Find t such that P(X = 0) <= 0.2
% For Poisson, P(X = 0) = exp(-mu) = exp(-lambda * t)
% Solve exp(-lambda * t) <= 0.2

% Rearranging: 
% -lambda * t <= ln(0.2)
% t >= -ln(0.2) / lambda

log_target = log(0.2);  % natural log of 0.2
t_required = -log_target / lambda;

fprintf('(c) Minimum time t so that P(0 loads) <= 0.2: %.4f years\n', t_required);