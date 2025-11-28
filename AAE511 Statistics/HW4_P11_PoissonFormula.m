clear; clc; close all

%% Parameters
lambda_per_hour = 8;       % Arrival rate
mu_1h = lambda_per_hour;   % For 1-hour period
mu_75min = lambda_per_hour * (75 / 60);     % For 75 minutes
mu_2_5h = lambda_per_hour * 2.5;            % For 2.5 hours

%% (a) Probabilities for 1-hour period
% a1: P(X = 7)
p_a1 = poisspdf(7, mu_1h);

% a2: P(X >= 7) = 1 - P(X <= 6)
p_a2 = 1 - poisscdf(6, mu_1h);

% a3: P(X >= 13) = 1 - P(X <= 12)
p_a3 = 1 - poisscdf(12, mu_1h);

%% (b) Expected value and std dev for 75-minute period
expected_b = mu_75min;
stddev_b = sqrt(mu_75min);

%% (c) Probabilities for 2.5-hour period
% c1: P(X >= 29) = 1 - P(X <= 28)
p_c1 = 1 - poisscdf(28, mu_2_5h);

% c2: P(X <= 14)
p_c2 = poisscdf(14, mu_2_5h);

%% Display results
fprintf('--- (a) 1-hour Poisson (mu = %.1f) ---\n', mu_1h);
fprintf('P(X = 7)       = %.3f\n', p_a1);
fprintf('P(X >= 7)      = %.3f\n', p_a2);
fprintf('P(X >= 13)     = %.3f\n', p_a3);

fprintf('\n--- (b) 75-minute Poisson (mu = %.2f) ---\n', mu_75min);
fprintf('Expected value = %.3f\n', expected_b);
fprintf('Std deviation  = %.3f\n', stddev_b);

fprintf('\n--- (c) 2.5-hour Poisson (mu = %.1f) ---\n', mu_2_5h);
fprintf('P(X >= 29)     = %.3f\n', p_c1);
fprintf('P(X <= 14)     = %.3f\n', p_c2);