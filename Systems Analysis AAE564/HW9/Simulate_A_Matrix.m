% Define the A matrix
A = [1i, 1; 0, 1i];  % Complex A matrix

% Define initial conditions (you can adjust these as needed)
initial_conditions = [1; 0];  % Example initial state

% Define the time span for the simulation
t_span = [0 60];  % Simulate from time t = 0 to t = 10

% Define the ODE function for the state-space system
odeFunc = @(t, x) A * x;

% Integrate the system using ode45
[t, X] = ode45(odeFunc, t_span, initial_conditions);

% Plot the real and imaginary parts of the solution
figure;
subplot(2, 1, 1);
plot(t, real(X(:, 1)), 'b', t, real(X(:, 2)), 'r');
xlabel('Time');
ylabel('Real Part');
legend('x_1', 'x_2');
title('Real Part of the State Variables');

subplot(2, 1, 2);
plot(t, imag(X(:, 1)), 'b', t, imag(X(:, 2)), 'r');
xlabel('Time');
ylabel('Imaginary Part');
legend('x_1', 'x_2');
title('Imaginary Part of the State Variables');