% =========================================================================
% 
% Filename:       HW2.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE590 - Applied Control in Astronautics
% Professor:      Dr. Kenshiro Oguri
% Assignment:     HW 2
% Semester:       Spring 2025
% 
% Description: Homework 1
%
% =========================================================================

%% Part a.1 Matlab vs. Handcalc compare
clear;clc;close all

syms x y z real

values = [x, y, z];  
numerical_values = [1 2 3]; 

r = sqrt(x^2 + y^2 + z^2);

r_vec = [x y z];

matlab_r = gradient(r,r_vec);
% pretty(simplify(matlab_r))

handcalc_r = [x y z]' / r;

% Evaluate
matlab_r_evaluated = double(subs(matlab_r, values, numerical_values));
fprintf('matlab_r_evaluated = [%g', matlab_r_evaluated(1));
fprintf('    %g', matlab_r_evaluated(2:end));
fprintf(']\n');

handcalc_r_evaluated = double(subs(handcalc_r, values, numerical_values));
fprintf('handcalc_r_evaluated = [%g', handcalc_r_evaluated(1));
fprintf('    %g', handcalc_r_evaluated(2:end));
fprintf(']\n');

%% a.2 Matlab vs. Handcalc compare
clear;clc;close all

syms x y z real

values = [x, y, z];  
numerical_values = [1 2 3]; 

r = sqrt(x^2 + y^2 + z^2);

r_vec = [x y z];

r_hat = r_vec'/r; % dimensions matter, r_vec needs to be column or answers will be different

matlab_r_hat = jacobian(r_hat,r_vec);
% pretty(simplify(matlab_r_hat))

handcalc_r_hat = 1/r*(eye(3) - r_hat*r_hat');

% Evaluate
matlab_r_hat_evaluated = double(subs(matlab_r_hat, values, numerical_values));
fprintf('matlab_r_hat_evaluated = \n');
for i = 1:size(matlab_r_hat_evaluated, 1)
    fprintf('   [');
    fprintf('%10.3f', matlab_r_hat_evaluated(i, :)); 
    fprintf(' ]\n');
end
handcalc_r_hat_evaluated = double(subs(handcalc_r_hat, values, numerical_values));
fprintf('handcalc_r_hat_evaluated = \n');
for i = 1:size(handcalc_r_hat_evaluated, 1)
    fprintf('   [');
    fprintf('%10.3f', handcalc_r_hat_evaluated(i, :)); 
    fprintf(' ]\n');
end



%% Part a.3 Matlab vs. Handcalc compare

clear;clc;close all

syms mu r_0 J2 x y z real

% Define r (magnitude of position vector)
r = sqrt(x^2 + y^2 + z^2);

C20 = -J2;
delta = asin(z/r);
P20 = (1/2) * (3*sin(delta)^2 - 1);
U_J2 = (mu/r) * (r_0/r)^2 * P20 * C20;

% Compute the acceleration as the gradient of U_J2
grad_UJ2 = gradient(U_J2, [x y z]);

% disp('J2 Perturbation Acceleration:')
% pretty(simplify(grad_UJ2));

% Matlab derivation =======================================================

% ---- SUBSTITUTE VALUES ----
values = [x, y, z, mu, r_0, J2];  
numerical_values = [7000, 1000, 2000, 3.986e14, 6378e3, 1.08263e-3]; 

% Hand derivation =========================================================
common_factor = (-3 * mu * r_0^2 * J2) / (2 * r^5);
ax = common_factor * ( (1 - 5 * (z^2 / r^2)) * x );
ay = common_factor * ( (1 - 5 * (z^2 / r^2)) * y );
az = common_factor * ( (3 - 5 * (z^2 / r^2)) * z );
% Combine into a vector
a_J2 = [ax; ay; az];

% Substitute into the gradient
matlab_grad_UJ2_evaluated = subs(grad_UJ2, values, numerical_values);
fprintf('matlab_grad_UJ2_evaluated = [%g', matlab_grad_UJ2_evaluated(1));
fprintf('    %g', matlab_grad_UJ2_evaluated(2:end));
fprintf(']\n');

% Substitute into the gradient
handcalc_grad_UJ2_evaluated = subs(a_J2, values, numerical_values);
fprintf('handcalc_grad_UJ2_evaluated = [%g', handcalc_grad_UJ2_evaluated(1));
fprintf('    %g', handcalc_grad_UJ2_evaluated(2:end));
fprintf(']\n');

%% Part a.4 Matlab vs. Handcalc compare
% Partial of a_J2 term

clear;clc;close all

syms mu r_0 J2 x y z real

% Define r (magnitude of position vector)
r = sqrt(x^2 + y^2 + z^2);

U_J2 = 1/r^5*(1-5*z^2/r^2)*[x y z];

% Compute the acceleration as the gradient of U_J2
jacobian_UJ2 = jacobian(U_J2, [x y z]);

% pretty(simplify(jacobian_UJ2))

% ---- SUBSTITUTE VALUES ----
values = [x, y, z, mu, r_0, J2];  
numerical_values = [7000, 1000, 2000, 3.986e14, 6378e3, 1.08263e-3]; 

A = 1/r^5 - 5*z^2/r^7;
r_vec = [x; y; z];

hand_derivation = r_vec *[ (35*z^2*x/r^9-5*x/r^7) ...
                           (35*z^2*y/r^9-5*y/r^7)...
                           (35*z^3/r^9 - 15*z/r^7)  ] + A*eye(3);

% pretty(simplify(hand_derivation))

% Substitute into the matlab computed jacobian
jacob_UJ2_evaluated = subs(jacobian_UJ2, values, numerical_values);
disp('Matlab Evaluated J2 Perturbation Acceleration:')
disp(vpa(jacob_UJ2_evaluated, 6))  % vpa() for numerical precision

% Substitute into the hand-computed jacobian
derek_UJ2_evaluated = subs(hand_derivation, values, numerical_values);
disp('Hand calculated J2 Perturbation Acceleration:')
disp(vpa(derek_UJ2_evaluated, 6))  % vpa() for numerical precision

%% ========================================================================


% =========================================================================


% =========================================================================

%% Finite difference a.1
clear;clc;close all

syms x y z real

values = [x, y, z];  
numerical_values = [1 2 3]; 

r = sqrt(x^2 + y^2 + z^2);

% pretty(simplify(matlab_r))

handcalc_r = [x y z]' / r;

% Convert to function handle for numerical evaluation
f_numeric = matlabFunction(r, 'Vars', {x, y, z});

% Define a reference point x0
x0 = [1; 1; 1];  

% Compute numerical gradient using finite differences
h_values = 10.^-(0:10);  % Step sizes h = 1, 10^-1, ..., 10^-8

grad_numeric = zeros(3, 1); 

for i = 1:length(h_values)

    h = h_values(i);

    for j = 1:3  % Iterate over x, y, z
        e_j = zeros(3,1);
        e_j(j) = 1;  % Unit vector in the j-th direction
    
        % Compute finite difference approximation for gradient
        grad_numeric(j) = (f_numeric(x0(1) + h*e_j(1), ...
                                     x0(2) + h*e_j(2), ...
                                     x0(3) + h*e_j(3)) ...
                                    - f_numeric(x0(1), x0(2), x0(3))) / h;
    end

    % Compute error norm between finite difference and hand calc
    % 'fro' for matrix norm
    error_norms(i) = norm(grad_numeric - double(subs(handcalc_r, [x, y, z], x0')), 'fro');

end

fig = figure('color','white');
loglog(h_values, error_norms, '-o', 'LineWidth', 2);
xlabel('Step size h');
ylabel('Error norm ||J_{Finite Diff} - J_{Analytical}||');
title('Finite Difference Approximation Error - $\frac{\partial r}{\partial \vec{r}}$', 'Interpreter', 'latex');
grid on;
ax = gca;
ax.FontSize = 20;
fig.Position = [0 0 1200 800];
set(gca, 'XDir', 'reverse')


%% Finite difference a.2
clear;clc;close all

syms x y z real

values = [x, y, z];  
numerical_values = [1 2 3]; 

r = sqrt(x^2 + y^2 + z^2);

r_vec = [x y z];

r_hat = r_vec'/r; % dimensions matter, r_vec needs to be column or answers will be different
% pretty(simplify(matlab_r_hat))

% Convert to function handle for numerical evaluation
f_numeric = matlabFunction(r_hat, 'Vars', {x, y, z});

handcalc_r_hat = 1/r*(eye(3) - r_hat*r_hat');

% Define a reference point x0
x0 = [1; 1; 1];  

% Compute numerical gradient using finite differences
h_values = 10.^-(0:10);  % Step sizes h = 1, 10^-1, ..., 10^-8

for i = 1:length(h_values)

    h = h_values(i);
    J_numerical = zeros(3, 3);

    for j = 1:3  % Iterate over x, y, z
        e_j = zeros(3,1);
        e_j(j) = 1;  % Unit vector in direction j

        % Compute finite difference approximation for column j
        J_numerical(:,j) = (f_numeric(x0(1) + h*e_j(1),... % eval x value
                                      x0(2) + h*e_j(2),... % eval y value
                                      x0(3) + h*e_j(3))... % eval z value
                                      - f_numeric(x0(1), x0(2), x0(3))) / h;
    end

    % Compute error norm between finite difference and hand calc
    % 'fro' for matrix norm
    error_norms(i) = norm(J_numerical - double(subs(handcalc_r_hat, [x, y, z], x0')), 'fro');
end

fig = figure('color','white');
loglog(h_values, error_norms, '-o', 'LineWidth', 2);
xlabel('Step size h');
ylabel('Error norm ||J_{Finite Diff} - J_{Analytical}||');
title('Finite Difference Approximation Error - $\frac{\partial \hat \mathbf{r}}{\partial \vec{r}}$', 'Interpreter', 'latex');
grid on;
ax = gca;
ax.FontSize = 20;
fig.Position = [0 0 1200 800];
set(gca, 'XDir', 'reverse')


%% Finite difference a.3
clear; clc; close all

syms x y z real

mu = 1;
r_0 = 1;
J2 = 1;

% Define r (magnitude of position vector)
r = sqrt(x^2 + y^2 + z^2);

C20 = -J2;
delta = asin(z/r);
P20 = (1/2) * (3*sin(delta)^2 - 1);
U_J2 = (mu/r) * (r_0/r)^2 * P20 * C20;

% Convert to function handle for numerical evaluation
f_numeric = matlabFunction(U_J2, 'Vars', {x, y, z});

% Hand derivation =========================================================
common_factor = (-3 * mu * r_0^2 * J2) / (2 * r^5);
ax = common_factor * ( (1 - 5 * (z^2 / r^2)) * x );
ay = common_factor * ( (1 - 5 * (z^2 / r^2)) * y );
az = common_factor * ( (3 - 5 * (z^2 / r^2)) * z );
% Combine into a vector
hand_derivation = [ax; ay; az];

% Define a reference point x0
x0 = [1; 1; 1];  

% Compute numerical gradient using finite differences
h_values = 10.^-(0:10);  % Step sizes h = 1, 10^-1, ..., 10^-8

grad_numeric = zeros(3, 1); 

for i = 1:length(h_values)

    h = h_values(i);

    for j = 1:3  % Iterate over x, y, z
        e_j = zeros(3,1);
        e_j(j) = 1;  % Unit vector in the j-th direction
    
        % Compute finite difference approximation for gradient
        grad_numeric(j) = (f_numeric(x0(1) + h*e_j(1), ...
                                     x0(2) + h*e_j(2), ...
                                     x0(3) + h*e_j(3)) ...
                                    - f_numeric(x0(1), x0(2), x0(3))) / h;
    end

    % Compute error norm between finite difference and hand calc
    % 'fro' for matrix norm
    error_norms(i) = norm(grad_numeric - double(subs(hand_derivation, [x, y, z], x0')), 'fro');

end

fig = figure('color','white');
loglog(h_values, error_norms, '-o', 'LineWidth', 2);
xlabel('Step size h');
ylabel('Error norm ||J_{Finite Diff} - J_{Analytical}||');
title('Finite Difference Approximation Error - $\frac{\partial U_{J2}}{\partial \vec{r}}$', 'Interpreter', 'latex');
grid on;
ax = gca;
ax.FontSize = 20;
fig.Position = [0 0 1200 800];
set(gca, 'XDir', 'reverse')


%% Finite difference a.4
clear; clc; close all

syms mu r_0 J2 x y z real

% Define r (magnitude of position vector)
r = sqrt(x^2 + y^2 + z^2);

U_J2 = 1/r^5 * (1 - 5*z^2/r^2) * [x; y; z];

% Convert to function handle for numerical evaluation
f_numeric = matlabFunction(U_J2, 'Vars', {x, y, z});

% Hand derivation
A = 1/r^5 - 5*z^2/r^7;
r_vec = [x; y; z];
hand_derivation = r_vec * [ (35*z^2*x/r^9 - 5*x/r^7), ...
    (35*z^2*y/r^9 - 5*y/r^7), (35*z^3/r^9 - 15*z/r^7) ] + A * eye(3);

% Define a reference point x0
x0 = [1; 1; 1];  

% Compute numerical gradient using finite differences
h_values = 10.^-(0:10);  % Step sizes h = 1, 10^-1, ..., 10^-8

for i = 1:length(h_values)

    h = h_values(i);
    J_numerical = zeros(3, 3);

    for j = 1:3  % Iterate over x, y, z
        e_j = zeros(3,1);
        e_j(j) = 1;  % Unit vector in direction j

        % Compute finite difference approximation for column j
        J_numerical(:,j) = (f_numeric(x0(1) + h*e_j(1),... % eval x value
                                      x0(2) + h*e_j(2),... % eval y value
                                      x0(3) + h*e_j(3))... % eval z value
                                      - f_numeric(x0(1), x0(2), x0(3))) / h;
    end

    % Compute error norm between finite difference and hand calc
    % 'fro' for matrix norm
    error_norms(i) = norm(J_numerical - double(subs(hand_derivation, [x, y, z], x0')), 'fro');
end

fig = figure('color','white');
loglog(h_values, error_norms, '-o', 'LineWidth', 2);
xlabel('Step size h');
ylabel('Error norm ||J_{Finite Diff} - J_{Analytical}||');
title('Finite Difference Approximation Error - $\frac{\partial}{\partial \mathbf{r}} \left[ \frac{1}{r^5} \left( 1 - \frac{5z^2}{r^2} \right) \mathbf{r} \right]$', 'Interpreter', 'latex');
grid on;
ax = gca;
ax.FontSize = 20;
fig.Position = [0 0 1200 800];
set(gca, 'XDir', 'reverse')


%% ========================================================================


% =========================================================================


% =========================================================================


%% 2.2.a Derive the full aj2 jacobian (just need the z term)
% Just needed the z term that has the 3 with it
% take previous hand calc:

clear;clc;close all

syms mu r_0 J2 x y z real

% Define r (magnitude of position vector)
r = sqrt(x^2 + y^2 + z^2);

% Compute Z term this time: (3 - 5*z^2/r^2)
U_J2 = 1/r^5 * (3 - 5*z^2/r^2)*[x y z];

% Compute the acceleration as the gradient of U_J2
jacobian_UJ2 = jacobian(U_J2, [x y z]);

% pretty(simplify(jacobian_UJ2))

% ---- SUBSTITUTE VALUES ----
values = [x, y, z, mu, r_0, J2];  
numerical_values = [7000, 1000, 2000, 3.986e14, 6378e3, 1.08263e-3]; 


A = 3/r^5 - 5*z^2/r^7;
r_vec = [x; y; z];

% Make new hand calc with the z 3 accounted for
hand_derivation = r_vec *[ (35*z^2*x/r^9-15*x/r^7) ...
                           (35*z^2*y/r^9-15*y/r^7)...
                           (35*z^3/r^9 - 25*z/r^7)  ] + A*eye(3);

% pretty(simplify(hand_derivation))

% Substitute into the matlab computed jacobian
jacob_UJ2_evaluated = subs(jacobian_UJ2, values, numerical_values);
disp('Matlab Evaluated J2 Perturbation Acceleration:')
disp(vpa(jacob_UJ2_evaluated, 6))  % vpa() for numerical precision

% Substitute into the hand-computed jacobian
derek_UJ2_evaluated = subs(hand_derivation, values, numerical_values);
disp('Hand calculated J2 Perturbation Acceleration:')
disp(vpa(derek_UJ2_evaluated, 6))  % vpa() for numerical precision

%% Verify correctness of full aJ2 jacobian
clear;clc;close all

syms mu J2 r_0 x y z

% mu = MU;
% J2 = 0.00108262668;
% r_0 = EARTH_RADIUS;

r = sqrt(x^2 + y^2 + z^2);
r_vec = [x y z]';

% jacobian aJ2 x y term
A_xy = 1/r^5 - 5*z^2/r^7;
jacob_aJ2_XY = r_vec *[ (35*z^2*x/r^9-5*x/r^7) ...
                           (35*z^2*y/r^9-5*y/r^7)...
                           (35*z^3/r^9 - 15*z/r^7)  ] + A_xy*eye(3);

% jacobian aJ2 z term
A_z = 3/r^5 - 5*z^2/r^7;
jacob_aJ2_Z = r_vec *[ (35*z^2*x/r^9-15*x/r^7) ...
                           (35*z^2*y/r^9-15*y/r^7)...
                           (35*z^3/r^9 - 25*z/r^7)  ] + A_z*eye(3);

% Create combined a_J2 jacobian
% jacobian_aJ2 = -3*mu*J2*r_0^2/2 * (jacob_aJ2_XY...
%                                   + jacob_aJ2_XY...
%                                   + jacob_aJ2_Z);

jacobian_aJ2 = -3*mu*J2*r_0^2/2 * ([jacob_aJ2_XY(1:2,:)
                                   + jacob_aJ2_Z(3,:)]);

% Matlab differentiation
a_J2 = -3*mu*J2*r_0^2/2/r^5*[ (1-5*z^2/r^2)*x  (1-5*z^2/r^2)*y  (3-5*z^2/r^2)*z ];

jacobian_UJ2 = jacobian(a_J2, [x y z]);

% ---- SUBSTITUTE VALUES ----
values = [x, y, z, mu, r_0, J2];  
numerical_values = [7000, 1000, 2000, 3.986e14, 6378e3, 1.08263e-3]; 

jacob_UJ2_evaluated = subs(jacobian_UJ2, values, numerical_values);
disp('Matlab Evaluated J2 Perturbation Acceleration:')
disp(vpa(jacob_UJ2_evaluated, 6))  % vpa() for numerical precision

% Substitute into the hand-computed jacobian
derek_UJ2_evaluated = subs(jacobian_aJ2, values, numerical_values);
disp('Hand calculated J2 Perturbation Acceleration:')
disp(vpa(derek_UJ2_evaluated, 6))  % vpa() for numerical precision


%% 2.2.b Verify correctness of hand derivation of a_2B (two body)

clear;clc;close all

syms mu r_0 J2 x y z real

% Define r (magnitude of position vector)
r = sqrt(x^2 + y^2 + z^2);

r_vec = [x y z]';

% two body dynamics
a_2B = -mu/r^3*r_vec;

% Compute the acceleration as the gradient of U_J2
jacobian_a2B = jacobian(a_2B, [x y z]);

% pretty(simplify(jacobian_UJ2))

% ---- SUBSTITUTE VALUES ----
values = [x, y, z, mu, r_0, J2];  
numerical_values = [7000, 1000, 2000, 3.986e14, 6378e3, 1.08263e-3]; 

% 2 Body accel jacobian hand calc
hand_derivation = -mu*(eye(3)/r^3 + r_vec*(-3/r^5*r_vec'));


% Substitute into the matlab computed jacobian
jacob_a2B_evaluated = subs(jacobian_a2B, values, numerical_values);
disp('Matlab Evaluated J2 Perturbation Acceleration:')
disp(vpa(jacob_a2B_evaluated, 6))  % vpa() for numerical precision

% Substitute into the hand-computed jacobian
derek_UJ2_evaluated = subs(hand_derivation, values, numerical_values);
disp('Hand calculated J2 Perturbation Acceleration:')
disp(vpa(derek_UJ2_evaluated, 6))  % vpa() for numerical precision


%% 2.2.b Part 2 Combine all derivations and verify correctness of A matrix
clear;clc;close all

syms mu r_0 J2 x y z real

% Define r (magnitude of position vector)
r = sqrt(x^2 + y^2 + z^2);
r_vec = [x y z]';

% Compute jacobian of two body dynamics
a_2B = -mu/r^3*r_vec;

% Compute the acceleration as the gradient of U_J2
jacobian_a2B = jacobian(a_2B, [x y z]);

% Compute a_J2 X and Y term
a_J2_XY = 1/r^5*(1-5*z^2/r^2)*[x y z];
jacobian_aJ2_XY = jacobian(a_J2_XY, [x y z]);

% Compute a_J2 Z term
a_J2_Z = 1/r^5 * (3 - 5*z^2/r^2)*[x y z];
jacobian_aJ2_Z = jacobian(a_J2_Z, [x y z]);

% Create combined a_J2 jacobian
jacobian_aJ2 = -3*mu*J2*r_0^2/2 * (jacobian_aJ2_XY...
                                    + jacobian_aJ2_XY...
                                    + jacobian_aJ2_Z);

% Form matlab computed A matrix
A_matlab = [zeros(3) eye(3); jacobian_a2B + jacobian_aJ2 zeros(3)];

% ---- SUBSTITUTE VALUES ----
values = [x, y, z, mu, r_0, J2];  
numerical_values = [7000, 1000, 2000, 3.986e14, 6378e3, 1.08263e-3]; 


% EVALUATE MATLAB A MATRIX TO COMPARE
matlab_A_evaluated = subs(A_matlab, values, numerical_values)


% 2 Body accel jacobian hand calc
jacob_a2B_handcalc = -mu*(eye(3)/r^3 + r_vec*(-3/r^5*r_vec'));

% jacobian aJ2 x y term hand calc
A_xy = 1/r^5 - 5*z^2/r^7;
jacob_aJ2_XY_handcalc = r_vec *[ (35*z^2*x/r^9-5*x/r^7) ...
                           (35*z^2*y/r^9-5*y/r^7)...
                           (35*z^3/r^9 - 15*z/r^7)  ] + A_xy*eye(3);

% jacobian aJ2 z term hand calc
A_z = 3/r^5 - 5*z^2/r^7;
jacob_aJ2_Z_handcalc = r_vec *[ (35*z^2*x/r^9-15*x/r^7) ...
                           (35*z^2*y/r^9-15*y/r^7)...
                           (35*z^3/r^9 - 25*z/r^7)  ] + A_z*eye(3);

% Create combined handcalc a_J2 jacobian
jacobian_aJ2_handcalc = -3*mu*J2*r_0^2/2 * (jacob_aJ2_XY_handcalc...
                                          + jacob_aJ2_XY_handcalc...
                                          + jacob_aJ2_Z_handcalc);

% Form handcalc computed A matrix
A_handcalc = [zeros(3), eye(3); (jacob_a2B_handcalc + jacobian_aJ2_handcalc) , zeros(3)];

% EVALUATE HANDCALC A MATRIX TO COMPARE
handcalc_A_evaluated = subs(A_handcalc, values, numerical_values)

% COMPARE. LOOKS GOOD MAN.
matlab_A_evaluated-handcalc_A_evaluated


%% Integrate STM to final time =============================================
clear; clc; close all

% Given:
a = 8000; %[km]
ecc = 0.15; 
incl = 80; % [°]
raan = 30; % [°]
argp = 20; % [°]
nu = 0; % [°]

[pos_0, vel_0] = keplerian2eci(a, ecc, incl, raan, argp, nu);

% Compute orbital period
T_orbit = 2 * pi * sqrt(a^3 / MU); 
numOrbits = 15;
t_final = numOrbits * T_orbit;

% Time span
tspan = [0, t_final];

% Define Initial STM as Identity Matrix
Phi0 = eye(6);
X0 = [pos_0; vel_0; Phi0(:)]; % Append STM as a column vector

% ODE Solver Options
options = odeset('AbsTol',1e-12,'RelTol',1e-12);

% Integrate the Equations of Motion and STM
[t, X] = ode45(@(t,X) stm_2B_J2(t, X), tspan, X0, options);

% Extract the Nominal State and STM Over Time
% X is a matrix with 42 columns: first 6 for [pos;vel] and then 36 for the STM.
x_nom = X(:, 1:6);         % nominal state (position and velocity) over time
Phi_all = X(:, 7:end);     

% Compute the Perturbed State Over Time Using the STM
% Define a small initial perturbation in the state
dx0 = [1; 0; 0; 0; 10E-3; 0];  % 6x1 perturbation vector

% Preallocate array for predicted perturbed state
x_pred = zeros(size(x_nom)); % will be  N x 6

for k = 1:length(t)
    % Extract the STM at time t(k) and reshape it to 6x6.
    Phi_k = reshape(Phi_all(k,:), 6, 6);

    % Predicted perturbed state = nominal state + STM * (initial perturbation)
    % x_pred(k,:) = (x_nom(k,:)' + Phi_k * dx0)';

    x_pert(k,:) = (Phi_k * dx0)'; % STM the perturbations

    % x_pred(k,:) = Phi_k * [pos_0 ; vel_0];
    
end

pos_nom = x_nom(:, 1:3);      % nominal positions over time
% pos_pred = x_pred(:, 1:3);    % predicted positions from the STM mapping
pos_pert = x_pert(:, 1:3);    % STM propagated perturbations

t_nom = t;
% Perturb inital states and simulate
pos_0 = pos_0 + dx0(1:3);
vel_0 = vel_0 + dx0(4:6);
X0 = [pos_0; vel_0; Phi0(:)]; % Append STM as a column vector
[t_pred, X_pert] = ode45(@(t,X) stm_2B_J2(t, X), tspan, X0, options);

pos_pred = X_pert(:, 1:3);  % propogated orbit with perturbed states


pos_pred_interp = interp1(t_pred, pos_pred, t_nom, 'linear', 'extrap');
nominalVsPerturbed = pos_pred_interp - pos_nom;

orb_fig = OrbitPlotSetup();
hold on; grid on; axis equal;
plot3(pos_nom(:,1), pos_nom(:,2), pos_nom(:,3), 'b-', 'LineWidth', 2);
plot3(pos_pred(:,1), pos_pred(:,2), pos_pred(:,3), 'g--', 'LineWidth', 2);

orb_fig.Position = [0 0 1000 1000];

xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('Orbit Trajectory over 15 Orbits');
legend('','Nominal Orbit', 'Perturbed Orbit (via STM)', 'Location', 'Best');

% Report the Final STM
X_final   = X(end,:)'; % final integrated state (42x1)
Phi_final = reshape(X_final(7:end), [6,6]);  % final STM (6x6 matrix)
disp('Final STM (State Transition Matrix) after 15 orbits:');

% Print variable name
fprintf('Phi_final = \n');

% Loop through rows and print each row nicely formatted
for i = 1:size(Phi_final,1)
    fprintf('   [');  % Start row
    fprintf(' %10.5f', Phi_final(i, :));  % Print each element with 5 digits precision
    fprintf(' ]\n');  % End row
end

hours = t/3600;

diff_fig = figure('Color','white');
hold on; grid minor;
plot(hours,nominalVsPerturbed)
title('\deltax(t) vs. Time')
xlabel('Time [hours]')
ylabel('Perturbed Trajectory - Nominal Trajectory [km]')
diff_fig.Position = [200 0 1500 1000];
ax = gca;
ax.FontSize = 20;

del_fig = figure('Color','white');
hold on; grid minor;
plot(hours,pos_pert)
title('\deltax(t) vs. Time')
xlabel('Time [hours]')
ylabel('\deltax(t) [km]')
del_fig.Position = [500 0 1500 1000];
ax = gca;
ax.FontSize = 20;


perturbation_compare = nominalVsPerturbed - pos_pert;

pertdiff_fig = figure('Color','white');
hold on; grid minor;
plot(hours,perturbation_compare)
title('(Perturbed Trajectory - Nominal Trajectory) - \deltax(t) vs. Time')
xlabel('Time [hours]')
ylabel('(Perturbed Trajectory - Nominal Trajectory) - \deltax(t) [km]')
pertdiff_fig.Position = [800 0 1500 1000];
ax = gca;
ax.FontSize = 20;


