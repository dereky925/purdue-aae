% =========================================================================
% 
% Filename:       HW8.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE533 - Space Traffic Management
% Professor:      Dr. Carolin Frueh
% Contact:        cfrueh@purdue.edu
% Assignment:     HW 8
% Semester:       Fall 2024
% 
% Description:
% 
%
%
% =========================================================================
clear;clc;close all

mu = 3.9860044e14; %m^3/s^2

mean_A = [153446.180; 41874155.872; 0; 3066.875; -11.374; 0]; %[m] [m/s]
P_A = [6494.080   -376.139   0     0.00160    -0.494     0;
       -376.139    22.560    0    -9.883E-4   0.0286     0;
          0          0     1.205     0          0     -6.071E-5;
        0.0160   -9.883E-4   0    4.437E-8  -1.212E-6    0;
       -0.494      0.0286    0   -1.212E-6  3.762E-5     0;
           0          0   -6.071E-5  0          0     3.39E-9];

mean_B = [153446.679; 41874156.372; 5; 3066.865; -11.364; -1.358E-6]; %[m] [m/s]
P_B = [6494.224  -376.156  -4.492E-5  0.016  -0.494  -5.902E-8;
       -376.156   22.561   2.550E-6    -9.885E-3  0.0286   3.419E-9;
       -4.492E-5  2.55E-6   1.205  -1.18E-10  3.419E-9  -6.072E-5;
       0.016  -9.885E-3  -1.18E-10  4.438E-8  -1.212E-6  -1.448E-13;
       -0.494  0.0286  3.419E-9  -1.212E-6  3.762E-5  4.492E-12;
       -5.902E-8  3.419E-9  -6.072E-5  -1.448E-13  4.492E-12  3.392E-9];

pho_AB = 15; % [m] combined hard-body radius

% Make Covariances Symmetric
P_A = (P_A + P_A')/2;
P_B = (P_B+ P_B')/2;

% Make covariance semi-positive definite
[V, D] = eig(P_A);        % Compute eigenvalues (D) and eigenvectors (V)
D(D < 0) = 0;            % Set negative eigenvalues to zero
P_A = V * D * V';     % Reconstruct the matrix

[V, D] = eig(P_B);        % Compute eigenvalues (D) and eigenvectors (V)
D(D < 0) = 0;            % Set negative eigenvalues to zero
P_B = V * D * V';     % Reconstruct the matrix

num_samples = 500;
num_samples_A = num_samples;
num_samples_B = num_samples;

total_samples = num_samples_A * num_samples_B; 

[a,ecc,incl,RAAN,argp,nu,truelon,arglat,lonper] = ...
                                  ijk2keplerian(mean_A(1:3), mean_A(4:6));

num_orbits = 0.5;
P_coast = num_orbits * (2*pi)/(sqrt(mu/(a^3))); %Coast time
t = 1:1:P_coast;

% Set random seed for reproducibility
rng('default')

% Generate 1000 samples from a 3D multivariate normal distribution
particles_1 = mvnrnd(mean_A', P_A, num_samples_A);
particles_2 = mvnrnd(mean_B', P_B, num_samples_B);

options = odeset('RelTol', 1e-10,'AbsTol',1e-15);

% Propogate all samples
k = 1;
for i = 1:num_samples
    [Tout_A(:,i), object_A(:,k:k+5)] = ode45(@two_body_ode,t,[particles_1(i,1:3), particles_1(i,4:6)],options);
    [Tout_B(:,i), object_B(:,k:k+5)] = ode45(@two_body_ode,t,[particles_2(i,1:3), particles_2(i,4:6)],options);
    k = k + 6;
end


% Initialize number of unique pairs at each time instance
n_unique_pairs = zeros(length(Tout_A),1); 
pos_ind_A = 1;
pos_ind_B = 1;

% If object A particle had collided, increment this to avoid checking that
% particle again
uncollided_particles = 1;

particle_B_collision_list = zeros(1,num_samples);

% i = Time instance
for i = 1:length(Tout_A)

    % j = loop through Object A particles
    for j = uncollided_particles:num_samples
        
        % j = loop through Object B particles, comparing to a single Object
        % A particle. Start loop at j because we only want unique pairs
        for k = j:num_samples

            % Check if distance is less than hardbody total and if this
            % particular particle B has already collided
            if norm( object_A(i,pos_ind_A:pos_ind_A+2) - object_B(i,pos_ind_B:pos_ind_B+2) ) <= 15 && ~sum(pos_ind_B == particle_B_collision_list)
                
                % fprintf("object_A(%d,%d:%d) Collision\n",i,pos_ind_A,pos_ind_A+2);
                % fprintf("object_B(%d,%d:%d) Collision\n",i,pos_ind_B,pos_ind_B+2);
                % disp(object_A(i,pos_ind_A:pos_ind_A+2))
                % disp(object_B(i,pos_ind_B:pos_ind_B+2))

                % fprintf("object_A position %.2f %.2f %.2f\n",object_A(i,pos_ind_A),object_A(i,pos_ind_A+1), object_A(i,pos_ind_A+2));
                % fprintf("object_B position %.2f %.2f %.2f\n",object_B(i,pos_ind_B),object_B(i,pos_ind_B+1), object_B(i,pos_ind_B+2));
                % fprintf("rho = %.2f\n", norm(object_A(i,pos_ind_A:pos_ind_A+2) - object_B(i,pos_ind_B:pos_ind_B+2)));
                % fprintf("\n")

                fprintf("object_A(%d,%d) Collision\n",i,(pos_ind_A+5)/6);
                fprintf("object_B(%d,%d) Collision\n",i,(pos_ind_B+5)/6);
                fprintf("\n")

                % If so, increment unique pairs
                n_unique_pairs(i) = n_unique_pairs(i) + 1;

                % If this particle of object A collided, don't check if it
                % collides again in the future
                uncollided_particles = uncollided_particles + 1;

                % Store particle B collision index
                particle_B_collision_list(j) = pos_ind_B;

                % Break out of loop to check rest of object B particles
                % since this object A particle already has a collision
                break

            end

            % After checking one object B particle, increment to the next
            % particle
            pos_ind_B = pos_ind_B + 6;

        end

        % Reset object B index
        pos_ind_B = 1;

        % After checking one object A particle, increment to the next
        % particle
        pos_ind_A = pos_ind_A + 6;

    end

    % Reset object A index to uncollided_particles
    pos_ind_A = uncollided_particles;

end

PC_cumulative = sum(n_unique_pairs) / total_samples;
disp('Chance of Collision:');
disp(string(PC_cumulative*100) + '%');


% Plot in 3D
figA = figure(1);
hold on
grid minor
scatter3(particles_1(:,1), particles_1(:,2), particles_1(:,3), '+')
scatter3(particles_2(:,1), particles_2(:,2), particles_2(:,3), '+')
figA.Position = [10 10 1000 1000];

% Set axis labels and title
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Scatter Plot Particle Clouds');
view(110, 10);  % Adjusts to a 3D viewing angle with azimuth 45° and elevation 30°
ax = gca;   
ax.FontSize = 16;
axis equal

legend('PDF A','PDF B','Location','eastoutside')


%% 2D Approximation ========================================================
clc
% Combine covariances
P = P_A + P_B;

% Form the encounter plane at first time instance at TCA
r = mean_B(1:3) - mean_A(1:3);
v = mean_B(4:6) - mean_A(4:6);

k_ = v/norm(v);
i_ = r/norm(r);
j_ = -(cross(i_,k_));
U = [i_ j_ k_];

P_pos = P(1:3,1:3); % Only position covariances
P_enc = U'*P_pos*U;

P_2x2 = [P_enc(1,1), P_enc(1,2); P_enc(2,1), P_enc(2,2)];
sig_x = sqrt(P_2x2(1,1));
sig_y = sqrt(P_2x2(2,2));
rho_xy = P_2x2(1,2) / sig_x / sig_y;
rho_xy = 15;

alpha = 0.5 * atan2( 2*P_2x2(1,2) , (P_enc(1,1) - P_enc(2,2)));

T = [cos(alpha), sin(alpha); -sin(alpha) cos(alpha)];

P_diag = T*P_2x2*T';

sig_1 = sqrt(P_diag(1,1));
sig_2 = sqrt(P_diag(2,2));

f = sig_1/sig_2;
sig = sig_1;
R = norm(r);


r_squared = @(theta) (R + rho_xy.*cos(theta)).^2 .*...
    (cos(alpha).^2 + f.^2.*sin(alpha).^2)...
    + rho_xy.^2.*sin(theta).^2.*(sin(alpha.^2) + f.^2.*cos(alpha).^2) ...
    + 2.*rho_xy.*(1-f.^2)*cos(alpha).*sin(alpha)...
    .* sin(theta).*(R+rho_xy.*cos(theta));


func = @(theta) 1/2*pi .* (f.*rho_xy^2 + R.*f.*rho_xy.*cos(theta))/r_squared(theta)...
    .* (1-exp(-r_squared(theta)/2.*sig^2));

q = integral(func,0,2*pi);

disp('2D Chance of Collision:');
disp(string(q*100) + '%');


