% =========================================================================
% 
% Filename:       HW1_P1.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE590 - Applied Control in Astronautics
% Professor:      Dr. Kenshiro Oguri
% Assignment:     HW 1
% Semester:       Spring 2025
% 
% Description: Homework 1
%
% =========================================================================

%% P1 
clear;clc;close all

mu = 398600.4418; % [km^3/s^2]

utc_datetime = datetime(2024, 10, 19, 17, 54, 19, 'TimeZone', 'UTC');

% % ISS pos vel vector 10 17 2024 12:00:00.00
% r_ijk = [656.0758  6596.8028 1473.0872]'; % [km]
% v_ijk = [-4.9631 -0.8127 5.7866]'; % [km/s]

% Given:
a = 8000; %[km]
ecc = 0.15; 
incl = 80; % [°]
raan = 30; % [°]
argp = 20; % [°]
nu = 0; % [°]

[r_ijk,v_ijk] = keplerian2eci(a, ecc, incl, raan, argp, nu);

% Form vector for orbital element propagation
element_vec = [a,ecc,incl,raan,argp,nu];

% Given: 15 orbits
num_orbits = 15;
T = num_orbits * (2*pi)*(sqrt(a^3/mu)); % Orbit period
t = 1:T;

options = odeset('RelTol',1E-12,'AbsTol',1E-12);

% Cartesian propagation
[~, eci_vec] = ode45(@(t,x) two_body_J2_ode(t,x, utc_datetime),...
                t, [r_ijk v_ijk], options);

% Propagate orbital elements
[~, orb_el_out] = ode45(@(t,x) orbital_elements_eom(t,x,0,utc_datetime),...
                t, element_vec, options);

% Convert orbital elements into equinoctial
[mod_eq_0(1),mod_eq_0(2),mod_eq_0(3),...
    mod_eq_0(4),mod_eq_0(5),mod_eq_0(6)] = ...
    keplerian2equinoctial(element_vec(1),element_vec(2),element_vec(3),...
                          element_vec(4),element_vec(5),element_vec(6));

% Propagate Modified equinoctial elements
[~, mod_eq_el_out] = ode45(@(t,x) ...
    modified_equinoctial_elements_eom(t,x,0,utc_datetime),...
                                        t, mod_eq_0, options);

% Compute specific angular momentum vector
h_vec = cross(r_ijk,v_ijk);

% Compute eccentricity vector
e_vec = 1/mu * cross(v_ijk,h_vec) - r_ijk/norm(r_ijk);

% Create initial Milankovitch vector
milankovitch_0 = [h_vec' e_vec' deg2rad(raan+argp+nu)];

% Propagate Milankovitch elements
[~, milankovitch_out] = ode45(@(t,x) ...
    milankovitch_elements_eom(t,x,0,utc_datetime),...
                              t, milankovitch_0, options);


% Convert orbital element propagation into ECI coordinates for plotting
for i = 1:length(orb_el_out)
    % Orbital elements ====================================================
    % Convert Mean anomaly to true anomaly
    [nu_temp,~] = m2nu(deg2rad(orb_el_out(i,6)), ecc);
    orb_el_out(i,6) = mod(rad2deg(nu_temp),360);

    % Convert orbital elements into cartesian
    [orb_el_ECI_pos(i,:),orb_el_ECI_vel(i,:)] = ...
        keplerian2eci(orb_el_out(i,1),orb_el_out(i,2),orb_el_out(i,3)...
        ,orb_el_out(i,4),orb_el_out(i,5),orb_el_out(i,6));

end

% Convert equinoctial element propagation into ECI for plotting
for i= 1:length(mod_eq_el_out)
    % Modified equinoctial elements =======================================
    [mee(i,1), mee(i,2), mee(i,3), mee(i,4), mee(i,5), mee(i,6)] = ...
        equinoctial2keplerian(mod_eq_el_out(i,1),mod_eq_el_out(i,2),...
                              mod_eq_el_out(i,3),mod_eq_el_out(i,4),...
                              mod_eq_el_out(i,5),mod_eq_el_out(i,6));

    % Keep true anomaly between 0 and 360°
    mee(i,6) = mod(mee(i,6),360);

    % Convert orbital elements into cartesian
    [mod_eq_el_ECI_pos(i,:),mod_eq_el_ECI_vel(i,:)] = ...
        keplerian2eci(mee(i,1), mee(i,2),...
                      mee(i,3) ,mee(i,4), ...
                      mee(i,5),mee(i,6));
end

% Convert milankovitch element propagation into ECI for plotting
for i = 1:length(milankovitch_out)
    % Milankovitch elements =======================================
    milankovitch_out(i,7) = mod(milankovitch_out(i,7),2*pi);

    [milankovitch_r(i,:), milankovitch_v(i,:)] = ...
                        milankovitch2eci(milankovitch_out(i,1:3)', ...
                        milankovitch_out(i,4:6)', milankovitch_out(i,7));

    [mlkvtch_orb_el(i,1),mlkvtch_orb_el(i,2),mlkvtch_orb_el(i,3),...
        mlkvtch_orb_el(i,4),mlkvtch_orb_el(i,5),mlkvtch_orb_el(i,6)] = ...
        eci2keplerian(milankovitch_r(i,:),milankovitch_v(i,:));

    % Compute Milankovitch constraint (h x e)
    milankovitch_constraint(i) = dot(milankovitch_out(i,1:3), ...
                                    milankovitch_out(i,4:6));
end

% Set up orbit plotting
orbit_plot_figure = OrbitPlotSetup();

% Plot stuff
plot3(eci_vec(:,1), eci_vec(:,2), eci_vec(:,3), 'r', 'LineWidth', 3);
plot3(orb_el_ECI_pos(:,1), orb_el_ECI_pos(:,2), orb_el_ECI_pos(:,3), 'b--', 'LineWidth', 3);
plot3(mod_eq_el_ECI_pos(:,1), mod_eq_el_ECI_pos(:,2), mod_eq_el_ECI_pos(:,3), 'g:', 'LineWidth', 3);
plot3(milankovitch_r(:,1), milankovitch_r(:,2), milankovitch_r(:,3), 'c-.', 'LineWidth', 3);

legend('','Cartesian','Classical Elements','Modified Equinoctial','Milankovitch');

orbit_plot_figure.Position = [0 0 1000 1000];


% =========================================================================
% Plot orbital elements over a period

% Compute orbital elements from cartesian
for i = 1:length(eci_vec)
    [a(i),ecc(i),incl(i),raan(i),argp(i),nu(i)] = ...
    eci2keplerian(eci_vec(i,1:3), eci_vec(i,4:6));
end

t_hours = t/3600;

% Create cell arrays to index entire variables
cartesian = {a, ecc, incl, raan, argp, real(nu)};
orb_el = {orb_el_out(:,1), orb_el_out(:,2), orb_el_out(:,3),...
          orb_el_out(:,4), orb_el_out(:,5), orb_el_out(:,6)};
mod_eq = {mee(:,1), mee(:,2), mee(:,3),...
          mee(:,4), mee(:,5), mee(:,6)};
milankovitch = {mlkvtch_orb_el(:,1), mlkvtch_orb_el(:,2),...
                mlkvtch_orb_el(:,3), mlkvtch_orb_el(:,4),...
                mlkvtch_orb_el(:,5), mlkvtch_orb_el(:,6)};

% orb_data_an = {a_a, ecc_a, incl_a, raan_a, argp_a, nu_a};
titles = {'a [km]', 'e', 'i [°]', 'RAAN [°]', '\omega [°]', '\nu [°]'};
ylabels = {'', '', '', '', '', ''}; 

orbital_elements = figure('color','white'); 
sgtitle('Orbital Elements over 15 Orbits','FontSize',30)
orbital_elements.Position = [1000 50 1800 1000];

% For legend only
plot(NaN,NaN,'r','LineWidth',2)
plot(NaN,NaN,'b--','LineWidth',1)
plot(NaN,NaN,'g:','LineWidth',2)
plot(NaN,NaN,'c-.','LineWidth',2)

for i = 1:6
    subplot(2,3,i);       % pick the i-th subplot
    hold on; grid minor;  % set up axes
    plot(t_hours, cartesian{i},'r', 'LineWidth', 4);
    plot(t_hours, orb_el{i},'b--', 'LineWidth', 4);
    plot(t_hours, mod_eq{i},'g:', 'LineWidth', 4);
    plot(t_hours, milankovitch{i},'c-.', 'LineWidth', 2);
    
    ax = gca;
    ax.FontSize = 15;

    title(titles{i}, 'FontSize', 20,'FontWeight', 'bold');
    ylabel(ylabels{i}, 'FontSize', 15);
    xlabel('Time [Hours]')
    
end

legend('Cartesian','Classical Elements','Modified Equinoctial','Milankovitch');

mil_constraint_fig = figure('color','white');
hold on 
grid minor
plot(t_hours,milankovitch_constraint,'LineWidth', 2)

title('Milankovitch Constraint vs. Time', 'FontSize', 20,'FontWeight', 'bold');
ylabel('dot(h,e)', 'FontSize', 15);
xlabel('Time [Hours]')
ax = gca;
ax.FontSize = 25;
mil_constraint_fig.Position = [800 500 1800 1000];


%% P1.c.
clear;clc;close all

mu = 398600.4418; % [km^3/s^2]
orbital_elements = figure('Color','white'); 
colors = {'r', 'b', 'g'};
tol_list = [1E-9, 1E-6, 1E-3];
line_size = [4, 3, 2];


for iter = 1:3
    
    % Change options, interested in lowering with each iteration
    options = odeset('RelTol',tol_list(iter),'AbsTol',tol_list(iter));

    utc_datetime = datetime(2024, 10, 19, 17, 54, 19, 'TimeZone', 'UTC');
    
    % % ISS pos vel vector 10 17 2024 12:00:00.00
    % r_ijk = [656.0758  6596.8028 1473.0872]'; % [km]
    % v_ijk = [-4.9631 -0.8127 5.7866]'; % [km/s]
    
    % Given:
    a = 8000; %[km]
    ecc = 0.15; 
    incl = 80; % [°]
    raan = 30; % [°]
    argp = 20; % [°]
    nu = 0; % [°]
    
    [r_ijk,v_ijk] = keplerian2eci(a, ecc, incl, raan, argp, nu);
    
    % Form vector for orbital element propagation
    element_vec = [a,ecc,incl,raan,argp,nu];
    
    % Given: 15 orbits
    num_orbits = 15;
    T = num_orbits * (2*pi)*(sqrt(a^3/mu)); % Orbit period
    t = 1:T;
   
    % Cartesian propagation
    [~, eci_vec] = ode45(@(t,x) two_body_J2_ode(t,x,utc_datetime),...
                    t, [r_ijk v_ijk], options);
    
    % Propagate orbital elements
    [~, orb_el_out] = ode45(@(t,x) orbital_elements_eom(t,x,0,utc_datetime),...
                    t, element_vec, options);
    
    % Convert orbital elements into equinoctial
    [mod_eq_0(1),mod_eq_0(2),mod_eq_0(3),...
        mod_eq_0(4),mod_eq_0(5),mod_eq_0(6)] = ...
        keplerian2equinoctial(element_vec(1),element_vec(2),element_vec(3),...
                              element_vec(4),element_vec(5),element_vec(6));
    
    % Propagate Modified equinoctial elements
    [~, mod_eq_el_out] = ode45(@(t,x) ...
        modified_equinoctial_elements_eom(t,x,0,utc_datetime),...
                                            t, mod_eq_0, options);
    
    % Compute specific angular momentum vector
    h_vec = cross(r_ijk,v_ijk);
    
    % Compute eccentricity vector
    e_vec = 1/mu * cross(v_ijk,h_vec) - r_ijk/norm(r_ijk);
    
    % Create initial Milankovitch vector
    milankovitch_0 = [h_vec' e_vec' deg2rad(raan+argp+nu)];
    
    % Propagate Milankovitch elements
    [~, milankovitch_out] = ode45(@(t,x) ...
        milankovitch_elements_eom(t,x,0,utc_datetime),...
                                  t, milankovitch_0, options);
    
    
    % Convert orbital element propagation into ECI coordinates for plotting
    for i = 1:length(orb_el_out)
    
        % Orbital elements ====================================================
        % Convert Mean anomaly to true anomaly
        [nu_temp,~] = m2nu(deg2rad(orb_el_out(i,6)), ecc);
        orb_el_out(i,6) = mod(rad2deg(nu_temp),360);
    
        % Convert orbital elements into cartesian
        [orb_el_ECI_pos(i,:),orb_el_ECI_vel(i,:)] = ...
            keplerian2eci(orb_el_out(i,1),orb_el_out(i,2),orb_el_out(i,3)...
            ,orb_el_out(i,4),orb_el_out(i,5),orb_el_out(i,6));
    
    end
    
    % Convert equinoctial element propagation into ECI for plotting
    for i= 1:length(mod_eq_el_out)
        % Modified equinoctial elements =======================================
        [mee(i,1), mee(i,2), mee(i,3), mee(i,4), mee(i,5), mee(i,6)] = ...
            equinoctial2keplerian(mod_eq_el_out(i,1),mod_eq_el_out(i,2),...
                                  mod_eq_el_out(i,3),mod_eq_el_out(i,4),...
                                  mod_eq_el_out(i,5),mod_eq_el_out(i,6));
    
        % Keep true anomaly between 0 and 360°
        mee(i,6) = mod(mee(i,6),360);
    
        % Convert orbital elements into cartesian
        [mod_eq_el_ECI_pos(i,:),mod_eq_el_ECI_vel(i,:)] = ...
            keplerian2eci(mee(i,1), mee(i,2),...
                          mee(i,3) ,mee(i,4), ...
                          mee(i,5),mee(i,6));
    end
    
    % Convert milankovitch element propagation into ECI for plotting
    for i = 1:length(milankovitch_out)
        % Milankovitch elements =======================================
    
        milankovitch_out(i,7) = mod(milankovitch_out(i,7),2*pi);
    
        [milankovitch_r(i,:), milankovitch_v(i,:)] = ...
                            milankovitch2eci(milankovitch_out(i,1:3), ...
                            milankovitch_out(i,4:6), milankovitch_out(i,7));
        [mlkvtch_orb_el(i,1),mlkvtch_orb_el(i,2),mlkvtch_orb_el(i,3),...
            mlkvtch_orb_el(i,4),mlkvtch_orb_el(i,5),mlkvtch_orb_el(i,6)] = ...
            eci2keplerian(milankovitch_r(i,:),milankovitch_v(i,:));
    
        % Compute Milankovitch constraint (h x e)
        milankovitch_constraint(i) = dot(milankovitch_out(i,1:3), ...
                                        milankovitch_out(i,4:6));
    
    end
    
    
    % =========================================================================
    % Plot orbital elements over a period
    
    % Compute orbital elements from cartesian
    for i = 1:length(eci_vec)
        [a(i),ecc(i),incl(i),raan(i),argp(i),nu(i)] = ...
        eci2keplerian(eci_vec(i,1:3), eci_vec(i,4:6));
    end
    
    t_hours = t/3600;
    
    % Create cell arrays to index entire variables
    cartesian = {a, ecc, incl, raan, argp, real(nu)};
    orb_el = {orb_el_out(:,1), orb_el_out(:,2), orb_el_out(:,3),...
              orb_el_out(:,4), orb_el_out(:,5), orb_el_out(:,6)};
    mod_eq = {mee(:,1), mee(:,2), mee(:,3),...
              mee(:,4), mee(:,5), mee(:,6)};
    milankovitch = {mlkvtch_orb_el(:,1), mlkvtch_orb_el(:,2),...
                    mlkvtch_orb_el(:,3), mlkvtch_orb_el(:,4),...
                    mlkvtch_orb_el(:,5), mlkvtch_orb_el(:,6)};
    
    % orb_data_an = {a_a, ecc_a, incl_a, raan_a, argp_a, nu_a};
    titles = {'a [km]', 'e', 'i [°]', 'RAAN [°]', '\omega [°]', '\nu [°]'};
    ylabels = {'', '', '', '', '', ''}; 
    
    
    sgtitle('Orbital Elements over 15 Orbits','FontSize',30)
    orbital_elements.Position = [1000 50 1800 1000];
    
    for i = 1:6
        subplot(2,3,i);       % pick the i-th subplot
        hold on; grid minor;  % set up axes

        % For legend only
        plot(NaN,NaN,'r','LineWidth',5)
        plot(NaN,NaN,'b','LineWidth',5)
        plot(NaN,NaN,'g','LineWidth',5)
        plot(NaN,NaN,'k-','LineWidth',2)
        plot(NaN,NaN,'k--','LineWidth',2)
        plot(NaN,NaN,'k:','LineWidth',2)
        plot(NaN,NaN,'k-.','LineWidth',2)

        plot(t_hours, cartesian{i}, 'LineWidth', line_size(iter), 'Color', colors{iter});
        plot(t_hours, orb_el{i},'--', 'LineWidth', line_size(iter), 'Color', colors{iter});
        plot(t_hours, mod_eq{i},':', 'LineWidth', line_size(iter), 'Color', colors{iter});
        plot(t_hours, milankovitch{i},'-.', 'LineWidth', line_size(iter), 'Color', colors{iter});
        
        ax = gca;
        ax.FontSize = 15;
    
        title(titles{i}, 'FontSize', 20,'FontWeight', 'bold');
        ylabel(ylabels{i}, 'FontSize', 15);
        xlabel('Time [Hours]')
        
    end


end

legend('Tol = 1E-9','Tol = 1E-6','Tol = 1E-3','Cartesian ECI',...
    'Classical Orbital Elements','Modified Equinoctial Elements',...
    'Milankovitch Elements')


%% P1.d. Add SRP to acceleration, same as P1.a otherwise
clear;clc;close all

A_M = 5.4E-6; % Given spacecraft ballistic coefficient

mu = 398600.4418; % [km^3/s^2]

utc_datetime = datetime(2024, 10, 19, 17, 54, 19, 'TimeZone', 'UTC');

% Given:
a = 8000; %[km]
ecc = 0.15; 
incl = 80; % [°]
raan = 30; % [°]
argp = 20; % [°]
nu = 0; % [°]

[r_ijk,v_ijk] = keplerian2eci(a, ecc, incl, raan, argp, nu);

% Form vector for orbital element propagation
element_vec = [a,ecc,incl,raan,argp,nu];

% Given: 15 orbits
num_orbits = 15;
T = num_orbits * (2*pi)*(sqrt(a^3/mu)); % Orbit period
t = 1:T;

options = odeset('RelTol',1E-12,'AbsTol',1E-12);

% Cartesian propagation
[~, eci_vec] = ode45(@(t,x) two_body_J2_SRP_ode(t,x,A_M,utc_datetime),...
                t, [r_ijk v_ijk], options);

% Propagate orbital elements
[~, orb_el_out] = ode45(@(t,x) orbital_elements_eom(t,x,A_M,utc_datetime),...
                t, element_vec, options);

% Convert orbital elements into equinoctial
[mod_eq_0(1),mod_eq_0(2),mod_eq_0(3),...
    mod_eq_0(4),mod_eq_0(5),mod_eq_0(6)] = ...
    keplerian2equinoctial(element_vec(1),element_vec(2),element_vec(3),...
                          element_vec(4),element_vec(5),element_vec(6));

% Propagate Modified Equinoctial Elements
[~, mod_eq_el_out] = ode45(@(t,x) ...
    modified_equinoctial_elements_eom(t,x,A_M,utc_datetime),...
                                        t, mod_eq_0, options);

% Compute specific angular momentum vector
h_vec = cross(r_ijk,v_ijk);

% Compute eccentricity vector
e_vec = 1/mu * cross(v_ijk,h_vec) - r_ijk/norm(r_ijk);

% Create initial Milankovitch vector
milankovitch_0 = [h_vec' e_vec' deg2rad(raan+argp+nu)];

% Propagate Milankovitch elements
[~, milankovitch_out] = ode45(@(t,x) ...
    milankovitch_elements_eom(t,x,A_M,utc_datetime),...
                              t, milankovitch_0, options);


% Convert orbital element propagation into ECI coordinates for plotting
for i = 1:length(orb_el_out)

    % Orbital elements ====================================================
    % Convert Mean anomaly to true anomaly
    [nu_temp,~] = m2nu(deg2rad(orb_el_out(i,6)), ecc);
    orb_el_out(i,6) = mod(rad2deg(nu_temp),360);

    % Convert orbital elements into cartesian
    [orb_el_ECI_pos(i,:),orb_el_ECI_vel(i,:)] = ...
        keplerian2eci(orb_el_out(i,1),orb_el_out(i,2),orb_el_out(i,3)...
        ,orb_el_out(i,4),orb_el_out(i,5),orb_el_out(i,6));

end

% Convert equinoctial element propagation into ECI for plotting
for i= 1:length(mod_eq_el_out)
    % Modified equinoctial elements =======================================
    [mee(i,1), mee(i,2), mee(i,3), mee(i,4), mee(i,5), mee(i,6)] = ...
        equinoctial2keplerian(mod_eq_el_out(i,1),mod_eq_el_out(i,2),...
                              mod_eq_el_out(i,3),mod_eq_el_out(i,4),...
                              mod_eq_el_out(i,5),mod_eq_el_out(i,6));

    % Keep true anomaly between 0 and 360°
    mee(i,6) = mod(mee(i,6),360);

    % Convert orbital elements into cartesian
    [mod_eq_el_ECI_pos(i,:),mod_eq_el_ECI_vel(i,:)] = ...
        keplerian2eci(mee(i,1), mee(i,2),...
                      mee(i,3) ,mee(i,4), ...
                      mee(i,5),mee(i,6));
end

% Convert milankovitch element propagation into ECI for plotting
for i = 1:length(milankovitch_out)
    % Milankovitch elements =======================================

    milankovitch_out(i,7) = mod(milankovitch_out(i,7),2*pi);

    [milankovitch_r(i,:), milankovitch_v(i,:)] = ...
                        milankovitch2eci(milankovitch_out(i,1:3), ...
                        milankovitch_out(i,4:6), milankovitch_out(i,7));
    [mlkvtch_orb_el(i,1),mlkvtch_orb_el(i,2),mlkvtch_orb_el(i,3),...
        mlkvtch_orb_el(i,4),mlkvtch_orb_el(i,5),mlkvtch_orb_el(i,6)] = ...
        eci2keplerian(milankovitch_r(i,:),milankovitch_v(i,:));

    % Compute Milankovitch constraint (h x e)
    milankovitch_constraint(i) = dot(milankovitch_out(i,1:3), ...
                                    milankovitch_out(i,4:6));

end

% Set up orbit plotting
orbit_plot_figure = OrbitPlotSetup();

% Plot stuff
plot3(eci_vec(:,1), eci_vec(:,2), eci_vec(:,3), 'r', 'LineWidth', 3);
plot3(orb_el_ECI_pos(:,1), orb_el_ECI_pos(:,2), orb_el_ECI_pos(:,3), 'b--', 'LineWidth', 3);
plot3(mod_eq_el_ECI_pos(:,1), mod_eq_el_ECI_pos(:,2), mod_eq_el_ECI_pos(:,3), 'g:', 'LineWidth', 3);
plot3(milankovitch_r(:,1), milankovitch_r(:,2), milankovitch_r(:,3), 'c-.', 'LineWidth', 3);


orbit_plot_figure.Position = [0 0 1000 1000];


% =========================================================================
% Plot orbital elements over a period

% Compute orbital elements from cartesian
for i = 1:length(eci_vec)
    [a(i),ecc(i),incl(i),raan(i),argp(i),nu(i)] = ...
    eci2keplerian(eci_vec(i,1:3), eci_vec(i,4:6));
end

t_hours = t/3600;

% Create cell arrays to index entire variables
cartesian = {a, ecc, incl, raan, argp, real(nu)};
orb_el = {orb_el_out(:,1), orb_el_out(:,2), orb_el_out(:,3),...
          orb_el_out(:,4), orb_el_out(:,5), orb_el_out(:,6)};
mod_eq = {mee(:,1), mee(:,2), mee(:,3),...
          mee(:,4), mee(:,5), mee(:,6)};
milankovitch = {mlkvtch_orb_el(:,1), mlkvtch_orb_el(:,2),...
                mlkvtch_orb_el(:,3), mlkvtch_orb_el(:,4),...
                mlkvtch_orb_el(:,5), mlkvtch_orb_el(:,6)};

% orb_data_an = {a_a, ecc_a, incl_a, raan_a, argp_a, nu_a};
titles = {'a [km]', 'e', 'i [°]', 'RAAN [°]', '\omega [°]', '\nu [°]'};
ylabels = {'', '', '', '', '', ''}; 

orbital_elements = figure('Color','white'); 
sgtitle('Orbital Elements over 15 Orbits','FontSize',30)
orbital_elements.Position = [1000 50 1800 1000];

for i = 1:6
    subplot(2,3,i);       % pick the i-th subplot
    hold on; grid minor;  % set up axes
    plot(t_hours, cartesian{i},'r', 'LineWidth', 4);
    plot(t_hours, orb_el{i},'b--', 'LineWidth', 4);
    plot(t_hours, mod_eq{i},'g:', 'LineWidth', 4);
    plot(t_hours, milankovitch{i},'c-.', 'LineWidth', 2);
    
    ax = gca;
    ax.FontSize = 15;

    title(titles{i}, 'FontSize', 20,'FontWeight', 'bold');
    ylabel(ylabels{i}, 'FontSize', 15);
    xlabel('Time [Hours]')
    
end

%% P1.2 3rd body in synodic frame
clear;clc;close all

% Given
mu      = 1.2151E-2;         % Mass ratio (dimensionless)
dEM     = 3.8475E5;          % Earth-Moon distance [km]
GMB     = 4.0350E5;          % GM of Earth+Moon [km^3/s^2]

n   = sqrt( GMB / dEM^3 );   % [1/s], angular velocity of the Earth-Moon pair
T0  = 1/n;                   % Fundamental time scale [s]

%       [ x0,   y0,    z0,   vx0,        vy0,         vz0,   tEnd ]
IC_data = [1.2   0       0         0     -1.06110124   0    6.20628;
           0.85  0    0.17546505   0      0.2628980369  0    2.5543991;
           0.05 -0.05    0         4.0    2.6           0    15];

em_fig = figure('Color','white'); 
hold on;
colors = {'r','b','g'};  

for k = 1:size(IC_data,1)
    % Extract ICs and final time (nondim)
    x0  = IC_data(k,1);
    y0  = IC_data(k,2);
    z0  = IC_data(k,3);
    vx0 = IC_data(k,4);
    vy0 = IC_data(k,5);
    vz0 = IC_data(k,6);
    tEnd= IC_data(k,7);

    X0 = [x0; y0; z0; vx0; vy0; vz0];

    % Integrate from t=0 to t=tEnd in nondimensional CR3BP time
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [t,X] = ode45(@(t,X) cr3bp_EOM(t,X,mu), [0, tEnd], X0, options);

    % X(:,1:3) are the nondimensional positions => rescale by dEM
    xDim = X(:,1)*dEM;  % [km]
    yDim = X(:,2)*dEM;  % [km]
    zDim = X(:,3)*dEM;  % [km]

    % Plot the 2D projection in the rotating (synodic) frame
    plot3(xDim, yDim, zDim, 'Color', colors{k}, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('IC-%d',k));
end

EarthPos = [-mu*dEM, 0];
MoonPos  = [(1 - mu)*dEM, 0];

plot(EarthPos(1), EarthPos(2), 'ko', 'MarkerFaceColor','blue',...
     'DisplayName','Earth');
plot(MoonPos(1),  MoonPos(2),  'ko', 'MarkerFaceColor',...
     [0.7, 0.7, 0.7],'DisplayName','Moon');

xlabel('x (km)');
ylabel('y (km)');
title('CR3BP Trajectories in Synodic Frame (Dimensional)');
legend('Location','northwest');
axis equal; 
grid minor;
ax = gca;
ax.FontSize = 20;
em_fig.Position = [0 0 1000 1000];
view(3)

%%
clc; clear; close all;

% Given
mu = 1.2151E-2;         % Mass ratio (dimensionless)
mu_E = 3.986004418e5;   % Earth GM [km^3/s^2]
mu_M = 4.9028e3;        % Moon GM [km^3/s^2]
dEM  = 3.8475e5;        % Earth-Moon distance [km]
muEM = mu_E + mu_M;     % total GM of Earth + Moon [km^3/s^2]
n    = sqrt(muEM / dEM^3);   % [rad/s], mean motion
Tn   = 1 / n;               % [s], fundamental period scale

%  IC_data: each row is [ x0,  y0,  z0,  vx0,  vy0,  vz0,  tEnd_nd ]
IC_data = [1.2   0       0         0     -1.06110124   0    6.20628;
           0.85  0    0.17546505   0      0.2628980369  0    2.5543991;
           0.05 -0.05    0         4.0    2.6           0    15];


plotColors = {'r','b','g'};

% Create Figures 
figECI = figure('Color','white'); 
hold on; grid minor; axis equal;
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title({'Trajectories under Earth+Moon Third-Body Perturbation','(ECI Frame)'},'FontSize',16);
figECI.Position = [0 0 1000 1000];

figSyn = figure('Color','white'); 
hold on; grid minor; axis equal;
xlabel('x_s [km]'); ylabel('y_s [km]'); zlabel('z_s [km]');
title({'Trajectories in the Synodic Frame','(Rotating with the Moon)'},'FontSize',16);
figSyn.Position = [1500 0 1000 1000];


% Loop over all 3 initial conditions
for k = 1:3

    nd_IC   = IC_data(k,1:6);  % [x0, y0, z0, vx0, vy0, vz0]
    tNonDim = IC_data(k,7);    % nondimensional final time

    x0_dim  = nd_IC(1)*dEM;
    y0_dim  = nd_IC(2)*dEM;
    z0_dim  = nd_IC(3)*dEM;
    vx0_dim = nd_IC(4)*dEM / Tn;
    vy0_dim = nd_IC(5)*dEM / Tn;
    vz0_dim = nd_IC(6)*dEM / Tn;

    X0 = [x0_dim; y0_dim; z0_dim; vx0_dim; vy0_dim; vz0_dim];

    % Convert nondim final time to seconds
    tEnd = tNonDim * Tn;

    % Integrate 
    options = odeset('RelTol',1e-12, 'AbsTol',1e-12);
    [tSol, eci_out] = ode45( ...
        @(t,y) eom_twobody_thirdbody(t, y, mu_E, mu_M, dEM, n), ...
        [0, tEnd], X0, options);

    x_eci = eci_out(:,1); 
    y_eci = eci_out(:,2); 
    z_eci = eci_out(:,3);

    % Plot in ECI
    figure(figECI); 
    plot3(x_eci, y_eci, z_eci, 'Color', plotColors{k}, ...
          'LineWidth', 2, 'DisplayName', sprintf('IC-%d',k));

    for i = 1:length(tSol)
        % Rotation angle about z
        theta = -n * tSol(i);

        spacecraft_eci  = [ x_eci(i) y_eci(i) z_eci(i)]';

        rSyn = R3(theta) * spacecraft_eci;

        xSyn(i) = rSyn(1);
        ySyn(i) = rSyn(2);
        zSyn(i) = rSyn(3);
    end

    % Plot in Synodic Frame
    figure(figSyn);
    plot3(xSyn, ySyn, zSyn, 'Color', plotColors{k}, ...
          'LineWidth', 2, 'DisplayName', sprintf('IC-%d',k));
end

% ECI
figure(figECI);
% Plot Earth at (0,0,0)
plot3(0,0,0,'ko','MarkerSize',10,'MarkerFaceColor','b','DisplayName','Earth');

% Plot moon
tFine = linspace(0,max(IC_data(:,7))*Tn,400);
rMoon = dEM*[cos(n*tFine); sin(n*tFine); zeros(1,length(tFine))];
plot3(rMoon(1,:), rMoon(2,:), rMoon(3,:),'k--','LineWidth',1.5, ...
      'DisplayName','Moon orbit');

legend('Location','best');
set(gca,'FontSize',12); 
view(3);

% Synodic
figure(figSyn);
% Earth in CR3BP synodic frame => at (-mu*dEM, 0) if t=0
plot3(-mu*dEM, 0, 0, 'bo', 'MarkerSize',10, 'MarkerFaceColor','b',...
      'DisplayName','Earth');

% Moon in CR3BP synodic frame => at ((1-mu)*dEM, 0) if t=0
plot3((1-mu)*dEM, 0, 0, 'ko', 'MarkerSize',10, 'MarkerFaceColor',[0.7,0.7,0.7],...
      'DisplayName','Moon');

legend('Location','best');
set(gca,'FontSize',12); 
view(3);

