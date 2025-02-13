function [C] = orbit_DCM(Omega, i, theta) 
    C = make_R(Omega, 3) * make_R(i, 1) * make_R(theta, 3);

function R = make_R(theta, axis)
    c_theta = cos(theta);
    s_theta = sin(theta);
    if axis == 3
        R = [c_theta, -s_theta, 0; 
             s_theta,  c_theta, 0; 
             0,        0,       1];
    elseif axis == 2
        R = [c_theta, 0, s_theta; 
             0,       1,       0;
             -s_theta,0, c_theta];
    elseif axis == 1
        R = [1,        0,        0;
             0,  c_theta, -s_theta; 
             0,  s_theta,  c_theta];
    end
end
end