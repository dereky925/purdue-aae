function a = jacobian(x)

    mu = 3.9860044e14;
    r = x(1:3,1);
    r_mag = norm(r); 
    X = mu * (-1/r_mag^3*eye(3) + 3/r_mag^5 * diag([r(1)^2, r(2)^2, r(3)^2]));
    X(1,2) = 3 * mu * r(1) * r(2) / r_mag^5;
    X(2,1) = X(1,2);
    X(1,3) = 3 * mu * r(1) * r(3) / r_mag^5;
    X(3,1) = X(1,3);
    X(2,3) = 3 * mu * r(2) * r(3) / r_mag^5;
    X(3,2) = X(2,3);
    
    a = [zeros(3), eye(3); X, zeros(3)];
end