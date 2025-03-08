function dxdt = two_body_STM(~, x)
   
    mu = 3.9860044e14;
    r = x(1:3);
    v = x(4:6);
    Phi = reshape(x(7:end),[6,6]);

    dxdt = zeros(size(x));

    dxdt(1:3) = v;
    dxdt(4:6) = -mu/norm(r)^3 * r;
    
    a = jacobian([r;v]);

    Phi_dot = a * Phi;
    dxdt(7:end) = Phi_dot(:);

end