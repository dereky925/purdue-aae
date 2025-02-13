function dxdt = two_body_STM(~, x)
   
    mu = 398600;
    r = x(1:3);
    v = x(4:6);
    Phi = reshape(x(7:end),[6,6]);
    dxdt = zeros(size(x));
    dxdt(1:3) = x(4:6);
    dxdt(4:6) = -mu/norm(r)^3 * r;
    
    F = jacobian([r;v]);

    Phi_dot = F * Phi;
    dxdt(7:end) = Phi_dot(:);

end