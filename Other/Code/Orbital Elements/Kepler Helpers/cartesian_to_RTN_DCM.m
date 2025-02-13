function C = cartesian_to_RTN_DCM(i,Omega,omega,nu)
    %C = angle2dcm(Omega, i, nu + omega, "ZXZ");
    C = orbit_DCM(Omega, i, nu + omega);
end