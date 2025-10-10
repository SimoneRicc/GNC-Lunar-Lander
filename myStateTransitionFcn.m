function xk1 = myStateTransitionFcn(xk, u)

    r = xk(1:3); v = xk(4:6);

    FI = u(1); Ts = u(2); m = u(3);

    gI = [0;0;-1.622];
    aI = (1/max(m,1e-6)) * FI + gI;
    
    r_next = r + Ts*v + 0.5*Ts^2*aI;
    v_next = v + Ts*aI;
    
    xk1 = [r_next; v_next];

end