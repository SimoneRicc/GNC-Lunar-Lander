function xk1 = TranslationStateFcn(xk, u)
% xk : [rI(3); vI(3)]
% u  : [FI(3); m; dt]  (stessa u va anche alla measFcn, se vuoi)
r  = xk(1:3);
v  = xk(4:6);

FI = u(1:3);
m  = max(u(4), eps);     % avoid division by zero
dt = 1/u(5);

% lunar gravitational acceleration 
gI = [0;0;-1.622];

aI = FI./m + gI;

r_next = r + v*dt + 0.5*aI*dt^2;
v_next = v + aI*dt;

xk1 = [r_next; v_next];
end
