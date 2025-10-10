function xk1 = stateFcn_attEuler(xk, u)
% xk = [phi; theta; psi; p; q; r]
% u  = [ M(3); J(:)(9); dt ]

phi=xk(1); theta=xk(2); psi=xk(3);
w = xk(4:6); 

M  = u(1:3);
J  = reshape(u(4:12),3,3);
dt = u(13);

% W(phi,theta)
cphi=cos(phi); sphi=sin(phi);
cth =cos(theta); sth=sin(theta);
cth = sign(cth)*max(abs(cth),1e-6);   % protezione singolarità

W = [ 1, sphi*sth/cth,  cphi*sth/cth;
      0, cphi,        -sphi;
      0, sphi/cth,     cphi/cth ];

% Dinamica rotazionale
wdot = J \ ( M - cross(w, J*w) );

% Propagazione discreta
eta_next = [phi;theta;psi] + W*w*dt;
w_next   = w + wdot*dt;

% wrap angoli
eta_next = arrayfun(@(a) atan2(sin(a),cos(a)), eta_next);

xk1 = [eta_next; w_next];
end
