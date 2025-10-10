% Guidance Law
% Initial conditions
r0 = [x0; y0; z0];
v0 = [vx0; vy0; vz0];

%% ===== General parameters =====
g_moon = 1.62;        % m/s^2

% Initial state (start of PHASE 1 - Braking). Adjust to your case.
x0  = 289448.1899567346;        % initial downrange to target (m)
z0  = 15e3; 
y0  = 0e3; 
% initial altitude (m)
vx0 = -1695;          % horizontal speed toward target (m/s) (negative if x decreases)
vz0 = 0;  
vy0 = 0 ;               % vertical descending speed (m/s, negative downward)


% Breaking phase
x1_t = 15000; y1_t = 0; z1_t = 8000;                                                                                                        
vx1_t = -202; vy1_t = 0; vz1_t = -30;
Tgo1  = max(10, (x1_t - x0) / (0.5*(vx0 + vx1_t)));
rt_breaking = [x1_t; y1_t; z1_t];
vt_breaking = [vx1_t; vy1_t; vz1_t];
tgo_breaking = Tgo1;

% Approach phase
x2_t = 0; z2_t = 100; y2_t = 0;
vx2_t = 0; vy2_t = 0; vz2_t = -8;
rt_approach = [x2_t; y2_t; z2_t];
vt_approach = [vx2_t; vy2_t; vz2_t];
at_approach = [ax2_t; ay2_t; az2_t];
tgo_approach = 150;

% Attitude Guidance
a_cmd = [5.1598; 0; -0.2963];
pitch = atan2(a_cmd(1,:),a_cmd(3,:));
acl_cmd_xz = sqrt(a_cmd(1,:).^2+a_cmd(3,:).^2);
bank = atan2(a_cmd(2,:),acl_cmd_xz);
qy = quaternion([0 pitch(1) 0],'euler','ZYX','frame');
qx = quaternion([0 0 bank(1)],'euler','ZYX','frame');
q0 = compact(qx*qy);
