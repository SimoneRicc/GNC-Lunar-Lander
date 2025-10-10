%% Lunar Powered Descent Guidance (3D) — braking, approach, vertical (jerk-safe)
% Jerk limit per axis in P1/P2, jerk + directional slew + tilt clamp only in P3.
% Implements polynomial guidance laws from Sostaric & Rea with online shaping.
% Fixed: preserves lateral accelerations ax, ay in P1/P2.

clear; clc; rng('default');

%% ===== General parameters =====
g_moon = 1.62;      % m/s^2
dt      = 0.2;      % integration step (s)

% Initial state (start of Phase 1 - Braking). Adjust to your case.
x  = 200e3;   % m (downrange)
y  = 0;                   % m (lateral)
z  = 15e3;                % m (altitude)
vx = -1695;               % m/s
vy = 0;                   % m/s
vz = 0;                   % m/s

%% ===== Optional hard acceleration limits (per-axis) =====
amax_x = 5;               % m/s^2  (set [] to disable)
amax_y = 5;               % m/s^2
amax_z = 5;               % m/s^2

%% ===== Acceleration shaping constraints =====
% Jerk limits (|da/dt| <= Jmax). With dt=0.2s, |Δa| <= Jmax*dt.
Jmax_x = 8;               % m/s^3  → Δa_x <= 1.6 m/s^2 per step
Jmax_y = 8;               % m/s^3
Jmax_z = 10;              % m/s^3  → Δa_z <= 2.0 m/s^2 per step

% Slew (angular) limit on acceleration vector (direction change per step) — P3 only
dtheta_max_deg = 8;       % deg/step (à 5 Hz ≈ 40 deg/s)
tilt_max_deg   = 30;      % deg vs vertical (+z), P3 only

alpha_acc = 0.0;          % extra LPF after jerk (0..1), typically 0

%% ===== Phase targets and nominal durations =====
% Phase 1 targets (Apollo-like approach start)
x1_t  = 25000;   y1_t  = 0;      z1_t  = 7936.75104830166;
vx1_t = -202;    vy1_t = 0;      vz1_t = -49.09;

Tgo1  = max(10, (x1_t - x) / (0.5*(vx + vx1_t)));

% Phase 2 targets (vertical over site, ~100 m altitude)
x2_t  = 0;     y2_t  = 0;     z2_t  = 100;
vx2_t = 0;     vy2_t = 0;     vz2_t = -8;
Tgo2  = max(10, (x2_t - x1_t) / (0.5*(vx1_t + vx2_t)));
ax2_t = 0;     ay2_t = 0;     az2_t = 1.2*g_moon;

% Phase 3 (Vertical) controller parameters
az3_cmd        = 1.2*g_moon;           % m/s^2 (nominal engine accel)
vx_tol_touch   = 0.5;                  % m/s
vy_tol_touch   = 0.5;                  % m/s
vz_target_touch= -1.0;                 % m/s

%% ===== Logs =====
t = 0; k = 0;
t_all  = []; x_all = []; y_all = []; z_all = [];
vx_all = []; vy_all = []; vz_all = [];
ax_all = []; ay_all = []; az_all = [];
ax_raw_all = []; ay_raw_all = []; az_raw_all = [];
jx_all = []; jy_all = []; jz_all = [];
phase_all = [];

% previous commanded acceleration (for jerk limits)
a_prev = [0 0 0];    % will be initialized to first raw each phase
a_lp   = a_prev;     % for optional extra smoothing

%% ===== Utils =====
clip = @(u,umin,umax) max(min(u,umax),umin);
hasAmax = @(lim) ~isempty(lim) && isfinite(lim);

% Single-axis jerk limiter
function [a_cmd, j_val] = jerk_limit_axis(a_raw, a_prev, Jmax, dt)
    if any(isnan([a_raw, a_prev])), a_raw = 0; end
    da    = a_raw - a_prev;
    daMax = Jmax * dt;
    da    = min(max(da, -daMax), daMax);
    a_cmd = a_prev + da;
    j_val = da / dt;
end

% Vector directional slew (keep magnitude), used ONLY in P3
function a_cmd = slew_dir(a_prev, a_cmd, dtheta_deg)
    mag_prev = max(norm(a_prev), 1e-12);
    mag_cmd  = max(norm(a_cmd),  1e-12);
    n_prev   = a_prev / mag_prev;
    n_cmd    = a_cmd  / mag_cmd;

    c = max(-1,min(1, dot(n_prev, n_cmd)));
    dtheta = acos(c);                             % [rad]
    dtheta_max = deg2rad(dtheta_deg);
    if dtheta > dtheta_max
        s = dtheta_max / dtheta;
        % slerp between n_prev and n_cmd
        n_lim = (sin((1-s)*dtheta)*n_prev + sin(s*dtheta)*n_cmd) / sin(dtheta);
        n_lim = n_lim / max(norm(n_lim),1e-12);
        a_cmd = mag_cmd * n_lim;
    end
end

% Tilt clamp vs +z (keep magnitude approx., reduce lateral if needed) — P3 only
function a_cmd = clamp_tilt(a_cmd, tilt_max_deg)
    if a_cmd(3) <= 1e-9, return; end   % meaningful only for upward (+z) engine accel
    tilt = atan2(hypot(a_cmd(1),a_cmd(2)), a_cmd(3));  % rad
    tilt_max = deg2rad(tilt_max_deg);
    if tilt > tilt_max
        % Reduce lateral keeping az, then renormalize to original magnitude
        mag = norm(a_cmd);
        scale = tan(tilt_max)/max(tan(tilt),1e-12);
        a_cmd(1:2) = a_cmd(1:2) * scale;
        % Recompute to keep magnitude close to original (numerically safe)
        lat2 = a_cmd(1)^2 + a_cmd(2)^2;
        a_cmd(3) = sqrt(max(mag^2 - lat2, 1e-12));
    end
end

%% ===== Coeff functions =====
function [c1, c2] = braking_coeffs(r0, v0, rt, vt, tgo)
    c1 = (-2*(vt+2*v0))/tgo + 6*(rt - r0)/(tgo^2);
    c2 = ( 6*(vt+v0))/(tgo^2) - 12*(rt - r0)/(tgo^3);
end
function [c0,c1,c2] = approach_coeffs(r0,v0,rt,vt,at,tgo)
    c0 = at - 6*(vt+v0)/tgo + 12*(rt - r0)/(tgo^2);
    c1 = -6*at/tgo + 6*(5*vt + 3*v0)/(tgo^2) - 48*(rt - r0)/(tgo^3);
    c2 =  6*at/(tgo^2) - 12*(2*vt + v0)/(tgo^3) + 36*(rt - r0)/(tgo^4);
end

%% ========================= PHASE 1: BRAKING =============================
tgo = Tgo1;
first_step = true;
while z > z1_t
    % raw polynomial command (a = c1 + c2 t → we use a_cmd = c1)
    [c1x,~] = braking_coeffs(x, vx, x1_t, vx1_t, tgo);
    [c1y,~] = braking_coeffs(y, vy, y1_t, vy1_t, tgo);
    [c1z,~] = braking_coeffs(z, vz, z1_t, vz1_t, tgo);
    a_raw = [c1x, c1y, c1z];

    % initialize a_prev at phase start to avoid crushing ax/ay at step 1
    if first_step, a_prev = a_raw; a_lp = a_raw; first_step = false; end

    % optional amplitude clips
    if hasAmax(amax_x), a_raw(1) = clip(a_raw(1), -amax_x, amax_x); end
    if hasAmax(amax_y), a_raw(2) = clip(a_raw(2), -amax_y, amax_y); end
    if hasAmax(amax_z), a_raw(3) = clip(a_raw(3), -amax_z, amax_z); end

    % jerk limit PER AXIS (no tilt/slew here)
    [ax_cmd, jx] = jerk_limit_axis(a_raw(1), a_prev(1), Jmax_x, dt);
    [ay_cmd, jy] = jerk_limit_axis(a_raw(2), a_prev(2), Jmax_y, dt);
    [az_cmd, jz] = jerk_limit_axis(a_raw(3), a_prev(3), Jmax_z, dt);
    a_cmd = [ax_cmd, ay_cmd, az_cmd];
    a_prev = a_cmd; a_lp = alpha_acc*a_lp + (1-alpha_acc)*a_cmd;

    % Integrate dynamics (gravity only on z)
    vx = vx + a_cmd(1)*dt;
    vy = vy + a_cmd(2)*dt;
    vz = vz + (a_cmd(3) - g_moon)*dt;
    x  = x  + vx*dt;
    y  = y  + vy*dt;
    z  = z  + vz*dt;

    % Log
    t  = t + dt;  tgo = tgo - dt;  k = k + 1;
    t_all(k,1)=t;  x_all(k,1)=x;  y_all(k,1)=y;  z_all(k,1)=z;
    vx_all(k,1)=vx; vy_all(k,1)=vy; vz_all(k,1)=vz;
    ax_all(k,1)=a_cmd(1); ay_all(k,1)=a_cmd(2); az_all(k,1)=a_cmd(3);
    ax_raw_all(k,1)=a_raw(1); ay_raw_all(k,1)=a_raw(2); az_raw_all(k,1)=a_raw(3);
    jx_all(k,1)=jx; jy_all(k,1)=jy; jz_all(k,1)=jz;
    phase_all(k,1)=1;

    if z <= 0, break; end
end
fprintf('stop phase 1 iteration : %d\n', size(z_all,1));

%% ========================= PHASE 2: APPROACH ============================
tgo = Tgo2;
first_step = true;
while z > z2_t
    % raw polynomial command (a = c0 + c1 t + c2 t^2 → we use a_cmd = c0)
    [c0x,~,~] = approach_coeffs(x, vx, x2_t, vx2_t, ax2_t, tgo);
    [c0y,~,~] = approach_coeffs(y, vy, y2_t, vy2_t, ay2_t, tgo);
    [c0z,~,~] = approach_coeffs(z, vz, z2_t, vz2_t, az2_t, tgo);
    a_raw = [c0x, c0y, c0z];

    if first_step, a_prev = a_raw; a_lp = a_raw; first_step = false; end

    % optional amplitude clips
    if hasAmax(amax_x), a_raw(1) = clip(a_raw(1), -amax_x, amax_x); end
    if hasAmax(amax_y), a_raw(2) = clip(a_raw(2), -amax_y, amax_y); end
    if hasAmax(amax_z), a_raw(3) = clip(a_raw(3), -amax_z, amax_z); end

    % jerk limit PER AXIS (no tilt/slew here)
    [ax_cmd, jx] = jerk_limit_axis(a_raw(1), a_prev(1), Jmax_x, dt);
    [ay_cmd, jy] = jerk_limit_axis(a_raw(2), a_prev(2), Jmax_y, dt);
    [az_cmd, jz] = jerk_limit_axis(a_raw(3), a_prev(3), Jmax_z, dt);
    a_cmd = [ax_cmd, ay_cmd, az_cmd];
    a_prev = a_cmd; a_lp = alpha_acc*a_lp + (1-alpha_acc)*a_cmd;

    % Integrate dynamics
    vx = vx + a_cmd(1)*dt;
    vy = vy + a_cmd(2)*dt;
    vz = vz + (a_cmd(3) - g_moon)*dt;
    x  = x  + vx*dt;
    y  = y  + vy*dt;
    z  = z  + vz*dt;

    % Log
    t = t + dt; tgo = tgo - dt; k = k + 1;
    t_all(k,1)=t;  x_all(k,1)=x;  y_all(k,1)=y;  z_all(k,1)=z;
    vx_all(k,1)=vx; vy_all(k,1)=vy; vz_all(k,1)=vz;
    ax_all(k,1)=a_cmd(1); ay_all(k,1)=a_cmd(2); az_all(k,1)=a_cmd(3);
    ax_raw_all(k,1)=a_raw(1); ay_raw_all(k,1)=a_raw(2); az_raw_all(k,1)=a_raw(3);
    jx_all(k,1)=jx; jy_all(k,1)=jy; jz_all(k,1)=jz;
    phase_all(k,1)=2;

    if z <= 0, break; end
end
fprintf('stop phase 2 iteration : %d\n', size(z_all,1));

%% ========================= PHASE 3: VERTICAL ============================
% Simple vertical descent controller + jerk + directional slew + tilt.
Kxy = 0.8;                     % horizontal damping
ax_max = 2.0; ay_max = 2.0;    % caps for lateral control
az_max = az3_cmd + 0.6*g_moon; % bound on commanded engine accel

v_td   = -0.1;                 % final target velocity near ground
a_up_max = max(az_max - g_moon, 0.1);
alpha_f  = 0.6;                % legacy LPF for az computation (kept)
az_cmd_filt = az3_cmd;

first_step = true;
while z > 0
    % (1) Horizontal “deadband” controller
    if z > 50
        vx_lim = 10; vy_lim = 10;
        ax_h = (abs(vx)>vx_lim) * (-Kxy*(abs(vx)-vx_lim)*sign(vx));
        ay_h = (abs(vy)>vy_lim) * (-Kxy*(abs(vy)-vy_lim)*sign(vy));
    else
        band = 0.1;
        ax_h = (vx> band)*(-Kxy*(vx-band)) + (vx<-band)*(-Kxy*(vx+band));
        ay_h = (vy> band)*(-Kxy*(vy-band)) + (vy<-band)*(-Kxy*(vy+band));
    end
    ax_h = clip(ax_h, -ax_max, ax_max);
    ay_h = clip(ay_h, -ay_max, ay_max);

    % (2) Vertical reference profile (stoppable envelope)
    v_env = -sqrt(max(v_td^2 + 2*a_up_max*max(z,0), 0));
    if     z > 200, v_cap = v_env;
    elseif z > 50,  v_cap = max(v_env, -2);
    elseif z > 5,   v_cap = max(v_env, -1);
    else            v_cap = max(v_env, v_td);
    end
    v_ref = v_cap;

    % base time constant
    if     z > 200, tau = 12.0;
    elseif z > 100, tau =  6.0;
    elseif z > 50,  tau =  4.0;
    elseif z > 20,  tau =  2.0;
    elseif z > 5,   tau =  1.5;
    else            tau =  2.0;
    end

    a_net_des  = (v_ref - vz)/tau;           % desired net (after gravity)
    az_raw     = g_moon + a_net_des;         % engine accel command
    az_raw     = clip(az_raw, 0, az_max);

    % a_raw vector (before shaping)
    a_raw = [ax_h, ay_h, az_raw];
    if first_step, a_prev = a_raw; a_lp = a_raw; first_step = false; end

    % optional amplitude clips by axis
    if hasAmax(amax_x), a_raw(1) = clip(a_raw(1), -amax_x, amax_x); end
    if hasAmax(amax_y), a_raw(2) = clip(a_raw(2), -amax_y, amax_y); end
    if hasAmax(amax_z), a_raw(3) = clip(a_raw(3), 0, amax_z);       end

    % legacy vertical LPF then jerk limit per axis
    az_cmd_filt = alpha_f*az_cmd_filt + (1-alpha_f)*a_raw(3);
    a_raw(3)    = az_cmd_filt;

    % jerk limit per axis
    [ax_cmd, jx] = jerk_limit_axis(a_raw(1), a_prev(1), Jmax_x, dt);
    [ay_cmd, jy] = jerk_limit_axis(a_raw(2), a_prev(2), Jmax_y, dt);
    [az_cmd, jz] = jerk_limit_axis(a_raw(3), a_prev(3), Jmax_z, dt);
    a_cmd = [ax_cmd, ay_cmd, az_cmd];

    % THEN vector slew + tilt clamp (now az_cmd >= 0)
    a_cmd = slew_dir(a_prev, a_cmd, dtheta_max_deg);
    a_cmd = clamp_tilt(a_cmd, tilt_max_deg);

    a_prev = a_cmd; a_lp = alpha_acc*a_lp + (1-alpha_acc)*a_cmd;

    % Integrate
    vx = vx + a_cmd(1)*dt;
    vy = vy + a_cmd(2)*dt;
    vz = vz + (a_cmd(3) - g_moon)*dt;
    x  = x  + vx*dt;
    y  = y  + vy*dt;
    z  = z  + vz*dt;

    % Log
    t = t + dt; k = k + 1;
    t_all(k,1)=t;  x_all(k,1)=x;  y_all(k,1)=y;  z_all(k,1)=max(z,0);
    vx_all(k,1)=vx; vy_all(k,1)=vy; vz_all(k,1)=vz;
    ax_all(k,1)=a_cmd(1); ay_all(k,1)=a_cmd(2); az_all(k,1)=a_cmd(3);
    ax_raw_all(k,1)=a_raw(1); ay_raw_all(k,1)=a_raw(2); az_raw_all(k,1)=a_raw(3);
    jx_all(k,1)=jx; jy_all(k,1)=jy; jz_all(k,1)=jz;
    phase_all(k,1)=3;
end

% Fix last ground point
if z_all(end) < 0, z_all(end) = 0; end

%% ===== ΔV & G-proper (based on engine accel) =====
a_mag  = sqrt(ax_all.^2 + ay_all.^2 + az_all.^2);
DV_total = sum(a_mag)*dt;
DV_P1 = sum(a_mag(phase_all==1))*dt;
DV_P2 = sum(a_mag(phase_all==2))*dt;
DV_P3 = sum(a_mag(phase_all==3))*dt;

g0 = 9.80665;
G_prop = a_mag / g0;
[maxG, iG] = max(G_prop);
fprintf('Total ΔV = %.1f m/s  (P1=%.1f, P2=%.1f, P3=%.1f)\n', DV_total, DV_P1, DV_P2, DV_P3);
fprintf('Max proper G = %.2f g at t = %.1f s (phase %d)\n', maxG, t_all(iG), phase_all(iG));

for p = 1:3
    idxp = (phase_all==p);
    if any(idxp)
        [mx, ii] = max(G_prop(idxp));
        tmx = t_all(find(idxp,1,'first') + ii - 1);
        fprintf('  Phase %d : max = %.2f g at t = %.1f s\n', p, mx, tmx);
    end
end

%% ===== Thrust/mass estimation (optional) =====
params.g_moon = g_moon;
params.Isp    = 450;
params.g0     = g0;
params.Tmax   = 100e3;      % N
params.m0     = 18000;      % kg
params.m_dry  = 11000;      % kg
params.tilt_max_deg = tilt_max_deg;

V = [vx_all(:), vy_all(:), vz_all(:)];
out = thrust_from_traj(t_all(:), V, params);
fprintf('Delta-V (integrated) = %.1f m/s, Propellant = %.1f kg\n', out.dV, out.m_prop);
fprintf('Max throttle = %.0f %%\n', 100*max(out.throttle));

%% ===== Plots =====
figure('Name','3D Trajectory & Time Series');

subplot(2,2,1); hold on; grid on;
plot3(x_all/1000, y_all/1000, z_all/1000, 'k-');
phc = lines(3);
for ph=1:3
    idxp = (phase_all==ph);
    plot3(x_all(idxp)/1000, y_all(idxp)/1000, z_all(idxp)/1000, '.', 'Color', phc(ph,:), 'MarkerSize', 6);
end
set(gca,'XDir','reverse'); xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)'); view(135,25);
title('3D trajectory'); legend('Full','P1','P2','P3','Location','best');

subplot(2,2,2); hold on; grid on;
plot(t_all, z_all/1000, 'LineWidth',1.2);
xlabel('t (s)'); ylabel('Altitude (km)'); title('Altitude vs time');

subplot(2,2,3); hold on; grid on;
plot(t_all, ax_raw_all,'--'); plot(t_all, ay_raw_all,'--'); plot(t_all, az_raw_all,'--');
plot(t_all, ax_all,'-','LineWidth',1.2);
plot(t_all, ay_all,'-','LineWidth',1.2);
plot(t_all, az_all,'-','LineWidth',1.2);
xlabel('t (s)'); ylabel('a_{cmd} (m/s^2)'); title('Acceleration: raw (--) vs shaped (-)');
legend('a_x raw','a_y raw','a_z raw','a_x','a_y','a_z','Location','best');

subplot(2,2,4); hold on; grid on;
plot(t_all, jx_all, t_all, jy_all, t_all, jz_all, 'LineWidth',1.1);
yline(Jmax_x,'--'); yline(-Jmax_x,'--');
xlabel('t (s)'); ylabel('jerk (m/s^3)'); title('Realized jerk (per axis)');
legend('j_x','j_y','j_z','\pmJmax_x','Location','best');

figure('Name','Speeds');
plot(t_all, vx_all, t_all, vy_all, t_all, vz_all, 'LineWidth',1.1); grid on;
yline(0,'--'); xlabel('t (s)'); ylabel('V (m/s)'); title('Speed'); legend('v_x','v_y','v_z','Location','best');

fprintf('--- Touchdown ---\n');
fprintf('t=%.1f s | x=%.1f m | y=%.1f m | z=%.2f m | vx=%.2f | vy=%.2f | vz=%.2f m/s\n', ...
    t_all(end), x_all(end), y_all(end), z_all(end), vx_all(end), vy_all(end), vz_all(end));

%% ====== Helper: thrust_from_traj =======================================
function out = thrust_from_traj(t, V, params)
    if size(V,2)==2, V = [V, zeros(size(V,1),1)]; end
    g    = getfielddef(params,'g_moon',1.62);
    Isp  = getfielddef(params,'Isp',450);
    g0   = getfielddef(params,'g0',9.80665);
    Tmax = getfielddef(params,'Tmax',1e5);
    m0   = params.m0;
    m_dry= getfielddef(params,'m_dry',0);

    vx = V(:,1); vy = V(:,2); vz = V(:,3);
    ax = gradient(vx, t);
    ay = gradient(vy, t);
    az = gradient(vz, t);

    % engine acceleration command ≈ a + [0 0 g]
    a_cmd = [ax, ay, az + g];
    a_cmd(:,3) = max(a_cmd(:,3), 0);

    amag  = sqrt(sum(a_cmd.^2,2));
    kappa = amag./(Isp*g0);
    int_kappa = cumtrapz(t, kappa);
    m = m0 .* exp(-int_kappa);
    if m_dry>0, m = max(m, m_dry); end

    T_vec = a_cmd .* m;
    T_mag = m .* amag;
    throttle = T_mag / Tmax;

    dV     = trapz(t, amag);
    m_prop = m0 - m(end);

    out.T_vec   = T_vec;
    out.T_mag   = T_mag;
    out.throttle= throttle;
    out.m       = m;
    out.dV      = dV;
    out.m_prop  = m_prop;
    out.a_cmd   = a_cmd;
    out.flags.sat_throttle = throttle>1;
end
function v = getfielddef(s, name, def)
    if isfield(s,name) && ~isempty(s.(name)), v = s.(name); else, v = def; end
end
