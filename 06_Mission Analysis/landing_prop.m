%% Lunar Powered Descent Guidance (3D) — braking, approach, vertical
% Implements the polynomial guidance laws from Sostaric & Rea
% - Phase 1 (Braking)
% - Phase 2 (Approach)
% - Phase 3 (Vertical

clear; clc;

%% ===== General parameters =====
g_moon = 1.62;        % m/s^2
dt      = 0.2;        % integration step (s) — can be increased if t_go is long
rng('default');

% Initial state (start of PHASE 1 - Braking). Adjust to your case.
x0  = 289448.1899567346;        % initial downrange to target (m)
z0  = 15e3; 
y0  = 0e3; 
% initial altitude (m)
vx0 = -1695;          % horizontal speed toward target (m/s) (negative if x decreases)
vz0 = 0;  
vy0=0 ;               % vertical descending speed (m/s, negative downward)

% Optional command limits (disable by setting [] )
amax_x = [];          % m/s^2  (e.g., 5)
amax_z = [];          % m/s^2  (e.g., 5)
amax_y = [];
difax_all = [];
%% ===== Phase targets and durations (adapt as needed) =====
% Phase 1 (Braking) -> conditions close to Apollo-like approach start
x1_t  = 25000;   y1_t  = 0;     z1_t  = 7936.75104830166   % 
vx1_t = -202;    vy1_t =  0;    vz1_t = -49.09;      % m/s

Tgo1  = max(10, (x1_t - x0) / (0.5*(vx0 + vx1_t)));         % s (allocated duration for braking phase)

% Phase 2 (Approach) -> bring vehicle vertical over site at 100 m
x2_t   = 0;           % m  (vertical over the site)
z2_t   = 100;        % m
y2_t=0;

vx2_t  = 0;           % m/s (cancel horizontal speed)
vy2_t=0;
vz2_t  = -8;        % m/s (target for entering vertical phase)
Tgo2  = max(10, (x2_t - x1_t) / (0.5*(vx1_t + vx2_t)));
ax2_t  = 0;           % m/s^2 (target zero horizontal acceleration)
ay2_t=0;
az2_t  = 1.2*g_moon;  % m/s^2 (target vertical engine accel ~1.2 g_moon)

% Phase 3 (Vertical) -> final constant descent
az3_cmd = 1.2*g_moon; % m/s^2
vx_tol_touch = 0.5;   % m/s (horizontal tolerance at touchdown)
vy_tol_touch=0.5;
vz_target_touch = -1.0; % m/s (target touchdown speed ~ -1 m/s)


%% ===== Utilities =====
clip = @(u,umin,umax) max(min(u,umax),umin);

function [c1, c2] = braking_coeffs(r0, v0, rt, vt, tgo)
    c1 = (-2*(vt+2*v0))/tgo + 6*(rt - r0)/(tgo^2);
    c2 = ( 6*(vt+v0))/(tgo^2) - 12*(rt - r0)/(tgo^3);
end

function [c0,c1,c2] = approach_coeffs(r0,v0,rt,vt,at,tgo)
    c0 = at - 6*(vt+v0)/tgo + 12*(rt - r0)/(tgo^2);
    c1 = -6*at/tgo + 6*(5*vt + 3*v0)/(tgo^2) - 48*(rt - r0)/(tgo^3);
    c2 =  6*at/(tgo^2) - 12*(2*vt + v0)/(tgo^3) + 36*(rt - r0)/(tgo^4);
end
%% 

% -------- PHASE 1: Braking (a = c1 + c2 t, we apply a_cmd = c1) ------
% Current state
t = 0; x = x0; y=y0; z = z0; vx = vx0; vy=vy0 ; vz = vz0;
% Write index
k = 0;
ax_cmd=-1
% -------- PHASE 1: Braking --------
tgo = Tgo1;
while (z > z1_t)
    [c1x,~] = braking_coeffs(x, vx, x1_t, vx1_t, tgo);
    [c1y,~] = braking_coeffs(y, vy, y1_t, vy1_t, tgo);
    [c1z,~] = braking_coeffs(z, vz, z1_t, vz1_t, tgo);
    ax_cmd = c1x;
    ay_cmd = c1y;
    az_cmd = c1z;

    

    if ~isempty(amax_x), ax_cmd = clip(ax_cmd,-amax_x,amax_x); end
    if ~isempty(amax_y), ay_cmd = clip(ay_cmd,-amax_y,amax_y); end
    if ~isempty(amax_z), az_cmd = clip(az_cmd,-amax_z,amax_z); end

    % Dynamics
    vx = vx + ax_cmd*dt;
    vy = vy + ay_cmd*dt;
    vz = vz + (az_cmd - g_moon)*dt;
    x  = x  + vx*dt;
    y= y + vy*dt;
    z  = z  + vz*dt;

    % Log (without end+1)
    t    = t + dt;  tgo = tgo - dt;  k = k + 1;
    t_all(k,1) = t;    x_all(k,1) = x; y_all(k,1) = y ;  z_all(k,1) = z;
    vx_all(k,1)= vx;  vy_all(k,1)=vy;  vz_all(k,1)= vz;
    ax_all(k,1)= ax_cmd; ay_all(k,1) = ay_cmd; az_all(k,1)= az_cmd;
    phase_all(k,1)=1;

    if z <= 0, break; end
end
fprintf('stop phase 1 iteration : %d\n', size(z_all,1));
%% 
% ----- Filtres d'orientation P2 (état et paramètres) -----
theta_rate_max_deg = 30;           % vitesse max d'angle (deg/s) — ajuste 10..60
alpha_track        = 0.2;          % 0..1, réactivité hors sursaut (0.2 à 0.5 ok)
jump_deg           = 12;           % seuil pour considérer un "gros sursaut" (deg)

theta_f_prev = [];                 % état interne (initialisé au 1er pas)

% -------- PHASE 2: Approach (a = c0 + c1 t + c2 t^2, we apply a_cmd = c0) --------
tgo = Tgo2;
jerk0=0
% Current state
while  (z > z2_t)
    % Coefficients per axis
    [c0x,~,~] = approach_coeffs(x, vx, x2_t, vx2_t, ax2_t, tgo);
    [c0y,~,~] = approach_coeffs(y, vy, y2_t, vy2_t, ay2_t, tgo);
    [c0z,~,~] = approach_coeffs(z, vz, z2_t, vz2_t, az2_t, tgo);
    
    ax_cmd = c0x;  ay_cmd = c0y;  az_cmd = c0z;
    

    
    if ~isempty(amax_x), ax_cmd = clip(ax_cmd,-amax_x,amax_x); end
    if ~isempty(amax_y), ay_cmd = clip(ay_cmd,-amax_y,amax_y); end
    if ~isempty(amax_z), az_cmd = clip(az_cmd,-amax_z,amax_z); end
    
    
    vx = vx + ax_cmd*dt;   vy = vy + ay_cmd*dt;
    vz = vz + (az_cmd - g_moon)*dt;
    x  = x  + vx*dt;       y  = y  + vy*dt;      z = z + vz*dt;

    t = t + dt; tgo = tgo - dt;
    t_all(end+1,1)=t; x_all(end+1,1)=x; y_all(end+1,1)=y; z_all(end+1,1)=z;
    vx_all(end+1,1)=vx; vy_all(end+1,1)=vy; vz_all(end+1,1)=vz;
    ax_all(end+1,1)=ax_cmd; ay_all(end+1,1)=ay_cmd; az_all(end+1,1)=az_cmd;
    phase_all(end+1,1)=2;
    if z <= 0, break; end
end

fprintf('stop phase 2 iteration : %d\n', size(z_all,1));
%% 

% -------- PHASE 3: Vertical + DIVERT strategy ----------------------------
% (1) Unchanged parameters
if ~exist('Kxy','var'),      Kxy = 0.8; end
if ~exist('ax_max','var'),   ax_max = 2.0; end
if ~exist('ay_max','var'),   ay_max = 2.0; end
if ~exist('daz_max','var'),  daz_max = 0.6*g_moon; end
if ~exist('az_cmd_filt','var'), az_cmd_filt = az3_cmd; end

v_td   = -0.1;                 % m/s target at impact
az_nom = az3_cmd;              % m/s^2 (1.2 g_moon)
az_max = az_nom + daz_max;     % thrust bound
a_up_max = max(az_max - g_moon, 0.1);  % net braking capability
alpha_f  = 0.6;                % smoothing 0..1

% (2) DIVERT table (from your table)
divert_D   = [25   50    75   100   125   150];      % m
divert_dT  = [-0.9 -0.6  4.0   8.3   11.1  15.4];    % s
divert_dDV = [-0.9  0.6  8.2  15.5   20.3  27.4];    % m/s (info)

% Additional remaining time (decreases over time)
dT_hold = 0;    % s
dDV_est = 0;    % m/s (optional, for log)

while z > 0
    % ================== 1) Horizontal correctors ==================
    if z > 50
        vx_lim = 10; vy_lim = 10;      % limits under 50 m
        ax_cmd = (abs(vx)>vx_lim) * (-Kxy*(abs(vx)-vx_lim)*sign(vx));
        ay_cmd = (abs(vy)>vy_lim) * (-Kxy*(abs(vy)-vy_lim)*sign(vy));
    else
        band = 0.1;                    % deadband ±0.1 m/s under 10 m
        ax_cmd = (vx> band)*(-Kxy*(vx-band)) + (vx<-band)*(-Kxy*(vx+band));
        ay_cmd = (vy> band)*(-Kxy*(vy-band)) + (vy<-band)*(-Kxy*(vy+band));
    end
    ax_cmd = clip(ax_cmd,-ax_max,ax_max);
    ay_cmd = clip(ay_cmd,-ay_max,ay_max);

    % ================== 2) DIVERT management (lateral distance) ========
    D = hypot(x,y);                                   % m
    dT_target = interp1(divert_D, divert_dT, D, 'linear','extrap');
    dDV_target = interp1(divert_D, divert_dDV, D, 'linear','extrap');

    % Never shorten the descent: keep only ΔT>0
    dT_target = max(0, dT_target);
    % Keep the maximum time need observed and consume it progressively
    dT_hold = max(dT_hold, dT_target);
    dDV_est = max(dDV_est, max(0,dDV_target));       % info (not used in loop)

    % Rough estimate of remaining time (used to normalize stretching)
    tgo_est = z / max(0.5, -vz + 0.1);               % s, robust if vz≈0
    slow_factor = 1 + min(0.8, dT_hold / max(5, tgo_est));   % ≤ +80% on τ

    % ================== 3) Vertical profile v_ref(z) ==================
    v_env = -sqrt(max(v_td^2 + 2*a_up_max*max(z,0), 0));   % stoppable envelope
    if     z > 200, v_cap = v_env;
    elseif z > 50,  v_cap = max(v_env, -2);   % intermediate speed limit
    elseif z > 5,  v_cap = max(v_env, -1);
    else            v_cap = max(v_env, v_td); % converge to -0.2 m/s
    end
    v_ref = v_cap;

    % Base time constant (tighter near ground), then “stretched” by DIVERT
    if     z > 200, tau = 12.0;
    elseif z > 100, tau = 6.0;
    elseif z > 50,  tau = 4.0;
    elseif z > 20,  tau = 2.0;
    elseif z > 5,  tau = 1.5;
    else            tau = 2;
    end
    tau = tau * slow_factor;     % <<< buy time for DIVERT

    % Desired net acceleration & commanded thrust
    a_net_des = (v_ref - vz)/tau;                   % m/s^2
    az_cmd_raw = g_moon + a_net_des;                % m->thrust
    az_cmd_raw = clip(az_cmd_raw, 0, az_max);
    az_cmd = alpha_f*az_cmd_filt + (1-alpha_f)*az_cmd_raw;
    az_cmd_filt = az_cmd;

    % ================== 4) Integration + log =========================
    vx = vx + ax_cmd*dt;   vy = vy + ay_cmd*dt;
    vz = vz + (az_cmd - g_moon)*dt;
    x  = x  + vx*dt;       y  = y  + vy*dt;        z = z + vz*dt;

    t = t + dt;
    t_all(end+1,1)=t; x_all(end+1,1)=x; y_all(end+1,1)=y; z_all(end+1,1)=z;
    vx_all(end+1,1)=vx; vy_all(end+1,1)=vy; vz_all(end+1,1)=vz;
    ax_all(end+1,1)=ax_cmd; ay_all(end+1,1)=ay_cmd; az_all(end+1,1)=az_cmd;
    phase_all(end+1,1)=3;

    % Consume the “bought” time stock
    dT_hold = max(0, dT_hold - dt);
end

% Fix last ground point
if z_all(end) < 0, z_all(end) = 0; end

% (optional) console trace of adopted strategy
fprintf('[DIVERT] D=%.1f m → target ΔT≈%.1f s, add. ΔV≈%.1f m/s (table)\n', ...
        D, dT_target, dDV_target);

%% ===== ΔV & proper G in 3D =====
idx = (phase_all > 0);
a_thrust = sqrt(ax_all(idx).^2 + ay_all(idx).^2 + az_all(idx).^2);
dv_step  = a_thrust * dt;
DV_total = sum(dv_step);
DV_P1 = sum(dv_step(phase_all(idx)==1));
DV_P2 = sum(dv_step(phase_all(idx)==2));
DV_P3 = sum(dv_step(phase_all(idx)==3));
fprintf('Total ΔV = %.1f m/s  (P1=%.1f, P2=%.1f, P3=%.1f)\n', DV_total, DV_P1, DV_P2, DV_P3);

g0 = 9.80665;
G_prop = sqrt(ax_all.^2 + ay_all.^2 + az_all.^2) / g0;
[maxG, iG] = max(G_prop);
fprintf('Max proper G = %.2f g at t = %.1f s (phase %d)\n', maxG, t_all(iG), phase_all(iG));

% Max per phase
for p = 1:3
    idx = (phase_all == p);
    if any(idx)
        [mx, ii] = max(G_prop(idx));
        tmx = t_all(find(idx,1,'first') + ii - 1);
        fprintf('  Phase %d : max = %.2f g at t = %.1f s\n', p, mx, tmx);
    end
end

%%
function out = thrust_from_traj(t, V, params)
% THRUST_FROM_TRAJ  Required thrust to follow a given trajectory.
% Axes: z = vertical (up +), x = horizontal (downrange), y = lateral.
%
% Inputs:
%   t      [Nx1]  time (s), not necessarily uniform
%   V      [Nx2] or [Nx3] velocities: columns = [vx vy (vz)]
%   params struct:
%       .g_moon       (default 1.62)
%       .Isp          (s) (default 450)
%       .g0           (default 9.80665)
%       .Tmax         (N) max thrust (e.g. 1e5)
%       .m0           (kg) initial (wet) mass   [required]
%       .m_dry        (kg) dry mass (optional, clamp)
%       .tilt_max_deg (optional) tilt limit (e.g. 25)
%
% Outputs:
%   T_vec   [Nx3]  thrust vector [Tx Ty Tz] (N)  (Tz = vertical)
%   T_mag   [Nx1]  thrust magnitude (N)
%   throttle[Nx1]  fraction of Tmax
%   m       [Nx1]  mass (kg)
%   dV      (1x1)  Delta-V thrust = ∫||a_cmd|| dt (m/s)
%   m_prop  (1x1)  propellant consumed (kg)
%   a_cmd   [Nx3]  required thrust acceleration (m/s^2)
%   tiltdeg [Nx1]  tilt vs vertical (deg)
%   flags   struct  indicators (saturation, tilt, Tvert<0)

    % --- Params & defaults
    if size(V,2)==2, V = [V, zeros(size(V,1),1)]; end  % assume vz=0 if missing
    g    = getfielddef(params,'g_moon',1.62);
    Isp  = getfielddef(params,'Isp',450);
    g0   = getfielddef(params,'g0',9.80665);
    Tmax = getfielddef(params,'Tmax',1e5);
    m0   = params.m0;
    m_dry= getfielddef(params,'m_dry',0);
    tiltMaxDeg = getfielddef(params,'tilt_max_deg',inf);

    % --- Accelerations (robust numerical derivative)
    vx = V(:,1); vy = V(:,2); vz = V(:,3);
    ax = gradient(vx, t);
    ay = gradient(vy, t);
    az = gradient(vz, t);

    % --- a_cmd = a + [0 0 g]  (lunar gravity downward)
    
    a_cmd = [ax, ay, az + g];  
    a_cmd(:,3) = max(a_cmd(:,3), 0);          % forbid "downward" thrust
    
    amag  = sqrt(sum(a_cmd.^2,2));
    
    % --- m(t) closed form: m = m0 * exp(-∫ (||a_cmd||/(Isp*g0)) dt )
    kappa = amag./(Isp*g0);
    int_kappa = cumtrapz(t, kappa);
    m = m0 .* exp(-int_kappa);
    if m_dry>0, m = max(m, m_dry); end

    % --- Required thrust
    T_vec = a_cmd .* m;              % [Tx Ty Tz]
    epsT = 1e-6 * params.Tmax;   % tolerance
    vertNeg = T_vec(:,3) < -epsT;
    T_mag = m .* amag;
    throttle = T_mag / Tmax;

    % --- Attitude/tilt checks (reference = vertical Tz)
    Tx = T_vec(:,1); Ty = T_vec(:,2); Tz = T_vec(:,3);
    Th = hypot(Tx,Ty);                               % in-plane component
    tiltdeg = atan2d(Th, max(Tz,0));                 % tilt vs vertical (0° = straight up)
    vertNeg = Tz < 0;                                 % negative vertical thrust (not realistic)
    overTilt= tiltdeg > tiltMaxDeg;                   % tilt beyond limit
    overTh  = throttle > 1;                           % exceeds Tmax

    % --- Totals
    dV     = trapz(t, amag);
    m_prop = m0 - m(end);

    % --- Pack
    out.T_vec   = T_vec;
    out.T_mag   = T_mag;
    out.throttle= throttle;
    out.m       = m;
    out.dV      = dV;
    out.m_prop  = m_prop;
    out.a_cmd   = a_cmd;
    out.tiltdeg = tiltdeg;
    out.flags.sat_throttle = overTh;
    out.flags.vert_negative= vertNeg;
    out.flags.over_tilt    = overTilt;
end

function v = getfielddef(s, name, def)
    if isfield(s,name) && ~isempty(s.(name)), v = s.(name); else, v = def; end
end

% Trajectory data:
params.g_moon = 1.62;
params.Isp    = 450;
params.g0     = 9.80665;
params.Tmax   = 100e3;
params.m0     = 18000;      % 
params.m_dry  = 11000;      % 11 t
params.tilt_max_deg = 25;   % attitude constraint (optional)

V = [vx_all(:), vy_all(:), vz_all(:)];   % [vx vy vz] with z vertical
out = thrust_from_traj(t_all(:), V, params);
out.throttle(end)=0

fprintf('Delta-V = %.1f m/s, Propellant = %.1f kg\n', out.dV, out.m_prop);
fprintf('Max throttle = %.0f %%\n', 100*max(out.throttle));

Massepropellant=(params.m0-params.m_dry)+(out.m-params.m0);
m_prop = params.m_dry * (exp(DV_total/(params.Isp*g0)) - 1)
if m_prop>(params.m0-params.m_dry);
    warning('Not enough propellant to land')
end
if any(out.flags.sat_throttle)
    warning('Trajectory requires >100%% thrust at some instants.');
end
if any(out.flags.vert_negative)
    warning('Negative vertical thrust component (unrealistic for single engine).');
end

%% ===== Plots =====



figure('Name','Throttle and propellant')
subplot(2,2,1);hold on; grid on;
plot(t_all,out.throttle*100)
xlabel('t (s)'); ylabel('Throttle (%)'); title('Throttle vs time');

subplot(2,2,2);hold on; grid on;
plot(t_all,Massepropellant)
xlabel('t (s)'); ylabel('Propellant Mass'); title('Propellant Mass vs time');
figure('Name','3D Trajectory');
subplot(2,2,1); hold on; grid on;
plot3(x_all/1000, y_all/1000, z_all/1000, 'Color',[0.4 0.4 0.4]); 
ph_colors = lines(3);
for ph=1:3
    idxp = (phase_all==ph);
    plot3(x_all(idxp)/1000, y_all(idxp)/1000, z_all(idxp)/1000, '.', 'Color', ph_colors(ph,:), 'MarkerSize', 6);
end
set(gca,'XDir','reverse'); xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
title('3D trajectory'); view(135,25);

subplot(2,2,2); hold on; grid on;
plot(t_all, z_all/1000, 'LineWidth',1.2);
xlabel('t (s)'); ylabel('Altitude z (km)'); title('Altitude vs time');

% 3) Commanded accelerations
subplot(2,2,3); hold on; grid on;
plot(t_all, ax_all, 'LineWidth',1.5);
plot(t_all, ay_all, 'LineWidth',1.5);
plot(t_all,az_all, 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('a_{cmd} (m/s^2)');
title('Acceleration Command');
legend('a_x','a_y','a_z (engine)','Location','best');


subplot(2,2,4); hold on; grid on;
plot(t_all, vx_all, t_all, vy_all, t_all, vz_all, 'LineWidth',1.1);
yline(0,'--'); xlabel('t (s)'); ylabel('V (m/s)'); title('Speed'); legend('v_x','v_y','v_z');

fprintf('--- Touchdown ---\n');
fprintf('t=%.1f s | x=%.1f m | y=%.1f m | z=%.2f m | vx=%.2f | vy=%.2f | vz=%.2f m/s\n', ...
    t_all(end), x_all(end), y_all(end), z_all(end), vx_all(end), vy_all(end), vz_all(end));

%% ===== Thrust angle vs surface / vertical =====
% Thrust vector in world coordinates (your recorded a_cmd)
Tnorm = sqrt(ax_all.^2 + ay_all.^2 + az_all.^2);

% Angle with vertical (surface normal = +z)
ang_vert = nan(size(Tnorm));                     % [deg]
mask = Tnorm > 1e-9;                             % avoid division by 0
ang_vert(mask) = acosd( az_all(mask) ./ Tnorm(mask) );

% Angle relative to surface (horizontal plane)
ang_surface = 90 - ang_vert;                     % [deg], 0° = flat, 90° = vertical

% (optional) thrust azimuth in horizontal plane
azimut_deg = atan2d(ay_all, ax_all);             % [-180,180]

% Useful stats
[maxTilt, iMax] = max(ang_surface);
fprintf('Max tilt vs surface = %.2f deg at t = %.1f s, z = %.1f m (phase %d)\n', ...
        maxTilt, t_all(iMax), z_all(iMax), phase_all(iMax));

% Save for reuse
ang_vert_all    = ang_vert;
ang_surface_all = ang_surface;
azimut_all      = azimut_deg;

%% === Correction du profil d'accélérations ===
raw = ang_vert;   % ton signal angle to vertical
t   = t_all;      

% --- Zone à corriger
t1 = 110; 
t2 = 171;

i1 = find(t >= t1, 1, 'first');
i2 = find(t <= t2, 1, 'last');

smooth = raw;
smooth(i1:i2) = linspace(raw(i1), raw(i2), i2-i1+1);

% === Reconstruction des accélérations corrigées ===
Tnorm = sqrt(ax_all.^2 + ay_all.^2 + az_all.^2); % norme originale
theta = deg2rad(smooth);                         % angle corrigé
phi   = atan2(ay_all, ax_all);                   % azimut inchangé

ax_filt = Tnorm .* sin(theta) .* cos(phi);
ay_filt = Tnorm .* sin(theta) .* sin(phi);
az_filt = Tnorm .* cos(theta);


% === On remplace les profils originaux par les corrigés ===
ax_all = ax_filt;
ay_all = ay_filt;
az_all = az_filt;

% (optionnel: mettre à jour aussi ang_vert si tu l'utilises plus tard)
ang_vert = smooth;


%% ================== 3D DESCENT ANIMATION and GAUGES (with lunar ground) ======

% Uses t_all, x_all, y_all, z_all, phase_all
% And for gauges: out.m, out.throttle (fallback m_all, thr_all)

Xkm = x_all/1000;  
Ykm = Xkm;     % <<< fixes old Ykm = Xkm
Zkm = z_all/1000;
N   = numel(t_all);

% -------- Gauge data (mass & throttle) ------------------------------
% m_dry must exist (otherwise use 0)
if ~exist('m_dry','var'), m_dry = 0; end

hasOut = exist('out','var') && isstruct(out);
if hasOut && isfield(out,'m'),        m_series = Massepropellant(:); 
elseif exist('m_all','var'),          m_series = m_all(:);
else,                                 m_series = repmat(m_dry, N, 1);
end

if hasOut && isfield(out,'throttle'), thr_series = out.throttle(:);thr_series(end)=0;
elseif exist('thr_all','var'),        thr_series = thr_all(:);
else,                                 thr_series = zeros(N,1);
end

% Interpolate if sizes differ
if numel(m_series) ~= N
    t_src = linspace(t_all(1), t_all(end), numel(m_series));
    m_series = interp1(t_src, m_series, t_all, 'linear','extrap');
end
if numel(thr_series) ~= N
    t_src = linspace(t_all(1), t_all(end), numel(thr_series));
    thr_series = interp1(t_src, thr_series, t_all, 'linear','extrap');
    thr_series(end)=0
end

prop_rem = max(m_series - m_dry, 0);                   % kg
prop0    = max(prop_rem(1), eps);
prop_pct = prop_rem / prop0;                           % 0..1
thr_pct  = min(max(thr_series,0), 1); 
thr_pct(end-1)=0;   % 0..1

% --- Animation parameters ---
realtime_factor = 55;
update_every    = 2;
z_zoom_start_km = 0.20;
zoom_frames     = 40;
tail_for_roi    = 300;

% --- Figure / axes (3D + 2 gauges on the right) ---
fig = figure('Name','3D Animation + Gauges');
% 3D area on the left
ax = axes('Parent',fig, 'Position',[0.06 0.08 0.62 0.84]); 
hold(ax,'on'); grid(ax,'on'); view(ax,135,25); set(ax,'XDir','reverse'); axis(ax,'vis3d');
xlabel(ax,'x (km)'); ylabel(ax,'y (km)'); zlabel(ax,'z (km)');

% PROPELLANT gauge (top right)
axFuel = axes('Parent',fig, 'Position',[0.74 0.60 0.22 0.30]);
axis(axFuel,[0 1 0 1]); axis(axFuel,'off'); title(axFuel,'Propellant','FontWeight','bold');
% Track
rectangle(axFuel,'Position',[0.35 0.08 0.30 0.84], ...
    'Curvature',0.15,'FaceColor',[0.90 0.90 0.90],'EdgeColor',[0.6 0.6 0.6]);
% Fill (green)
hFuelFill = rectangle(axFuel,'Position',[0.35 0.08 0.30 0.84*prop_pct(1)], ...
    'Curvature',0.15,'FaceColor',[0.20 0.80 0.30],'EdgeColor','none');
% Text
hFuelTxt  = text(axFuel,0.50,0.96, ...
    sprintf('%.1f t  (%.0f%%)', prop_rem(1)/1000, 100*prop_pct(1)), ...
    'HorizontalAlignment','center','FontWeight','bold');

% THROTTLE gauge (bottom right)
axThr = axes('Parent',fig, 'Position',[0.74 0.15 0.22 0.40]);
axis(axThr,[0 1 0 1]); axis(axThr,'off'); title(axThr,'Throttle','FontWeight','bold');
rectangle(axThr,'Position',[0.35 0.08 0.30 0.84], ...
    'Curvature',0.15,'FaceColor',[0.90 0.90 0.90],'EdgeColor',[0.6 0.6 0.6]);
hThrFill = rectangle(axThr,'Position',[0.35 0.08 0.30 0.84*thr_pct(1)], ...
    'Curvature',0.15,'FaceColor',[0.20 0.55 0.95],'EdgeColor','none');
hThrTxt  = text(axThr,0.50,0.96, ...
    sprintf('%.0f %%', 100*thr_pct(1)), ...
    'HorizontalAlignment','center','FontWeight','bold');
text(0.66,0.84,0.8,'Safe Throttle (90%)', ...
    'HorizontalAlignment','left','FontSize',8,'FontWeight','bold','Color',[0 0.8 0]);
text(0.05,0.88,0.8,'Limit Safe Throttle (95%)', ...
    'HorizontalAlignment','left','FontSize',8,'FontWeight','bold','Color',[1 0.5 0.2]);
text(0.66,0.925,1,'Ultimate Throttle (100%)', ...
    'FontWeight','bold','FontSize',8,'Color',[0.8 0.0 0.2]);
% Visual mark 90%
line(axThr,[0.35 0.655],[0.08+0.84*0.90 0.08+0.84*0.90],'Color',[0 0.8 0],'LineStyle','--','LineWidth',1);
line(axThr,[0.35 0.655],[0.08+0.887*0.90 0.08+0.887*0.90],'Color',[1 0.5 0],'LineStyle','--','LineWidth',1);
line(axThr,[0.35 0.655],[0.08+0.933*0.90 0.08+0.933*0.90],'Color',[0.8 0 0],'LineStyle','--','LineWidth',1);

%% ====== Warning "Not enough propellant" (red box) or Mission Success ======
% Small safeguards to recover params if needed
if ~exist('g0','var'),            g0 = 9.80665; end
if ~exist('params','var') || ~isstruct(params)
    params = struct; 
end
if ~isfield(params,'m0')     && exist('m0','var'),     params.m0     = m0;     end
if ~isfield(params,'m_dry')  && exist('m_dry','var'),  params.m_dry  = m_dry;  end
if ~isfield(params,'Isp')    && exist('Isp','var'),    params.Isp    = Isp;    end

% “Massepropellant” series (remaining prop) like you already use
Massepropellant = (params.m0 - params.m_dry) + (m_series - params.m0); % = m(t) - m_dry

% Choose a ΔV to check: DV_total else integral ΔV, else out.dV
if exist('DV_total','var')
    DV_chk = DV_total;
elseif exist('DV_int','var')
    DV_chk = DV_int;
elseif exist('out','var') && isfield(out,'dV')
    DV_chk = out.dV;
else
    DV_chk = 0;  % default
end

% Required propellant by rocket equation (your chosen formula)
m_prop_required  = params.m_dry * (exp(DV_chk/(params.Isp*g0)) - 1);
m_prop_available = params.m0 - params.m_dry;

notEnough = (m_prop_required > m_prop_available);
NoThrottle= (max(out.throttle(:))>0.98);
Enough=(m_prop_required < m_prop_available) && (vx_all(end)<0.3) && (vx_all(end)>-0.3) && NoThrottle==0;


% Create the box (normalized to figure coordinates)
warnStr = sprintf('NOT ENOUGH PROPELLANT\n Need %.1f t  >  Available %.1f t', ...
                  m_prop_required/1000, m_prop_available/1000);
warnStrThrottle = sprintf('Too High Throttle needed');
warnStr2 = sprintf('MISSION SUCCESS');
hWarn = annotation(fig,'textbox',[0.72 0.02 0.26 0.10], ...
    'String', warnStr, ...
    'Color',[0.85 0 0], 'FontWeight','bold','FontSize',18, ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'EdgeColor',[0.85 0 0], 'LineWidth',1.8, ...
    'BackgroundColor',[1.00 0.95 0.95], 'Margin',10, ...
    'Visible', ternary(notEnough,'on','off'));
hWarn = annotation(fig,'textbox',[0.72 0.02 0.26 0.10], ...
    'String', warnStrThrottle, ...
    'Color',[0.85 0 0], 'FontWeight','bold','FontSize',18, ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'EdgeColor',[0.85 0 0], 'LineWidth',1.8, ...
    'BackgroundColor',[1.00 0.95 0.95], 'Margin',10, ...
    'Visible', ternary(NoThrottle,'on','off'));
hWarn = annotation(fig,'textbox',[0.72 0.02 0.26 0.10], ...
    'String', warnStr2, ...
    'Color',[0 0.55 0], 'FontWeight','bold','FontSize',18, ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'EdgeColor',[0 0.85 0], 'LineWidth',1.8, ...
    'BackgroundColor',[0.95 1 0.95], 'Margin',10, ...
    'Visible', ternary(Enough,'on','off'));

% (small ternary helper)
function out = ternary(cond, A, B), if cond, out=A; else, out=B; end, end


% Initial global limits (robust, floor at z=0)
xlim(ax, safeLimitsKM(Xkm, 0.5));
ylim(ax, safeLimitsKM(Ykm, 0.5));
zlim(ax, safeZLimitsKM(Zkm, 0.2));    % <-- never < 0

% ================= Textured ground (z=0) =================
img = [];
for ext = ["png","jpg","jpeg"]
    fn = ['Lunar_surface.' + ext];
    if exist(fn,'file'), img = imread(fn); break; end
end
if isempty(img)
    cb  = uint8(255*checkerboard(16,4,4));
    img = cat(3,cb,cb,cb);
end
xg = safeLimitsKM(Xkm, 0.5);
yg = safeLimitsKM(Ykm, 0.5);
dx = diff(xg); dy = diff(yg);
xg = xg + 0.25*[-dx dx];  yg = yg + 0.25*[-dy dy];
[Xt,Yt] = meshgrid([xg(1) xg(2)], [yg(1) yg(2)]);
Zt = zeros(2);   % z=0 (km)
surface(ax, Xt, Yt, Zt, img, 'FaceColor','texturemap', 'EdgeColor','none', 'HandleVisibility','off');

% Animated traces per phase
cols = lines(3);
hAL(1) = animatedline('Parent',ax,'Color',cols(1,:),'LineWidth',2);
hAL(2) = animatedline('Parent',ax,'Color',cols(2,:),'LineWidth',2);
hAL(3) = animatedline('Parent',ax,'Color',cols(3,:),'LineWidth',2);

% Vehicle (marker)
[shipT, shipPatch, shipScaleKM] = loadShipModel(ax, 'Apollo Lunar Module.stl');  
lighting(ax,'gouraud'); material(shipPatch,'metal'); camlight(ax,'headlight');
vis_scale = 5000;   % visual magnification
z_eps     = 0.003;  % 3 m in km

% Zoom targets (computed on the end)
k_zoom   = find(Zkm <= z_zoom_start_km, 1, 'first');  if isempty(k_zoom), k_zoom = N; end
roi_idx  = max(1, N-tail_for_roi):N;
xlim_tgt = safeLimitsKM(Xkm(roi_idx), 0.2);
ylim_tgt = safeLimitsKM(Ykm(roi_idx), 0.2);
zlim_tgt = safeZLimitsKM(Zkm(roi_idx), 0.35);

did_zoom = false;  tic;
for k = 1:update_every:N
    ph = phase_all(k);  if ph<1 || ph>3, ph = 1; end
    addpoints(hAL(ph), Xkm(k), Ykm(k), Zkm(k));

    % Model pose: translation + yaw (from horizontal velocity)
    yaw = atan2(vy_all(k), vx_all(k));
    T  = makehgtform('translate', [Xkm(k) Ykm(k) Zkm(k)]);
    Rz = makehgtform('zrotate',   yaw);
    S  = makehgtform('scale',     shipScaleKM*vis_scale);
    shipT.Matrix = T * Rz * S;        

    % ------- Gauge updates (propellant & throttle) ------
    % (values are indexed by k; can also interpolate if needed)
    rFuel = min(max(prop_pct(k),0),1);
    set(hFuelFill,'Position',[0.35 0.08 0.30 0.84*rFuel]);
    set(hFuelTxt, 'String', sprintf('%.1f t  (%.0f%%)', prop_rem(k)/1000, 100*rFuel));

    rThr  = min(max(thr_pct(k),0),1);
    set(hThrFill, 'Position',[0.35 0.08 0.30 0.84*rThr]);
    set(hThrTxt,  'String', sprintf('%.0f %%', 100*rThr));
    % throttle color (green <=90%, orange >90%)
    if rThr <= 0.90
        set(hThrFill,'FaceColor',[0.20 0.55 0.95]);
    else
        set(hThrFill,'FaceColor',[0.95 0.45 0.20]);
    end

    title(ax, sprintf('t = %.1f s   z = %.0f m', t_all(k), z_all(k)));

    % Single zoom when passing below z_zoom_start_km
    if ~did_zoom && k >= k_zoom
        did_zoom = true;
        x0 = xlim(ax); y0 = ylim(ax); z0 = zlim(ax);
        for s = 1:zoom_frames
            a = s/zoom_frames;
            xlim(ax, (1-a)*x0 + a*xlim_tgt);
            ylim(ax, (1-a)*y0 + a*ylim_tgt);
            zNew = (1-a)*z0 + a*zlim_tgt; zNew(1) = max(zNew(1), 0);
            zlim(ax, zNew);
            drawnow limitrate;
        end
    end

    drawnow limitrate;

    % animation cadence ~ simulated time / realtime_factor
    if k>1
        dt_vis = (t_all(k) - t_all(k-update_every))/realtime_factor;
        if dt_vis > 0, pause(min(dt_vis, 0.05)); end
    end
end


% Start / Touchdown markers
plot3(ax, Xkm(1),Ykm(1),Zkm(1), 'ko','MarkerFaceColor','k');
plot3(ax, Xkm(end),Ykm(end),Zkm(end), 'ks','MarkerFaceColor','k');
text(Xkm(end),Ykm(end),Zkm(end), '  Touchdown','Parent',ax,'VerticalAlignment','top');
legend(ax,{'Braking','Approach','Vertical'},'Location','best');

%% ================== END 3D ANIMATION + GAUGES ===========================


% ---------- Utils ---------------------------------------------------------
function lim = safeLimitsKM(v, minSpan)
% Robust [lo hi] (km); adds margin and handles min==max
    v = v(isfinite(v));
    if isempty(v), lim = [-minSpan/2, minSpan/2]; return; end
    lo = min(v); hi = max(v); span = hi-lo;
    if span < max(minSpan, eps)
        mid = 0.5*(lo+hi); half = max(minSpan/2, 1e-6);
        lim = [mid-half, mid+half];
    else
        pad = 0.10*span; lim = [lo-pad, hi+pad];
    end
end

function lim = safeZLimitsKM(vz, minSpan)
% Like safeLimitsKM but ENFORCES lim(1)=0 (do not show z<0)
    vz = vz(isfinite(vz));
    if isempty(vz), lim = [0, minSpan]; return; end
    hi = max([0; vz]);    % includes 0
    span = hi - 0;
    if span < max(minSpan, eps)
        hi = max(minSpan, 1e-6);
    else
        hi = hi*(1+0.10); % 10% margin
    end
    lim = [0, hi];
end

function [hT, hP, scale_km] = loadShipModel(ax, filename)
% Load a 3D model (STEP/STP or STL) and place it in an hgtransform.
% Returns: hT (transform), hP (patch), scale_km (km scale consistent with your axes).
%
% Common assumption: STEP/STL files are in METERS -> convert to km.

    [~,~,ext] = fileparts(filename);
    if isempty(ext)
        cand = ["step","stp","stl"];
        for e = cand
            f = filename + "." + e;
            if exist(f,'file'), filename = f; ext = "."+e; break; end
        end
    end
    ext = lower(ext);
    hT = hgtransform('Parent', ax);
    scale_km = 1e-3;                         % m -> km

    switch ext
        case {'.step','.stp'}
            if ~exist('importGeometry','file')
                error(['importGeometry (PDE Toolbox) required for STEP/STP. ', ...
                       'Convert to STL if you do not have it.']);
            end
            mdl = createpde();
            importGeometry(mdl, filename);    % read STEP
            msh = generateMesh(mdl,'GeometricOrder','linear');  % surface mesh
            V = (msh.Nodes.' ) * 1e-3;        % m -> km
            F = (msh.Elements(1:3,:).').';
            hP = patch('Parent',hT, 'Faces',F, 'Vertices',V, ...
                       'FaceColor',[0.82 0.84 0.88], 'EdgeColor','none');
        case '.stl'
            if exist('stlread','file')
                TR = stlread(filename);       % R2018b+ or File Exchange
                V  = TR.Points * 1e-3;        % m -> km
                F  = TR.ConnectivityList;
                hP = patch('Parent',hT, 'Faces',F, 'Vertices',V, ...
                           'FaceColor',[0.82 0.84 0.88], 'EdgeColor','none');
            else
                error('stlread unavailable. Use a STEP (PDE) or add stlread.');
            end
        otherwise
            error('Unsupported extension: %s (use .step/.stp/.stl)', ext);
    end
end

function [ax_cmd,ay_cmd,az_cmd] = limit_tilt_and_slew(ax_cmd,ay_cmd,az_cmd, opts)
% Limite l'angle de tilt et la vitesse de rotation du vecteur a_cmd
%  - conserve (quasi) la norme -> trajectoire peu affectée
%  - agit seulement sur la direction
    arguments
        ax_cmd (1,1) double
        ay_cmd (1,1) double
        az_cmd (1,1) double
        opts.tilt_max_deg (1,1) double = 25      % tilt max vs vertical
        opts.dtheta_max_deg (1,1) double = 2     % rotation max par pas
    end

    persistent n_prev
    v = [ax_cmd ay_cmd az_cmd];
    amag = max(norm(v), 1e-9);
    n    = v/amag;                     % direction instantanée

    if isempty(n_prev), n_prev = n; end

    % --- 1) Rate limit sur la rotation (slew de direction)
    c = max(-1,min(1, dot(n_prev,n)));             % cos(angle)
    dtheta = acos(c);                               % [rad]
    dtheta_max = deg2rad(opts.dtheta_max_deg);
    if dtheta > dtheta_max
        % interpolation sphérique (slerp) vers n avec pas limité
        s = dtheta_max / dtheta;
        n = (sin((1-s)*dtheta)*n_prev + sin(s*dtheta)*n) / sin(dtheta);
        n = n / max(norm(n),1e-9);
    end

    % --- 2) Clamp du tilt vs vertical: tan(theta)=Th/Tz
    tilt_max = deg2rad(opts.tilt_max_deg);
    Th = hypot(n(1),n(2));
    tilt = atan2(Th, max(n(3),1e-9));
    if tilt > tilt_max
        scale = tan(tilt_max)/max(tan(tilt),1e-9);
        n(1:2) = n(1:2)*scale;
        n(3)   = sqrt(max(1 - sum(n(1:2).^2), 1e-9));
    end

    n_prev = n;
    ax_cmd = amag*n(1); ay_cmd = amag*n(2); az_cmd = amag*n(3);
end

% Série brute (avec sursauts)
raw = ang_vert;

% Paramètre de lissage (0 < alpha < 1, plus petit = plus lisse)
alpha = 0.9;

smooth = zeros(size(raw));
smooth(1) = raw(1);
for k = 2:length(raw)
    smooth(k) = alpha*smooth(k-1) + (1-alpha)*raw(k);
end

% Tracé comparatif
figure; hold on; grid on;
plot(t_all, out.tiltdeg, 'Color',[0.6 0.6 0.6]);   % courbe brute (gris)
plot(t_all, smooth, 'b','LineWidth',1.5); % courbe lissée
xlabel('t (s)'); ylabel('Angle to vertical (deg)');
legend('raw','smoothed');


%writematrix([x_all(:), y_all(:), z_all(:)], 'C:\Users\mperr\OneDrive - ISAE-SUPAERO\Bureau\SEEDS\Phase 3\Codes tests\positions.csv');

