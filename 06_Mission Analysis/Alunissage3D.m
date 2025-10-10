%% Lunar Powered Descent Guidance (2D) — braking, approach, vertical
% Implement the laws of guidance polynomials of the paper: Sostaric & Rea
% - Phase 1 (Braking)
% - Phase 2 (Approach)
% - Phase 3 (Vertical)

clear; clc;

%% ===== General Parameters =====
g_moon = 1.62;        % m/s^2
dt      = 0.2;        % pas d'intégration (s) — peut être augmenté si t_go est long
rng('default');

% Initial State (Start of PHASE 1 - Braking). Adjust for your case.
x0  = 100e3;%m (100,000)        % initial place of the target (m)
z0  = 15e3;%m  (15,000)
y0  = 5;%m 

% Initial Speeds
%vx0 = -1695 %m/s (this is the OG value, updated one from FDS is below)
vx0 = -1633.62;%m/s      % horizontal speed towards the target (m/s) (negative if x decreases)
vz0 = 0;    %m  
vy0 = 0;    %m/s      % vertical descending speed (m/s, negative towards the ground)

% Optional Limits (disable by putting [])
amax_x = [];          % m/s^2  (ex: 5)
amax_z = [];          % m/s^2  (ex: 5)
amax_y = [];

%% ===== Targets and durations for the phase (to adapt) =====
% Phase 1 (Braking) -> conditions relative to the start of the Apollo-like approach
x1_t  = 3e3;   % m (3,000)
y1_t  = 0e3;   % m (0) 
z1_t  = 1.2e3; % m (1,200)

vx1_t = -220;  % m/s  
vy1_t =  0;    % m/s
vz1_t = -45;   % m/s

Tgo1  = max(10, (x1_t - x0) / (0.5*(vx0 + vx1_t))); % s (durée allouée to the phase de freinage)

% Phase 2 (Approach) -> To bring the engine to the vertical of the site at 100 m
x2_t = 0;    % m  (to the vertical of the site)
y2_t = 0;    % m
z2_t = 100;  % m

vx2_t = 0;   % m/s (horizontal cancellation speed)
vy2_t = 0;   % m/s
vz2_t = -8;  % m/s (Target for the entry in the vertical phase)

Tgo2  = max(10, (x2_t - x1_t) / (0.5*(vx1_t + vx2_t)));

ax2_t = 0;           % m/s^2 (Target of the horizontal null acceleration)
ay2_t = 0;
az2_t = 1.2*g_moon;  % m/s^2 (Target of the engine's vertical acceleration ~1.2 g_Lune)

% Phase 3 (Vertical) -> Final Descent (constant)
az3_cmd = 1.2*g_moon; % m/s^2

vx_tol_touch    = 0.5;  % m/s (Horizontal tolerance at landing)
vy_tol_touch    = 0.5;  % m/s
vz_target_touch = -1.0; % m/s (Speed of the Landing Target ~ -1 m/s)

%% ===== Utilitaires =====
clip = @(u,umin,umax) max(min(u,umax),umin);

function [c1, c2] = braking_coeffs(r0, v0, rt, vt, tgo) %Function to calculate the braking coefficients
    c1 = (-2*(vt+2*v0))/tgo + 6*(rt - r0)/(tgo^2);
    c2 = ( 6*(vt+v0))/(tgo^2) - 12*(rt - r0)/(tgo^3);
end

function [c0,c1,c2] = approach_coeffs(r0,v0,rt,vt,at,tgo) %Function to calculate the approach coefficients
    c0 = at - 6*(vt+v0)/tgo + 12*(rt - r0)/(tgo^2);
    c1 = -6*at/tgo + 6*(5*vt + 3*v0)/(tgo^2) - 48*(rt - r0)/(tgo^3);
    c2 =  6*at/(tgo^2) - 12*(2*vt + v0)/(tgo^3) + 36*(rt - r0)/(tgo^4);
end
%% 

% -------- PHASE 1: Braking (a = c1 + c2 t, on applique a_cmd = c1) ------
% Current State
t = 0;

x = x0; 
y = y0; 
z = z0; 

vx = vx0; 
vy = vy0 ; 
vz = vz0;

% Index d'écriture
k = 0;
ax_cmd=-1;

% -------- PHASE 1: Braking --------
tgo = Tgo1;
while (z > z1_t) %Whilst the altitude is greater than the altitude at the end of Phase 1
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
    y  = y  + vy*dt;
    z  = z  + vz*dt;

    % Log (sans end+1)
    t = t + dt;  
    tgo = tgo - dt;  
    k = k + 1;

    t_all(k,1) = t;    
    x_all(k,1) = x; 
    y_all(k,1) = y ;  
    z_all(k,1) = z;

    vx_all(k,1)= vx;  
    vy_all(k,1)= vy;  
    vz_all(k,1)= vz;

    ax_all(k,1)= ax_cmd; 
    ay_all(k,1)= az_cmd; 
    az_all(k,1)= az_cmd;
    phase_all(k,1)=1;

    if z <= 0, break; end
end
fprintf('Phase 1 Iterations: %d\n', size(z_all,1));
%% 

% -------- PHASE 2: Approach (a = c0 + c1 t + c2 t^2, on applique a_cmd = c0) --------
tgo = Tgo2;
% Etat courant
while  (z > z2_t)
    % Coeffs séparés par axe
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

fprintf('Phase 2 Iterations: %d\n', size(z_all,1));
%% 

% -------- PHASE 3: Vertical + stratégie DIVERT ----------------------------
% (1) Paramètres inchangés
if ~exist('Kxy','var'),      Kxy = 0.8; end
if ~exist('ax_max','var'),   ax_max = 2.0; end
if ~exist('ay_max','var'),   ay_max = 2.0; end
if ~exist('daz_max','var'),  daz_max = 0.6*g_moon; end
if ~exist('az_cmd_filt','var'), az_cmd_filt = az3_cmd; end

v_td   = -0.1;                 % m/s Target at impact
az_nom = az3_cmd;              % m/s^2 (1.2 g_lune)
az_max = az_nom + daz_max;     % borne de poussée
a_up_max = max(az_max - g_moon, 0.1);  % capacité nette de freinage
alpha_f  = 0.6;                % lissage 0..1

% (2) Table DIVERT (depuis ton tableau)
divert_D   = [25   50    75   100   125   150];      % m
divert_dT  = [-0.9 -0.6  4.0   8.3   11.1  15.4];    % s
divert_dDV = [-0.9  0.6  8.2  15.5   20.3  27.4];    % m/s (info)

% Temps additionnel restant (décroît avec le temps)
dT_hold = 0;    % s
dDV_est = 0;    % m/s (optionnel, pour log)

while z > 0
    % ================== 1) Correcteurs horizontaux ==================
    if z > 50
        vx_lim = 10; vy_lim = 10;      % limites sous 50 m
        ax_cmd = (abs(vx)>vx_lim) * (-Kxy*(abs(vx)-vx_lim)*sign(vx));
        ay_cmd = (abs(vy)>vy_lim) * (-Kxy*(abs(vy)-vy_lim)*sign(vy));
    else
        band = 0.1;                    % bande morte ±0.1 m/s sous 10 m
        ax_cmd = (vx> band)*(-Kxy*(vx-band)) + (vx<-band)*(-Kxy*(vx+band));
        ay_cmd = (vy> band)*(-Kxy*(vy-band)) + (vy<-band)*(-Kxy*(vy+band));
    end
    ax_cmd = clip(ax_cmd,-ax_max,ax_max);
    ay_cmd = clip(ay_cmd,-ay_max,ay_max);

    % ================== 2) Gestion DIVERT (distance latérale) ========
    D = hypot(x,y);                                   % m
    dT_target = interp1(divert_D, divert_dT, D, 'linear','extrap');
    dDV_target = interp1(divert_D, divert_dDV, D, 'linear','extrap');

    % On ne “raccourcit” jamais la descente : on ne retient que ΔT>0
    dT_target = max(0, dT_target);
    % On mémorise le besoin de temps max observé et on le consomme peu à peu
    dT_hold = max(dT_hold, dT_target);
    dDV_est = max(dDV_est, max(0,dDV_target));       % info (non utilisé en boucle)

    % Estimation grossière du temps-restant (sert à normaliser l’allongement)
    tgo_est = z / max(0.5, -vz + 0.1);               % s, robuste si vz≈0
    slow_factor = 1 + min(0.8, dT_hold / max(5, tgo_est));   % ≤ +80% de τ

    % ================== 3) Profil vertical v_ref(z) ==================
    v_env = -sqrt(max(v_td^2 + 2*a_up_max*max(z,0), 0));   % enveloppe stoppable
    if     z > 200, v_cap = v_env;
    elseif z > 50,  v_cap = max(v_env, -2);   % limite de vitesse intermédiaire
    elseif z > 10,  v_cap = max(v_env, -1);
    else            v_cap = max(v_env, v_td); % converge vers -0.2 m/s
    end
    v_ref = v_cap;

    % Constante de temps de base (plus ferme en bas), puis “étirée” par DIVERT
    if     z > 200, tau = 12.0;
    elseif z > 100, tau = 6.0;
    elseif z > 50,  tau = 4.0;
    elseif z > 20,  tau = 2.0;
    elseif z > 10,  tau = 1.0;
    else            tau = 0.5;
    end
    tau = tau * slow_factor;     % <<< on achète du temps pour le DIVERT

    % Accélération nette désirée & poussée commandée
    a_net_des = (v_ref - vz)/tau;                   % m/s^2
    az_cmd_raw = g_moon + a_net_des;                % m->poussée
    az_cmd_raw = clip(az_cmd_raw, 0, az_max);
    az_cmd = alpha_f*az_cmd_filt + (1-alpha_f)*az_cmd_raw;
    az_cmd_filt = az_cmd;

    % ================== 4) Intégration + log =========================
    vx = vx + ax_cmd*dt;   vy = vy + ay_cmd*dt;
    vz = vz + (az_cmd - g_moon)*dt;
    x  = x  + vx*dt;       y  = y  + vy*dt;        z = z + vz*dt;

    t = t + dt;
    t_all(end+1,1)=t; x_all(end+1,1)=x; y_all(end+1,1)=y; z_all(end+1,1)=z;
    vx_all(end+1,1)=vx; vy_all(end+1,1)=vy; vz_all(end+1,1)=vz;
    ax_all(end+1,1)=ax_cmd; ay_all(end+1,1)=ay_cmd; az_all(end+1,1)=az_cmd;
    phase_all(end+1,1)=3;

    % On consomme le stock de temps “acheté”
    dT_hold = max(0, dT_hold - dt);
end

% Corriger dernier point au sol
if z_all(end) < 0, z_all(end) = 0; end

% (optionnel) trace console de la stratégie adoptée
fprintf('\n[DIVERT] D = %.1f m → ΔT Target ≈ %.1f s, ΔV add. ≈ %.1fm/s (tab)\n', ...
        D, dT_target, dDV_target);

%% ===== ΔV & G en 3D =====
idx = (phase_all > 0);
a_thrust = sqrt(ax_all(idx).^2 + ay_all(idx).^2 + az_all(idx).^2);
dv_step  = a_thrust * dt;
DV_total = sum(dv_step);
DV_P1 = sum(dv_step(phase_all(idx)==1));
DV_P2 = sum(dv_step(phase_all(idx)==2));
DV_P3 = sum(dv_step(phase_all(idx)==3));
fprintf('\nΔV total = %.1fm/s  (P1 = %.1f, P2 = %.1f, P3 = %.1f)\n', DV_total, DV_P1, DV_P2, DV_P3);

g0 = 9.80665;
G_prop = sqrt(ax_all.^2 + ay_all.^2 + az_all.^2) / g0;
[maxG, iG] = max(G_prop);
fprintf('\nMax proper G = %.2f g at t = %.1f s (phase %d)\n', maxG, t_all(iG), phase_all(iG));

% Max par phase
for p = 1:3
    idx = (phase_all == p);
    if any(idx)
        [mx, ii] = max(G_prop(idx));
        tmx = t_all(find(idx,1,'first') + ii - 1);
        fprintf('  Phase %d : max = %.2fg at t = %.1fs\n', p, mx, tmx);
    end
end

%stop;
%% ===== Plots =====
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

% 3) Accélérations commandées
subplot(2,2,3); hold on; grid on;
plot(t_all, ax_all, 'LineWidth',1.5);
plot(t_all, ay_all, 'LineWidth',1.5);
plot(t_all,az_all, 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('a_{cmd} (m/s^2)');
title('Acceleration Command');
legend('a_x','a_y','a_z (moteur)','Location','best');


subplot(2,2,4); hold on; grid on;
plot(t_all, vx_all, t_all, vy_all, t_all, vz_all, 'LineWidth',1.1);
yline(0,'--'); xlabel('t (s)'); ylabel('V (m/s)'); title('Speed'); legend('v_x','v_y','v_z');

fprintf('\n--- Touchdown ---\n');
fprintf('t=%.1f s | x=%.1f m | y=%.1f m | z=%.2f m | vx=%.2f | vy=%.2f | vz=%.2f m/s\n', ...
    t_all(end), x_all(end), y_all(end), z_all(end), vx_all(end), vy_all(end), vz_all(end));

%% ===== Angle de poussée vs surface / verticale =====
% Vecteur de poussée en coordonnées monde (tes a_cmd enregistrés)
Tnorm = sqrt(ax_all.^2 + ay_all.^2 + az_all.^2);

% Angle avec la verticale (surface normale = +z)
ang_vert = nan(size(Tnorm));                     % [deg]
mask = Tnorm > 1e-9;                             % évite division par 0
ang_vert(mask) = acosd( az_all(mask) ./ Tnorm(mask) );

% Angle par rapport à la surface (plan horizontal)
ang_surface = 90 - ang_vert;                     % [deg], 0° = à plat, 90° = vertical

% (option) azimut de la poussée dans le plan horizontal
azimut_deg = atan2d(ay_all, ax_all);             % [-180,180]

% Stats utiles
[maxTilt, iMax] = max(ang_surface);
fprintf('Tilt max vs surface = %.2f deg à t = %.1f s, z = %.1f m (phase %d)\n', ...
        maxTilt, t_all(iMax), z_all(iMax), phase_all(iMax));

% Sauvegarde pour réutilisation
ang_vert_all    = ang_vert;
ang_surface_all = ang_surface;
azimut_all      = azimut_deg;

% Tracés rapides
figure('Name','Orientation de poussée');
subplot(3,1,1); grid on; plot(t_all, ang_vert, 'LineWidth',1.2);
ylabel('Angle à la verticale (deg)'); title('Poussée vs verticale/surface/azimut');

subplot(3,1,2); grid on; plot(t_all, ang_surface, 'LineWidth',1.2);
ylabel('Angle vs surface (deg)');

subplot(3,1,3); grid on; plot(t_all, azimut_deg, 'LineWidth',1.2);
ylabel('Azimut (deg)'); xlabel('t (s)');

%% ================== ANIMATION 3D DE LA DESCENTE (avec sol luniare) ======
% Utilise t_all, x_all, y_all, z_all, phase_all

Xkm = x_all/1000;  Ykm = Xkm;  Zkm = z_all/1000;
N   = numel(t_all);

% --- Paramètres d'animation ---
realtime_factor = 6;
update_every    = 2;
z_zoom_start_km = 0.20;
zoom_frames     = 40;
tail_for_roi    = 300;

% --- Figure / axes ---
fig = figure('Name','Animation 3D'); ax = axes(fig); hold(ax,'on'); grid(ax,'on');
view(ax,135,25); set(ax,'XDir','reverse'); axis(ax,'vis3d');
xlabel(ax,'x (km)'); ylabel(ax,'y (km)'); zlabel(ax,'z (km)');

% Limites globales initiales (robustes, ET plancher z=0)
xlim(ax, safeLimitsKM(Xkm, 0.5));
ylim(ax, safeLimitsKM(Ykm, 0.5));
zlim(ax, safeZLimitsKM(Zkm, 0.2));    % <-- jamais < 0

% ================= Sol texturé (z=0) =================
% Essaie 'Lunar_surface.png' puis .jpg/.jpeg ; sinon damier de secours
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

dx = diff(xg);                % = xg(2)-xg(1)
dy = diff(yg);                % = yg(2)-yg(1)

xg = xg + 0.25*[-dx dx];      % marge 25%
yg = yg + 0.25*[-dy dy];
[Xt,Yt] = meshgrid([xg(1) xg(2)], [yg(1) yg(2)]);
Zt = zeros(2);   % z=0 (km)
hGround = surface(ax, Xt, Yt, Zt, img, ...
    'FaceColor','texturemap', 'EdgeColor','none', 'HandleVisibility','off');

% Traces animées par phase
cols = lines(3);
hAL(1) = animatedline('Parent',ax,'Color',cols(1,:),'LineWidth',2);
hAL(2) = animatedline('Parent',ax,'Color',cols(2,:),'LineWidth',2);
hAL(3) = animatedline('Parent',ax,'Color',cols(3,:),'LineWidth',2);

% Véhicule (marqueur)
[shipT, shipPatch, shipScaleKM] = loadShipModel(ax, 'Apollo Lunar Module.stl');  
lighting(ax,'gouraud'); material(shipPatch,'metal'); camlight(ax,'headlight');
vis_scale = 5000;        % <<< AGRANDISSEMENT VISUEL (x5000)
z_eps     = 0.003;     % 3 m en km : évite la fusion avec le sol (z-fighting)

% Cibles de zoom (calculées sur la fin)
k_zoom = find(Zkm <= z_zoom_start_km, 1, 'first');  if isempty(k_zoom), k_zoom = N; end
roi_idx   = max(1, N-tail_for_roi):N;
xlim_tgt  = safeLimitsKM(Xkm(roi_idx), 0.2);
ylim_tgt  = safeLimitsKM(Ykm(roi_idx), 0.2);
zlim_tgt  = safeZLimitsKM(Zkm(roi_idx), 0.35);   % <-- jamais < 0

did_zoom = false;  tic;
for k = 1:update_every:N
    ph = phase_all(k);  if ph<1 || ph>3, ph = 1; end

    addpoints(hAL(ph), Xkm(k), Ykm(k), Zkm(k));
    % Pose du modèle : translation + yaw (direction de la vitesse horizontale)
    yaw = atan2(vy_all(k), vx_all(k));           % orientation simple
    T  = makehgtform('translate', [Xkm(k) Ykm(k) Zkm(k)]);
    Rz = makehgtform('zrotate',   yaw);
    S  = makehgtform('scale',     shipScaleKM*vis_scale);  % m -> km (échelle du modèle)
    shipT.Matrix = T * Rz * S;        
    title(ax, sprintf('t = %.1f s   z = %.0f m', t_all(k), z_all(k)));

    % Zoom une seule fois quand on passe sous z_zoom_start_km
    if ~did_zoom && k >= k_zoom
        did_zoom = true;
        x0 = xlim(ax); y0 = ylim(ax); z0 = zlim(ax);
        for s = 1:zoom_frames
            a = s/zoom_frames;
            xlim(ax, (1-a)*x0 + a*xlim_tgt);
            ylim(ax, (1-a)*y0 + a*ylim_tgt);
            zNew = (1-a)*z0 + a*zlim_tgt;    % clamp plancher 0 au cas où
            zNew(1) = max(zNew(1), 0);
            zlim(ax, zNew);
            drawnow limitrate;
        end
    end

    drawnow limitrate;

    % cadence d’animation ~ temps simulé / realtime_factor
    if k>1
        dt_vis = (t_all(k) - t_all(k-update_every))/realtime_factor;
        if dt_vis > 0, pause(min(dt_vis, 0.05)); end
    end
end

% Marqueurs Start / Touchdown
plot3(ax, Xkm(1),Ykm(1),Zkm(1), 'ko','MarkerFaceColor','k');
plot3(ax, Xkm(end),Ykm(end),Zkm(end), 'ks','MarkerFaceColor','k');
text(Xkm(end),Ykm(end),Zkm(end), '  Touchdown','Parent',ax,'VerticalAlignment','top');
legend(ax,{'Braking','Approach','Vertical'},'Location','best');


%% ================== FIN ANIMATION 3D =====================================

% ---------- Utils ---------------------------------------------------------
function lim = safeLimitsKM(v, minSpan)
% [lo hi] (km) robuste ; ajoute marge et gère cas min==max
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
% Comme safeLimitsKM mais IMPOSE lim(1)=0 (on ne montre pas z<0)
    vz = vz(isfinite(vz));
    if isempty(vz), lim = [0, minSpan]; return; end
    hi = max([0; vz]);    % inclut 0
    span = hi - 0;
    if span < max(minSpan, eps)
        hi = max(minSpan, 1e-6);
    else
        hi = hi*(1+0.10); % 10% marge
    end
    lim = [0, hi];
end

function [hT, hP, scale_km] = loadShipModel(ax, filename)
% Charge un modèle 3D (STEP/STP ou STL) et le place dans un hgtransform.
% Retourne: hT (transform), hP (patch), scale_km (échelle km/cohérente avec tes axes).
%
% Hypothèse courante : les fichiers STEP/STL sont en METRES -> on convertit en km.

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
                error(['importGeometry (PDE Toolbox) requis pour STEP/STP. ', ...
                       'Convertis en STL si tu ne l’as pas.']);
            end
            mdl = createpde();
            importGeometry(mdl, filename);    % lit le STEP
            msh = generateMesh(mdl,'GeometricOrder','linear');  % maillage surface
            V = (msh.Nodes.' ) * 1e-3;        % m -> km
            F = (msh.Elements(1:3,:).').';
            hP = patch('Parent',hT, 'Faces',F, 'Vertices',V, ...
                       'FaceColor',[0.82 0.84 0.88], 'EdgeColor','none');
        case '.stl'
            if exist('stlread','file')
                TR = stlread(filename);       % R2018b+ ou File Exchange
                V  = TR.Points * 1e-3;        % m -> km
                F  = TR.ConnectivityList;
                hP = patch('Parent',hT, 'Faces',F, 'Vertices',V, ...
                           'FaceColor',[0.82 0.84 0.88], 'EdgeColor','none');
            else
                error('stlread indisponible. Utilise un STEP (PDE) ou ajoute stlread.');
            end
        otherwise
            error('Extension non supportée: %s (utilise .step/.stp/.stl)', ext);
    end
end
