%% Monte-Carlo (Tmax & Isp ±3σ) sur la trajectoire jerk-limitée (3 phases)
% Guidance = polynômes (Sostaric & Rea) + limite rotation vecteur a (slew)
% Évaluation propulsive via thrust_from_traj avec Tmax/Isp dispersés

clear; clc; rng('default');

%% ===== Constantes / pas =====
g_moon = 1.62;        % m/s^2
dt      = 0.2;        % s (guidance pas fin pour la traje)
v_td    = -0.1;       % m/s à l'atterrissage (phase 3)

% Propulsion / masses (nominales)
params = struct( ...
  'g_moon', g_moon, 'Isp', 450, 'g0', 9.80665, ...
  'Tmax', 100e3, 'm0', 18000, 'm_dry', 11000, 'tilt_max_deg', 25);

%% ===== Cibles FIXES =====
z0_fix  = 15e3;               % FIXED
vx0_fix = -1695;              % FIXED
x2_t = 0;  y2_t = 0;  z2_t = 100;
vx2_t = 0; vy2_t = 0; vz2_t = -8;
ax2_t = 0; ay2_t = 0; az2_t = 1.2*g_moon;
az3_cmd = 1.2*g_moon;

%% ===== Cas nominal (targets de la traje) =====
nom = struct();
nom.x0=289000; nom.y0=0; nom.z0=z0_fix;
nom.vx0=vx0_fix; nom.vy0=0; nom.vz0=0;
nom.x1_t=25000; nom.y1_t=0; nom.z1_t=7936.75104830166;
nom.vx1_t=-202; nom.vy1_t=0; nom.vz1_t=-49.09;
nom.x2_t=x2_t; nom.y2_t=y2_t; nom.z2_t=z2_t;
nom.vx2_t=vx2_t; nom.vy2_t=vy2_t; nom.vz2_t=vz2_t;
nom.ax2_t=ax2_t; nom.ay2_t=ay2_t; nom.az2_t=az2_t; nom.az3_cmd=az3_cmd;

[~, ~] = phase_times(nom);  %#ok<ASGLU>

%% ===== Dispersions (±3σ tronquées autour de 1.0) =====
THR_3SIG_PCT = 0.09;   % ±9% à 3σ sur Tmax
ISP_3SIG_PCT = 0.06;   % ±6% à 3σ sur Isp
SIG_T   = THR_3SIG_PCT/3;
SIG_ISP = ISP_3SIG_PCT/3;
draw_scale_pm3sig = @(sigma) min(1+3*sigma, max(1-3*sigma, 1 + sigma*randn()));

%% ===== Simulation nominale (Tscale=IspScale=1) =====
[sim_nom, out_nom] = simulate_jerk_guided_traj(nom, params, dt, v_td, 1.0, 1.0);

% Critères de succès (0.3%) + throttle <= 95% + pas de panne sèche
THR_MAX_OK = 0.95;
Sx=max([abs(nom.vx0),abs(nom.vx1_t),300]);
Sy=max([abs(nom.vy1_t),100]);
Sz=max([abs(nom.vz1_t),abs(nom.vz2_t),20]);
prop_rem_nom = max(out_nom.m - params.m_dry, 0);
fuelOK_nom   = all( prop_rem_nom(1:end-1) > 0 );

ok_nom = (abs(sim_nom.vx(end)) <= 0.003*Sx) && ...
         (abs(sim_nom.vy(end)) <= 0.003*Sy) && ...
         (abs(sim_nom.vz(end)-v_td) <= 0.003*Sz) && ...
         (max(out_nom.throttle) <= THR_MAX_OK) && ...
         fuelOK_nom;

fprintf('Nominal -> ok=%d | vx=%.3f vy=%.3f vz=%.3f | thrMax=%.1f%% | ΔV=%.1f\n', ...
    ok_nom, sim_nom.vx(end), sim_nom.vy(end), sim_nom.vz(end), 100*max(out_nom.throttle), out_nom.dV);

%% ===== Monte-Carlo =====
nMC = 100000;

% Variables dispersées autour du nominal (d'autres FIXES)
PCT = struct('x0',0.02,'x1',0.0,'z1',0.05,'vx1',0.10,'vz1',0.10); % variations modestes
ABS = struct('y0',5e3,'y1',0,'vy1',0);

cases = cell(nMC+1,1); sims = cell(nMC+1,1); outs = cell(nMC+1,1);
res = struct('success',false(nMC+1,1),'dv',nan(nMC+1,1),'thrmax',nan(nMC+1,1), ...
             'vx_td',nan(nMC+1,1),'vy_td',nan(nMC+1,1),'vz_td',nan(nMC+1,1), ...
             't',nan(nMC+1,1),'Tscale',nan(nMC+1,1),'IspScale',nan(nMC+1,1));

% Run #1 = nominal
cases{1}=nom; sims{1}=sim_nom; outs{1}=out_nom;
[ok1, dv1, thr1] = score_run(nom, sims{1}, outs{1}, v_td, THR_MAX_OK, params.m_dry);
res.success(1)=ok1; res.dv(1)=dv1; res.thrmax(1)=thr1;
res.vx_td(1)=sim_nom.vx(end); res.vy_td(1)=sim_nom.vy(end); res.vz_td(1)=sim_nom.vz(end); res.t(1)=sim_nom.t(end);
res.Tscale(1)=1.0; res.IspScale(1)=1.0;

for i=2:nMC+1
    s = nom;
    s.x0   = vary_pct(nom.x0,  PCT.x0);
    s.y0   = vary_abs(nom.y0, ABS.y0);
    s.z0   = z0_fix;                 % FIXED
    s.vx0  = vx0_fix;                % FIXED
    s.vy0  = nom.vy0; s.vz0 = nom.vz0;

    s.x1_t = vary_pct(nom.x1_t, PCT.x1);
    s.y1_t = nom.y1_t + (2*rand-1)*ABS.y1;
    s.z1_t = vary_pct(nom.z1_t, PCT.z1);

    s.vx1_t= vary_pct(nom.vx1_t, PCT.vx1);
    s.vy1_t= nom.vy1_t + (2*rand-1)*ABS.vy1;
    s.vz1_t= vary_pct(nom.vz1_t, PCT.vz1);

    % Dispersions propulsives
    Tscale_i   = draw_scale_pm3sig(SIG_T);
    IspScale_i = draw_scale_pm3sig(SIG_ISP);

    % Trajectoire jerk-limitée + évaluation poussée avec Tmax/Isp dispersés
    [sim, out] = simulate_jerk_guided_traj(s, params, dt, v_td, Tscale_i, IspScale_i);

    cases{i}=s; sims{i}=sim; outs{i}=out;
    [ok_i, dv_i, thr_i] = score_run(s, sim, out, v_td, THR_MAX_OK, params.m_dry);
    res.success(i)=ok_i; res.dv(i)=dv_i; res.thrmax(i)=thr_i;
    res.vx_td(i)=sim.vx(end); res.vy_td(i)=sim.vy(end); res.vz_td(i)=sim.vz(end); res.t(i)=sim.t(end);
    res.Tscale(i) = Tscale_i; res.IspScale(i) = IspScale_i;
end

okIdx_all = find(res.success);
fprintf('MC: %d/%d successes (%.1f%%)\n', numel(okIdx_all), numel(res.success), 100*numel(okIdx_all)/numel(res.success));

if ~isempty(okIdx_all)
    [~,k] = min(res.dv(okIdx_all)); best = okIdx_all(k);
else
    [~,best] = min(res.dv);
end

fprintf(['Best case: run %d | ΔV=%.1f | thrMax=%.1f%% | v_td=[%.3f %.3f %.3f] | t=%.1fs | ' ...
         'Tscale=%.3f | IspScale=%.3f\n'], ...
    best, res.dv(best), 100*res.thrmax(best), res.vx_td(best), res.vy_td(best), res.vz_td(best), ...
    res.t(best), res.Tscale(best), res.IspScale(best));

% ---- Petits plots sur le best case
bs = sims{best}; bo = outs{best};
figure('Name','Best case'); tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
nexttile; hold on; grid on;
C=lines(3);
for ph=1:3, idx=(bs.phase==ph); plot3(bs.x(idx)/1000,bs.y(idx)/1000,bs.z(idx)/1000,'.','Color',C(ph,:)); end
set(gca,'XDir','reverse'); xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)'); view(135,25); title('Trajectory');

nexttile; hold on; grid on; plot(bs.t, 100*bo.throttle);
if exist('yline','file'), yline(95,'--r'); end
xlabel('t (s)'); ylabel('Throttle (%)'); title('Throttle');

nexttile; hold on; grid on; plot(bs.t, max(bo.m-params.m_dry,0)/1000);
xlabel('t (s)'); ylabel('Propellant (t)'); title('Remaining prop');

nexttile; hold on; grid on; plot(bs.t, bs.vx, bs.t, bs.vy, bs.t, bs.vz);
if exist('yline','file'), yline(0,'--k'); end
legend('vx','vy','vz'); title('Velocities');

%% ===== Intervalles de confiance sur ΔV (succès uniquement) =====
okIdx = find(res.success);
if ~isempty(okIdx)
    dv_ok = res.dv(okIdx);
    q68   = prctile(dv_ok,[16 50 84]);
    q95   = prctile(dv_ok,[2.5 50 97.5]);
    q997  = prctile(dv_ok,[0.15 99.85]);

    fprintf('ΔV median = %.1f m/s\n', q95(2));
    fprintf('  68%% CI   = [%.1f, %.1f] m/s\n', q68(1),  q68(3));
    fprintf('  95%% CI   = [%.1f, %.1f] m/s\n', q95(1),  q95(3));
    fprintf('  99.7%% CI = [%.1f, %.1f] m/s\n', q997(1), q997(2));

    figure('Name','ΔV vs dispersions — 95% CI & normal bell');
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    dv_best     = res.dv(best);
    dv95_lo     = q95(1); dv95_med = q95(2); dv95_hi = q95(3);
    delta95_abs = max(0, dv95_hi - dv_best);
    delta95_pct = 100 * delta95_abs / max(dv_best, eps);

    nexttile; hold on; grid on;
    scatter(res.Tscale(okIdx), dv_ok, 8, 'filled');
    xlabel('T_{scale}'); ylabel('\DeltaV (m/s)'); title('\DeltaV vs T_{scale}');
    xl = xlim; patch([xl(1) xl(2) xl(2) xl(1)],[dv95_lo dv95_lo dv95_hi dv95_hi],[1 0 0],'FaceAlpha',0.08,'EdgeColor','none');
    if exist('yline','file'), yline(dv_best,'--k'); yline(dv95_med,':k'); end
    yl = ylim; text(xl(1)+0.02*(xl(2)-xl(1)), yl(2)-0.08*(yl(2)-yl(1)), ...
        sprintf('95%%: [%.1f, %.1f] — Best %.1f\n+%.1f m/s (%.1f%%) for 95%%', dv95_lo,dv95_hi,dv_best,delta95_abs,delta95_pct), ...
        'FontWeight','bold','BackgroundColor','w');

    nexttile; hold on; grid on;
    scatter(res.IspScale(okIdx), dv_ok, 8, 'filled');
    xlabel('Isp_{scale}'); ylabel('\DeltaV (m/s)'); title('\DeltaV vs Isp_{scale}');
    xl = xlim; patch([xl(1) xl(2) xl(2) xl(1)],[dv95_lo dv95_lo dv95_hi dv95_hi],[1 0 0],'FaceAlpha',0.08,'EdgeColor','none');
    if exist('yline','file'), yline(dv_best,'--k'); yline(dv95_med,':k'); end
    yl = ylim; text(xl(1)+0.02*(xl(2)-xl(1)), yl(2)-0.08*(yl(2)-yl(1)), ...
        sprintf('95%%: [%.1f, %.1f] — Best %.1f\n+%.1f m/s (%.1f%%) for 95%%', dv95_lo,dv95_hi,dv_best,delta95_abs,delta95_pct), ...
        'FontWeight','bold','BackgroundColor','w');

    fprintf(['\n[95%% Summary] add +%.1f m/s (%.1f%%) to best %.1f m/s -> %.1f m/s for 95%% coverage.\n'], ...
            delta95_abs, delta95_pct, dv_best, dv95_hi);
end

%% ===== Nuage de TOUTES les trajectoires =====
okPlot = find(res.success);
koPlot = find(~res.success);

MAXPLOT = 1000;
rng(0);
if numel(okPlot) + numel(koPlot) > MAXPLOT
    nOK = min(numel(okPlot), ceil(0.60*MAXPLOT));
    nKO = min(numel(koPlot), MAXPLOT - nOK);
    okPlot = okPlot(randperm(numel(okPlot), nOK));
    koPlot = koPlot(randperm(numel(koPlot), nKO));
    fprintf('[TrajPlot] Sub-sampling: %d successes, %d failures displayed.\n', nOK, nKO);
end

figure('Name','All trajectories — green=success, red=failure'); hold on; grid on;
xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
title('All trajectories (green=success, red=failure)');
set(gca,'XDir','reverse'); view(135,25);

colOK = [0.0 0.6 0.0]; colKO = [0.85 0.15 0.15];
for idx = okPlot(:)' , si = sims{idx}; plot3(si.x/1000, si.y/1000, si.z/1000, 'Color', colOK, 'LineWidth', 0.5); end
for idx = koPlot(:)' , si = sims{idx}; plot3(si.x/1000, si.y/1000, si.z/1000, 'Color', colKO, 'LineWidth', 0.5); end

td_ok = []; td_ko = [];
for idx = okPlot(:)' , si = sims{idx}; td_ok(end+1,:) = [si.x(end)/1000, si.y(end)/1000, si.z(end)/1000]; end %#ok<SAGROW>
for idx = koPlot(:)' , si = sims{idx}; td_ko(end+1,:) = [si.x(end)/1000, si.y(end)/1000, si.z(end)/1000]; end %#ok<SAGROW>
if ~isempty(td_ok),  plot3(td_ok(:,1), td_ok(:,2), td_ok(:,3), 'g.', 'MarkerSize', 8); end
if ~isempty(td_ko),  plot3(td_ko(:,1), td_ko(:,2), td_ko(:,3), 'r.', 'MarkerSize', 8); end

h1 = plot3(nan,nan,nan,'-','Color',colOK,'LineWidth',1.5);
h2 = plot3(nan,nan,nan,'-','Color',colKO,'LineWidth',1.5);
legend([h1 h2], {'Success','Fail'}, 'Location','best');

%% ===== Rapport BEST CASE =====
bestCase = cases{best};
[Tgo1_best, Tgo2_best] = phase_times(bestCase);
metrics = struct( ...
    'run',           best, ...
    'DeltaV_mps',    res.dv(best), ...
    'ThrottleMax',   res.thrmax(best), ...
    't_touch_s',     res.t(best), ...
    'vx_td_mps',     res.vx_td(best), ...
    'vy_td_mps',     res.vy_td(best), ...
    'vz_td_mps',     res.vz_td(best), ...
    'Tscale',        res.Tscale(best), ...
    'IspScale',      res.IspScale(best) );
fixed = struct( ...
    'dt_s',      dt, ...
    'g_moon',    params.g_moon, ...
    'v_td_mps',  v_td, ...
    'Propulsion', struct( ...
        'Isp_s',        params.Isp, ...
        'g0',           params.g0, ...
        'Tmax_N',       params.Tmax, ...
        'm0_kg',        params.m0, ...
        'm_dry_kg',     params.m_dry, ...
        'tilt_max_deg', params.tilt_max_deg), ...
    'Dispersion', struct( ...
        'THR_3SIG_PCT', THR_3SIG_PCT, ...
        'ISP_3SIG_PCT', ISP_3SIG_PCT));

fprintf('\n===== BEST CASE — Parameters used =====\n');
fprintf('PHASE 0 (initial)\n');
fprintf('  x0=%.3e, y0=%.3e, z0=%.3e\n', bestCase.x0, bestCase.y0, bestCase.z0);
fprintf('  vx0=%.3f, vy0=%.3f, vz0=%.3f\n', bestCase.vx0, bestCase.vy0, bestCase.vz0);
fprintf('PHASE 1 (braking targets)\n');
fprintf('  x1_t=%.3e, y1_t=%.3e, z1_t=%.3e | vx1_t=%.3f, vy1_t=%.3f, vz1_t=%.3f | Tgo1=%.2f s\n', ...
        bestCase.x1_t,bestCase.y1_t,bestCase.z1_t,bestCase.vx1_t,bestCase.vy1_t,bestCase.vz1_t,Tgo1_best);
fprintf('PHASE 2 (approach -> fixed targets)\n');
fprintf('  x2_t=%.3e, y2_t=%.3e, z2_t=%.3e | vx2_t=%.3f, vy2_t=%.3f, vz2_t=%.3f | Tgo2=%.2f s\n', ...
        bestCase.x2_t,bestCase.y2_t,bestCase.z2_t,bestCase.vx2_t,bestCase.vy2_t,bestCase.vz2_t,Tgo2_best);
fprintf('PHASE 3 (vertical)  az3_cmd=%.3f  v_td=%.3f\n', bestCase.az3_cmd, v_td);
fprintf('ΔV=%.1f m/s | throttle_max=%.1f %% | t_touch=%.1f s | Tscale=%.3f | IspScale=%.3f\n', ...
    metrics.DeltaV_mps, 100*metrics.ThrottleMax, metrics.t_touch_s, metrics.Tscale, metrics.IspScale);

% Sauvegarde optionnelle
best_report = struct('BestCase',bestCase,'Metrics',metrics,'Fixed',fixed);
json = jsonencode(best_report,'PrettyPrint',true);
fid = fopen('best_case_report.json','w'); fwrite(fid, json); fclose(fid);
save('best_case_report.mat','best_report');

%% ===================== FONCTIONS =====================

function [ok,dv,thrmax] = score_run(s, sim, out, v_td, THR_MAX_OK, m_dry)
    Sx = max([abs(s.vx0), abs(s.vx1_t), 300]);
    Sy = max([abs(s.vy1_t), 100]);
    Sz = max([abs(s.vz1_t), abs(s.vz2_t), 20]);
    prop_rem = max(out.m - m_dry, 0);
    prop_ok  = all( prop_rem(1:end-1) > 0 );
    ok = (abs(sim.vx(end)) <= 0.003*Sx) && ...
         (abs(sim.vy(end)) <= 0.003*Sy) && ...
         (abs(sim.vz(end)-v_td) <= 0.003*Sz) && ...
         (max(out.throttle) <= THR_MAX_OK) && prop_ok;
    dv     = out.dV;
    thrmax = max(out.throttle);
end

function y = vary_pct(x,p), y = x*(1 + p*(2*rand-1)); end
function y = vary_abs(x,a), y = x + (2*rand-1)*a; end

function [Tgo1,Tgo2] = phase_times(s)
    Tgo1 = max(10, (s.x1_t - s.x0) / (0.5*(s.vx0 + s.vx1_t))) + 5; % petit buffer
    Tgo2 = max(10, (s.x2_t - s.x1_t) / (0.5*(s.vx1_t + s.vx2_t))) + 5;
end

function [sim, out] = simulate_jerk_guided_traj(S, params, dt, v_td, Tscale, IspScale)
% Génère la trajectoire 3 phases (polynomial + jerk limiting) puis
% calcule la poussée requise (avec Tmax/Isp dispersés).
g_moon = params.g_moon;

% --- helpers poly
    function [c1,c2]=braking_coeffs(r0,v0,rt,vt,tgo)
        c1=(-2*(vt+2*v0))/tgo + 6*(rt-r0)/(tgo^2);
        c2=( 6*(vt+v0))/(tgo^2) - 12*(rt-r0)/(tgo^3);
    end
    function [c0,c1,c2]=approach_coeffs(r0,v0,rt,vt,at,tgo)
        c0 = at - 6*(vt+v0)/tgo + 12*(rt - r0)/(tgo^2);
        c1 = -6*at/tgo + 6*(5*vt + 3*v0)/(tgo^2) - 48*(rt - r0)/(tgo^3);
        c2 =  6*at/(tgo^2) - 12*(2*vt + v0)/(tgo^3) + 36*(rt - r0)/(tgo^4);
    end
    function [ax_cmd,ay_cmd,az_cmd,prev]=slew_limit(axd,ayd,azd,prev,rate_deg_s,dt)
        if any(isnan(prev)), ax_cmd=axd; ay_cmd=ayd; az_cmd=azd; prev=[ax_cmd ay_cmd az_cmd]; return; end
        v0=prev; v1=[axd ayd azd]; n0=norm(v0); n1=norm(v1);
        if n0<1e-9 || n1<1e-9, ax_cmd=axd; ay_cmd=ayd; az_cmd=azd; prev=[ax_cmd ay_cmd az_cmd]; return; end
        u0=v0/n0; u1=v1/n1; c=max(-1,min(1,dot(u0,u1))); ang=acosd(c);
        max_step = rate_deg_s*dt;
        if ang<=max_step, ax_cmd=axd; ay_cmd=ayd; az_cmd=azd;
        else
            s = max_step/ang; u = (1-s)*u0 + s*u1; u = u/max(norm(u),1e-12);
            v = u*n1; ax_cmd=v(1); ay_cmd=v(2); az_cmd=v(3);
        end
        prev=[ax_cmd ay_cmd az_cmd];
    end

% --- init
tilt_rate_max_deg = 5;  % limite rotation du vecteur accélération
Kxy=0.8; ax_max=2.0; ay_max=2.0; daz_max=0.6*g_moon;
az_nom=1.2*g_moon; az_max=az_nom+daz_max; a_up_max=max(az_max-g_moon,0.1);
alpha_f=0.6; az_cmd_filt=az_nom;

t=0; x=S.x0; y=S.y0; z=S.z0; vx=S.vx0; vy=S.vy0; vz=S.vz0;
t_all=t; x_all=x; y_all=y; z_all=z; vx_all=vx; vy_all=vy; vz_all=vz;
ax_all=0; ay_all=0; az_all=0; phase_all=0;

prev_acc = [NaN NaN NaN];

% ===== Phase 1 =====
Tgo1 = max(10, (S.x1_t - S.x0) / (0.5*(S.vx0 + S.vx1_t))) + 5; tgo=Tgo1;
while (z > S.z1_t)
    [c1x,~]=braking_coeffs(x,vx,S.x1_t,S.vx1_t,tgo);
    [c1y,~]=braking_coeffs(y,vy,S.y1_t,S.vy1_t,tgo);
    [c1z,~]=braking_coeffs(z,vz,S.z1_t,S.vz1_t,tgo);
    [ax_cmd,ay_cmd,az_cmd,prev_acc] = slew_limit(c1x,c1y,c1z,prev_acc,tilt_rate_max_deg,dt);

    vx=vx+ax_cmd*dt; vy=vy+ay_cmd*dt; vz=vz+(az_cmd - g_moon)*dt;
    x=x+vx*dt; y=y+vy*dt; z=z+vz*dt; t=t+dt; tgo=tgo-dt;

    t_all(end+1,1)=t; x_all(end+1,1)=x; y_all(end+1,1)=y; z_all(end+1,1)=z;
    vx_all(end+1,1)=vx; vy_all(end+1,1)=vy; vz_all(end+1,1)=vz;
    ax_all(end+1,1)=ax_cmd; ay_all(end+1,1)=ay_cmd; az_all(end+1,1)=az_cmd; phase_all(end+1,1)=1;
    if z<=0, break; end
end
if z<=0, finalize(); return; end

% ===== Phase 2 =====
Tgo2 = max(10, (S.x2_t - S.x1_t) / (0.5*(S.vx1_t + S.vx2_t))) + 5; tgo=Tgo2;
while (z > S.z2_t)
    [c0x,~,~]=approach_coeffs(x,vx,S.x2_t,S.vx2_t,S.ax2_t,tgo);
    [c0y,~,~]=approach_coeffs(y,vy,S.y2_t,S.vy2_t,S.ay2_t,tgo);
    [c0z,~,~]=approach_coeffs(z,vz,S.z2_t,S.vz2_t,S.az2_t,tgo);
    [ax_cmd,ay_cmd,az_cmd,prev_acc] = slew_limit(c0x,c0y,c0z,prev_acc,tilt_rate_max_deg,dt);

    vx=vx+ax_cmd*dt; vy=vy+ay_cmd*dt; vz=vz+(az_cmd - g_moon)*dt;
    x=x+vx*dt; y=y+vy*dt; z=z+vz*dt; t=t+dt; tgo=tgo-dt;

    t_all(end+1,1)=t; x_all(end+1,1)=x; y_all(end+1,1)=y; z_all(end+1,1)=z;
    vx_all(end+1,1)=vx; vy_all(end+1,1)=vy; vz_all(end+1,1)=vz;
    ax_all(end+1,1)=ax_cmd; ay_all(end+1,1)=ay_cmd; az_all(end+1,1)=az_cmd; phase_all(end+1,1)=2;
    if z<=0, break; end
end
if z<=0, finalize(); return; end

% ===== Phase 3 =====
while z>0
    if z>50
        vx_lim=10; vy_lim=10;
        ax_cmd_des=(abs(vx)>vx_lim)*(-Kxy*(abs(vx)-vx_lim)*sign(vx));
        ay_cmd_des=(abs(vy)>vy_lim)*(-Kxy*(abs(vy)-vy_lim)*sign(vy));
    else
        band=0.1;
        ax_cmd_des=(vx> band)*(-Kxy*(vx-band)) + (vx<-band)*(-Kxy*(vx+band));
        ay_cmd_des=(vy> band)*(-Kxy*(vy-band)) + (vy<-band)*(-Kxy*(vy+band));
    end
    ax_cmd_des=max(min(ax_cmd_des,ax_max),-ax_max);
    ay_cmd_des=max(min(ay_cmd_des,ay_max),-ay_max);

    v_env=-sqrt(max(v_td^2 + 2*a_up_max*max(z,0),0));
    if     z>200, v_ref=v_env; elseif z>50, v_ref=max(v_env,-2);
    elseif z>5, v_ref=max(v_env,-1); else, v_ref=max(v_env, v_td); end
    if     z>200, tau=12; elseif z>100, tau=6; elseif z>50, tau=4;
    elseif z>20, tau=2; elseif z>5, tau=1.5; else, tau=2; end

    a_net_des=(v_ref - vz)/tau;
    az_cmd_raw=g_moon + a_net_des; az_cmd_raw=max(min(az_cmd_raw,az_max),0);
    az_cmd=alpha_f*az_cmd_filt + (1-alpha_f)*az_cmd_raw; az_cmd_filt=az_cmd;

    [ax_cmd,ay_cmd,az_cmd,prev_acc] = slew_limit(ax_cmd_des,ay_cmd_des,az_cmd,prev_acc,tilt_rate_max_deg,dt);

    vx=vx+ax_cmd*dt; vy=vy+ay_cmd*dt; vz=vz+(az_cmd - g_moon)*dt;
    x=x+vx*dt; y=y+vy*dt; z=z+vz*dt; t=t+dt;

    t_all(end+1,1)=t; x_all(end+1,1)=x; y_all(end+1,1)=y; z_all(end+1,1)=z;
    vx_all(end+1,1)=vx; vy_all(end+1,1)=vy; vz_all(end+1,1)=vz;
    ax_all(end+1,1)=ax_cmd; ay_all(end+1,1)=ay_cmd; az_all(end+1,1)=az_cmd; phase_all(end+1,1)=3;
end
if z_all(end)<0, z_all(end)=0; end

finalize();

    function finalize()
        sim.t=t_all; sim.x=x_all; sim.y=y_all; sim.z=z_all;
        sim.vx=vx_all; sim.vy=vy_all; sim.vz=vz_all;
        sim.ax=ax_all; sim.ay=ay_all; sim.az=az_all; sim.phase=phase_all;
        sim.dv_total = trapz(sim.t, sqrt(sim.ax.^2+sim.ay.^2+sim.az.^2));

        % ==== Évaluation propulsive avec dispersions ====
        p_eval = params;
        p_eval.Tmax = params.Tmax * Tscale;
        p_eval.Isp  = params.Isp  * IspScale;
        V = [sim.vx(:), sim.vy(:), sim.vz(:)];
        out = thrust_from_traj(sim.t(:), V, p_eval);
        out.throttle(end)=0; % évite un spike final
        out.dV = sim.dv_total;  % cohérent avec intégrale |a_cmd|
    end
end

function out = thrust_from_traj(t, V, params)
% Requiert accel = dV/dt + [0 0 g]; calcule T, throttle, m(t), ΔV, etc.
    if size(V,2)==2, V = [V, zeros(size(V,1),1)]; end
    g    = getfielddef(params,'g_moon',1.62);
    Isp  = getfielddef(params,'Isp',450);
    g0   = getfielddef(params,'g0',9.80665);
    Tmax = getfielddef(params,'Tmax',1e5);
    m0   = params.m0;
    m_dry= getfielddef(params,'m_dry',0);
    tiltMaxDeg = getfielddef(params,'tilt_max_deg',inf);

    vx = V(:,1); vy = V(:,2); vz = V(:,3);
    ax = gradient(vx, t);
    ay = gradient(vy, t);
    az = gradient(vz, t);

    a_cmd = [ax, ay, az + g];  
    a_cmd(:,3) = max(a_cmd(:,3), 0);     % pas de poussée vers le bas
    amag  = sqrt(sum(a_cmd.^2,2));

    kappa = amag./(Isp*g0); int_kappa = cumtrapz(t, kappa);
    m = m0 .* exp(-int_kappa); if m_dry>0, m = max(m, m_dry); end

    T_vec = a_cmd .* m;  T_mag = m .* amag;
    throttle = T_mag / Tmax;

    Tx=T_vec(:,1); Ty=T_vec(:,2); Tz=T_vec(:,3);
    Th = hypot(Tx,Ty);
    tiltdeg = atan2d(Th, max(Tz,0));
    overTh  = throttle > 1;
    overTilt= tiltdeg  > tiltMaxDeg;
    vertNeg = Tz < 0;

    dV = trapz(t, amag);
    m_prop = m0 - m(end);

    out.T_vec=T_vec; out.T_mag=T_mag; out.throttle=throttle; out.m=m;
    out.dV=dV; out.m_prop=m_prop; out.a_cmd=a_cmd; out.tiltdeg=tiltdeg;
    out.flags = struct('sat_throttle',overTh,'vert_negative',vertNeg,'over_tilt',overTilt);
end

function v = getfielddef(s, name, def)
    if isfield(s,name) && ~isempty(s.(name)), v = s.(name); else, v = def; end
end
% ==== Tilt angle plot for best case ====
bs = sims{best};
bo = outs{best};

figure('Name','Best Case Thrust Tilt');
hold on; grid on;

plot(bs.t, bo.tiltdeg, 'Color',[0.2 0.6 0.9], 'LineWidth',1.5); % courbe principale


xlabel('t (s)'); ylabel('Tilt vs vertical (deg)');
title('Thrust vector tilt angle (best case)');
