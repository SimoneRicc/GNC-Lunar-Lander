%% Monte-Carlo around the nominal case (uses your integrator as-is)
clear; rng('default');

%% ===== General parameters =====
g_moon = 1.62;          % m/s^2
dt      = 0.5;          % s
v_td    = -0.1;         % target m/s in phase 3

% Propulsion/mass
params = struct( ...
  'g_moon', g_moon, 'Isp', 450, 'g0', 9.80665, ...
  'Tmax', 100e3, 'm0', 17500, 'm_dry', 11000, 'tilt_max_deg', 150);

%% ===== Requested FIXES =====
z0_fix  = 15e3;            % FIXED
vx0_fix = -1695;           % FIXED
x2_t = 0;  y2_t = 0;  z2_t = 100;    % FIXED for phase 2
vx2_t = 0; vy2_t = 0; vz2_t = -8;    % FIXED for phase 2
ax2_t = 0; ay2_t = 0; az2_t = 1.2*g_moon;
az3_cmd = 1.2*g_moon;

%% ===== Nominal case =====
nom = struct();
nom.x0=300e3; nom.y0=0; nom.z0=z0_fix;
nom.vx0=vx0_fix; nom.vy0=0; nom.vz0=0;
nom.x1_t=25e3; nom.y1_t=0; nom.z1_t=5e3;
nom.vx1_t=-220; nom.vy1_t=0; nom.vz1_t=-45;
nom.x2_t=x2_t; nom.y2_t=y2_t; nom.z2_t=z2_t;
nom.vx2_t=vx2_t; nom.vy2_t=vy2_t; nom.vz2_t=vz2_t;
nom.ax2_t=ax2_t; nom.ay2_t=ay2_t; nom.az2_t=az2_t; nom.az3_cmd=az3_cmd;

[Tgo1_nom, Tgo2_nom] = phase_times(nom);
sim_nom = simulate_case(nom, params, dt, v_td);
out_nom = thrust_from_traj(sim_nom.t, [sim_nom.vx sim_nom.vy sim_nom.vz], params);

% === success criteria (0.3%) ===
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


%% ===== Monte-Carlo (around nominal) =====
nMC = 10000;

% Variables TO VARY (only z0, x2_t, y2_t, z2_t, vx0 remain FIXED)
% Percent variations are relative; absolute variations are additive.
PCT = struct('x0',0.2,'x1',0.0,'z1',1,'vx1',0.5,'vz1',0.3);
ABS = struct('y0',0,'y1',0e3,'vy1',000);   % NEW: y0 varies ±50 km; y1/vy1 kept 0 here

cases = cell(nMC+1,1); sims = cell(nMC+1,1); outs = cell(nMC+1,1);
res = struct('success',false(nMC+1,1),'dv',nan(nMC+1,1),'thrmax',nan(nMC+1,1), ...
             'vx_td',nan(nMC+1,1),'vy_td',nan(nMC+1,1),'vz_td',nan(nMC+1,1), 't',nan(nMC+1,1));

% Insert nominal as #1 (guarantees at least one physically “working” case)
cases{1}=nom; sims{1}=sim_nom; outs{1}=out_nom;
[ok1, dv1, thr1] = score_run(nom, sims{1}, outs{1}, v_td, THR_MAX_OK, params.m_dry);
res.success(1)=ok1; res.dv(1)=dv1; res.thrmax(1)=thr1;
res.vx_td(1)=sim_nom.vx(end); res.vy_td(1)=sim_nom.vy(end); res.vz_td(1)=sim_nom.vz(end); res.t(1)=sim_nom.t(end);

for i=2:nMC+1
    s = nom;                                   % base case
    s.x0   = vary_pct(nom.x0,  PCT.x0);
    s.y0   = vary_abs(nom.y0, ABS.y0);         % <<< vary initial lateral offset y0
    s.z0   = z0_fix;                            % FIXED
    s.vx0  = vx0_fix;                           % FIXED
    s.vy0  = nom.vy0; s.vz0 = nom.vz0;

    s.x1_t = vary_pct(nom.x1_t, PCT.x1);
    s.y1_t = nom.y1_t + (2*rand-1)*ABS.y1;
    s.z1_t = vary_pct(nom.z1_t, PCT.z1);

    s.vx1_t= vary_pct(nom.vx1_t, PCT.vx1);
    s.vy1_t= nom.vy1_t + (2*rand-1)*ABS.vy1;
    s.vz1_t= vary_pct(nom.vz1_t, PCT.vz1);

    % Phase 2 and 3 remain FIXED (x2_t,y2_t,z2_t,vx2_t,vy2_t,vz2_t,az2_t,az3_cmd)

    sim = simulate_case(s, params, dt, v_td);
    out = thrust_from_traj(sim.t, [sim.vx sim.vy sim.vz], params);

    cases{i}=s; sims{i}=sim; outs{i}=out;
    [ok_i, dv_i, thr_i] = score_run(s, sim, out, v_td, THR_MAX_OK, params.m_dry);
    res.success(i)=ok_i; res.dv(i)=dv_i; res.thrmax(i)=thr_i;
    res.vx_td(i)=sim.vx(end); res.vy_td(i)=sim.vy(end); res.vz_td(i)=sim.vz(end); res.t(i)=sim.t(end);
end

okIdx = find(res.success);
fprintf('MC: %d/%d successes (%.1f%%)\n', numel(okIdx), numel(res.success), 100*numel(okIdx)/numel(res.success));

if ~isempty(okIdx)
    [~,k] = min(res.dv(okIdx)); best = okIdx(k);
else
    [~,best] = min(res.dv);   % no success: “closest” ΔV
end

fprintf('Best case: run %d | ΔV=%.1f | thrMax=%.1f%% | v_td=[%.3f %.3f %.3f] | t=%.1fs\n', ...
    best, res.dv(best), 100*res.thrmax(best), res.vx_td(best), res.vy_td(best), res.vz_td(best), res.t(best));

% ---- Small plots on the best case
bs = sims{best}; bo = outs{best};
figure('Name','Best case'); tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
nexttile; hold on; grid on;
C=lines(3);
for ph=1:3, idx=(bs.phase==ph); plot3(bs.x(idx)/1000,bs.y(idx)/1000,bs.z(idx)/1000,'.','Color',C(ph,:)); end
set(gca,'XDir','reverse'); xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)'); view(135,25); title('Trajectory');
nexttile; hold on; grid on; plot(bs.t, 100*bo.throttle); yline(95,'--r'); xlabel('t'); ylabel('Throttle (%)'); title('Throttle');
nexttile; hold on; grid on; plot(bs.t, max(bo.m-params.m_dry,0)/1000); xlabel('t'); ylabel('Propellant (t)'); title('Remaining prop');
nexttile; hold on; grid on; plot(bs.t, bs.vx, bs.t, bs.vy, bs.t, bs.vz); yline(0,'--k'); legend('vx','vy','vz'); title('Velocities');

%% ===================== FUNCTIONS =====================

function [ok,dv,thrmax] = score_run(s, sim, out, v_td, THR_MAX_OK, m_dry)
    % scaling 0.3 %
    Sx = max([abs(s.vx0), abs(s.vx1_t), 300]);
    Sy = max([abs(s.vy1_t), 100]);
    Sz = max([abs(s.vz1_t), abs(s.vz2_t), 20]);

    % remaining propellant > 0 before touchdown
    prop_rem = max(out.m - m_dry, 0);
    prop_ok  = all( prop_rem(1:end-1) > 0 );

    ok = (abs(sim.vx(end)) <= 0.003*Sx) && ...
         (abs(sim.vy(end)) <= 0.003*Sy) && ...
         (abs(sim.vz(end)-v_td) <= 0.003*Sz) && ...
         (max(out.throttle) <= THR_MAX_OK) && ...
         prop_ok;

    dv     = out.dV;
    thrmax = max(out.throttle);
end

function y = vary_pct(x,p), y = x*(1 + p*(2*rand-1)); end
function y = vary_abs(x,a), y = x + (2*rand-1)*a; end   % <<< helper for absolute ±a variation

function [Tgo1,Tgo2] = phase_times(s)
    Tgo1 = max(10, (s.x1_t - s.x0) / (0.5*(s.vx0 + s.vx1_t)));
    Tgo2 = max(10, (s.x2_t - s.x1_t) / (0.5*(s.vx1_t + s.vx2_t)));
end

function sim = simulate_case(S, params, dt, v_td)
% === Your original integrator, unchanged in spirit ===
g_moon=params.g_moon; clip=@(u,umin,umax) max(min(u,umax),umin);
% ---- coefficients
    function [c1,c2]=braking_coeffs(r0,v0,rt,vt,tgo)
        c1=(-2*(vt+2*v0))/tgo + 6*(rt-r0)/(tgo^2);
        c2=( 6*(vt+v0))/(tgo^2) - 12*(rt-r0)/(tgo^3);
    end
    function [c0,c1,c2]=approach_coeffs(r0,v0,rt,vt,at,tgo)
        c0 = at - 6*(vt+v0)/tgo + 12*(rt - r0)/(tgo^2);
        c1 = -6*at/tgo + 6*(5*vt + 3*v0)/(tgo^2) - 48*(rt - r0)/(tgo^3);
        c2 =  6*at/(tgo^2) - 12*(2*vt + v0)/(tgo^3) + 36*(rt - r0)/(tgo^4);
    end

% ---- initialization
t=0; x=S.x0; y=S.y0; z=S.z0; vx=S.vx0; vy=S.vy0; vz=S.vz0;
t_all=t; x_all=x; y_all=y; z_all=z; vx_all=vx; vy_all=vy; vz_all=vz;
ax_all=0; ay_all=0; az_all=0; phase_all=0;

% ---- PHASE 1
Tgo1 = max(10, (S.x1_t - S.x0) / (0.5*(S.vx0 + S.vx1_t)));
tgo=Tgo1;
while (z > S.z1_t)
    [c1x,~]=braking_coeffs(x,vx,S.x1_t,S.vx1_t,tgo);
    [c1y,~]=braking_coeffs(y,vy,S.y1_t,S.vy1_t,tgo);
    [c1z,~]=braking_coeffs(z,vz,S.z1_t,S.vz1_t,tgo);
    ax_cmd=c1x; ay_cmd=c1y; az_cmd=c1z;
    vx=vx+ax_cmd*dt; vy=vy+ay_cmd*dt; vz=vz+(az_cmd-g_moon)*dt;
    x=x+vx*dt; y=y+vy*dt; z=z+vz*dt;
    t=t+dt; tgo=tgo-dt;
    t_all(end+1,1)=t; x_all(end+1,1)=x; y_all(end+1,1)=y; z_all(end+1,1)=z;
    vx_all(end+1,1)=vx; vy_all(end+1,1)=vy; vz_all(end+1,1)=vz;
    ax_all(end+1,1)=ax_cmd; ay_all(end+1,1)=ay_cmd; az_all(end+1,1)=az_cmd; phase_all(end+1,1)=1;
    if z<=0, break; end
end
if z<=0, finalize(); return; end

% ---- PHASE 2
Tgo2 = max(10, (S.x2_t - S.x1_t) / (0.5*(S.vx1_t + S.vx2_t)));
tgo=Tgo2;
while (z > S.z2_t)
    [c0x,~,~]=approach_coeffs(x,vx,S.x2_t,S.vx2_t,S.ax2_t,tgo);
    [c0y,~,~]=approach_coeffs(y,vy,S.y2_t,S.vy2_t,S.ay2_t,tgo);
    [c0z,~,~]=approach_coeffs(z,vz,S.z2_t,S.vz2_t,S.az2_t,tgo);
    ax_cmd=c0x; ay_cmd=c0y; az_cmd=c0z;
    vx=vx+ax_cmd*dt; vy=vy+ay_cmd*dt; vz=vz+(az_cmd-g_moon)*dt;
    x=x+vx*dt; y=y+vy*dt; z=z+vz*dt;
    t=t+dt; tgo=tgo-dt;
    t_all(end+1,1)=t; x_all(end+1,1)=x; y_all(end+1,1)=y; z_all(end+1,1)=z;
    vx_all(end+1,1)=vx; vy_all(end+1,1)=vy; vz_all(end+1,1)=vz;
    ax_all(end+1,1)=ax_cmd; ay_all(end+1,1)=ay_cmd; az_all(end+1,1)=az_cmd; phase_all(end+1,1)=2;
    if z<=0, break; end
end
if z<=0, finalize(); return; end

% ---- PHASE 3 (vertical + XY correctors)
Kxy=0.8; ax_max=2.0; ay_max=2.0; daz_max=0.6*g_moon;
az_nom=1.2*g_moon; az_max=az_nom+daz_max; a_up_max=max(az_max-g_moon,0.1);
alpha_f=0.6; az_cmd_filt=az_nom;
while z>0
    if z>50
        vx_lim=10; vy_lim=10;
        ax_cmd=(abs(vx)>vx_lim)*(-Kxy*(abs(vx)-vx_lim)*sign(vx));
        ay_cmd=(abs(vy)>vy_lim)*(-Kxy*(abs(vy)-vy_lim)*sign(vy));
    else
        band=0.1;
        ax_cmd=(vx> band)*(-Kxy*(vx-band)) + (vx<-band)*(-Kxy*(vx+band));
        ay_cmd=(vy> band)*(-Kxy*(vy-band)) + (vy<-band)*(-Kxy*(vy+band));
    end
    ax_cmd=max(min(ax_cmd,ax_max),-ax_max);
    ay_cmd=max(min(ay_cmd,ay_max),-ay_max);

    v_env=-sqrt(max(v_td^2 + 2*a_up_max*max(z,0),0));
    if     z>200, v_ref=v_env;
    elseif z>50,  v_ref=max(v_env,-2);
    elseif z>5,   v_ref=max(v_env,-1);
    else          v_ref=max(v_env, v_td);
    end
    if     z>200, tau=12; elseif z>100, tau=6; elseif z>50, tau=4;
    elseif z>20,  tau=2;  elseif z>5,  tau=1.5; else, tau=2; end

    a_net_des=(v_ref - vz)/tau;
    az_cmd_raw=g_moon + a_net_des; az_cmd_raw=max(min(az_cmd_raw,az_max),0);
    az_cmd=alpha_f*az_cmd_filt + (1-alpha_f)*az_cmd_raw; az_cmd_filt=az_cmd;

    vx=vx+ax_cmd*dt; vy=vy+ay_cmd*dt; vz=vz+(az_cmd-g_moon)*dt;
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
        amag=sqrt(sim.ax.^2 + sim.ay.^2 + sim.az.^2);
        sim.dv_total = trapz(sim.t, amag);
    end
end


function out = thrust_from_traj(t, V, params)
    if size(V,2)==2, V=[V, zeros(size(V,1),1)]; end
    g=params.g_moon; Isp=params.Isp; g0=params.g0; Tmax=params.Tmax; m0=params.m0; m_dry=params.m_dry;
    vx=V(:,1); vy=V(:,2); vz=V(:,3);
    ax=gradient(vx,t); ay=gradient(vy,t); az=gradient(vz,t);
    a_cmd=[ax, ay, az + g]; a_cmd(:,3)=max(a_cmd(:,3),0);
    amag=sqrt(sum(a_cmd.^2,2));
    m = m0 .* exp(-cumtrapz(t, amag./(Isp*g0))); if m_dry>0, m=max(m,m_dry); end
    T_mag = m .* amag; throttle = T_mag / Tmax;
    T_vec = a_cmd.*m; Tx=T_vec(:,1); Ty=T_vec(:,2); Tz=T_vec(:,3);
    tiltdeg = atan2d(hypot(Tx,Ty), max(Tz,0));
    out.T_vec=T_vec; out.T_mag=T_mag; out.throttle=throttle; out.m=m;
    out.dV=trapz(t,amag); out.m_prop=m0-m(end); out.a_cmd=a_cmd; out.tiltdeg=tiltdeg;
    out.flags.sat_throttle = throttle>1; out.flags.vert_negative = Tz<0; out.flags.over_tilt = tiltdeg>params.tilt_max_deg;
end

% small utility
function v = getfielddef(s,name,def), if isfield(s,name)&&~isempty(s.(name)), v=s.(name); else, v=def; end, end

%% ===== Full report of BEST CASE parameters =====
bestCase = cases{best};              % all parameters drawn for this run
[Tgo1_best, Tgo2_best] = phase_times(bestCase);

% Associated metrics
metrics = struct( ...
    'run',           best, ...
    'DeltaV_mps',    res.dv(best), ...
    'ThrottleMax',   res.thrmax(best), ...      % fraction 0..1
    't_touch_s',     res.t(best), ...
    'vx_td_mps',     res.vx_td(best), ...
    'vy_td_mps',     res.vy_td(best), ...
    'vz_td_mps',     res.vz_td(best) );

% Fixed context (propulsion, integration)
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
        'tilt_max_deg', params.tilt_max_deg));

% Readable output
fprintf('\n===== BEST CASE — Parameters used =====\n');
fprintf('PHASE 0 (initial state)\n');
fprintf('  x0=%.3e m, y0=%.3e m, z0=%.3e m\n', bestCase.x0, bestCase.y0, bestCase.z0);
fprintf('  vx0=%.3f m/s, vy0=%.3f m/s, vz0=%.3f m/s\n', bestCase.vx0, bestCase.vy0, bestCase.vz0);

fprintf('PHASE 1 (braking -> targets)\n');
fprintf('  x1_t=%.3e m, y1_t=%.3e m, z1_t=%.3e m\n', bestCase.x1_t, bestCase.y1_t, bestCase.z1_t);
fprintf('  vx1_t=%.3f m/s, vy1_t=%.3f m/s, vz1_t=%.3f m/s\n', bestCase.vx1_t, bestCase.vy1_t, bestCase.vz1_t);
fprintf('  Tgo1=%.2f s\n', Tgo1_best);

fprintf('PHASE 2 (approach -> FIXED targets)\n');
fprintf('  x2_t=%.3e m, y2_t=%.3e m, z2_t=%.3e m\n', bestCase.x2_t, bestCase.y2_t, bestCase.z2_t);
fprintf('  vx2_t=%.3f m/s, vy2_t=%.3f m/s, vz2_t=%.3f m/s\n', bestCase.vx2_t, bestCase.vy2_t, bestCase.vz2_t);
fprintf('  ax2_t=%.3f m/s^2, ay2_t=%.3f m/s^2, az2_t=%.3f m/s^2\n', bestCase.ax2_t, bestCase.ay2_t, bestCase.az2_t);
fprintf('  Tgo2=%.2f s\n', Tgo2_best);

fprintf('PHASE 3 (vertical)\n');
fprintf('  az3_cmd=%.3f m/s^2, v_td=%.3f m/s\n', bestCase.az3_cmd, v_td);

fprintf('\nMETRICS (best run #%d)\n', metrics.run);
fprintf('  ΔV=%.1f m/s | throttle_max=%.1f %% | t_touch=%.1f s\n', ...
    metrics.DeltaV_mps, 100*metrics.ThrottleMax, metrics.t_touch_s);
fprintf('  Touchdown v = [vx=%.3f, vy=%.3f, vz=%.3f] m/s\n', ...
    metrics.vx_td_mps, metrics.vy_td_mps, metrics.vz_td_mps);

fprintf('\nCONSTANTS & PROPULSION\n');
fprintf('  dt=%.3f s, g_moon=%.3f m/s^2, v_td=%.3f m/s\n', fixed.dt_s, fixed.g_moon, fixed.v_td_mps);
fprintf('  Isp=%.0f s, g0=%.5f, Tmax=%.0f N, m0=%.0f kg, m_dry=%.0f kg, tilt_max=%d°\n', ...
    fixed.Propulsion.Isp_s, fixed.Propulsion.g0, fixed.Propulsion.Tmax_N, ...
    fixed.Propulsion.m0_kg, fixed.Propulsion.m_dry_kg, fixed.Propulsion.tilt_max_deg);


% (optional) Save best case report
best_report = struct('BestCase',bestCase,'Metrics',metrics,'Fixed',fixed);
json = jsonencode(best_report,'PrettyPrint',true);
fid = fopen('best_case_report.json','w'); fwrite(fid, json); fclose(fid);
save('best_case_report.mat','best_report');  % MATLAB save


%% ===== Histogram of ΔV for successful cases (bins of 50 m/s) =====
okIdx = find(res.success);    % indices of successful runs
dv_ok = res.dv(okIdx);        % ΔV values of successful runs

% bin edges (from min to max, step of 50 m/s)
edges = floor(min(dv_ok)/10)*10 : 10 : ceil(max(dv_ok)/10)*10;

figure('Name','ΔV Histogram — Success only');
histogram(dv_ok, edges, 'FaceColor',[0.2 0.6 0.8], 'EdgeColor','k');
xlabel('ΔV (m/s)');
ylabel('Number of cases');
title(sprintf('ΔV distribution (bin = 50 m/s) for %d successes out of %d runs', numel(dv_ok), nMC));
grid on;


%% ===== Plot of ALL trajectories (green = success, red = fail) =====
okIdx = find(res.success);
koIdx = find(~res.success);

% Performance option: cap the number of plotted trajectories
MAXPLOT = 2000;   
rng(0);
okPlot = okIdx;
koPlot = koIdx;
if numel(okPlot) + numel(koPlot) > MAXPLOT
    % keep ~60% successes if possible
    nOK = min(numel(okIdx), ceil(0.60*MAXPLOT));
    nKO = min(numel(koIdx), MAXPLOT - nOK);
    okPlot = okIdx(randperm(numel(okIdx), nOK));
    koPlot = koIdx(randperm(numel(koIdx), nKO));
    fprintf('[TrajPlot] Sub-sampling for performance: %d successes, %d failures displayed.\n', nOK, nKO);
end

figure('Name','All trajectories — green=success, red=failure'); hold on; grid on;
xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
title('All trajectories (green=success, red=failure)');
set(gca,'XDir','reverse'); view(135,25);

colOK = [0.0 0.6 0.0];    % green
colKO = [0.85 0.15 0.15]; % red

% Plot successful runs
for idx = okPlot(:)'
    si = sims{idx};
    plot3(si.x/1000, si.y/1000, si.z/1000, 'Color', colOK, 'LineWidth', 0.5);
end

% Plot failed runs
for idx = koPlot(:)'
    si = sims{idx};
    plot3(si.x/1000, si.y/1000, si.z/1000, 'Color', colKO, 'LineWidth', 0.5,'');
end

% Mark touchdown points (z=0) for displayed runs
td_ok = []; td_ko = [];
for idx = okPlot(:)'
    si = sims{idx};
    td_ok(end+1,:) = [si.x(end)/1000, si.y(end)/1000, si.z(end)/1000]; 
end
for idx = koPlot(:)'
    si = sims{idx};
    td_ko(end+1,:) = [si.x(end)/1000, si.y(end)/1000, si.z(end)/1000]; 
end
if ~isempty(td_ok),  plot3(td_ok(:,1), td_ok(:,2), td_ok(:,3), 'g.', 'MarkerSize', 8); end
if ~isempty(td_ko),  plot3(td_ko(:,1), td_ko(:,2), td_ko(:,3), 'r.', 'MarkerSize', 8); end

% Add only 2 legend entries
h1 = plot3(nan,nan,nan,'-','Color',colOK,'LineWidth',1.5);
h2 = plot3(nan,nan,nan,'-','Color',colKO,'LineWidth',1.5);
legend([h1 h2], {'Success','Fail'}, 'Location','best');
