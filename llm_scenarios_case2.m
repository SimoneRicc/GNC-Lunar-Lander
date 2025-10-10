clc
clear all
close all

%% Constants
g0 = 9.81;              % [m/s^2] Earth gravitational acceleration at sea level

%% Delta-v 
% Lander
dv_sk = 150/365;        % [m/s/day] Delta-v per day for station-keeping in LLO
dv_land = 1.9*1e3;      % [m/s] Delta-v for landing (from LLO to lunar surface)
dv_launch = 1.9*1e3;    % [m/s] Delta-v for lanunching (from lunar surface to LLO)
% Space Tug
dv_llo2geo = 2.05*1e3;  % [m/s] Delta-v for transfer from GEO to LLO and viceversa
dv_rdvz = 100;          % [m/s] Delta-v for randez-vous
% Margin
dv_margin = 1.2;

%% Propulsion system
isp_main = 450;         % [s] Main thruster Specific Impulse
c_main = isp_main*g0;   % [m/s] Main thruster Effective Velocity
mix_ratio = 6;      % [-] Oxygen-to-Hydrogen mass ratio
boil_off_hyd = 0.01;    % [1/day] Hydrogen boil-off per day
boil_off_ox = 0.0025;   % [1/day] Oxygen boil-off per day

%% System Mass
m_dry_llm = 2700;       % [kg] Dry mass of the Launch-Landing Module (LLM)
m_dry_st = 1500;        % [kg] Dry mass of the Space Tug (ST)

%% Model parameters
mission_duration = linspace(0,800,365);  % [days] Duration of the mission from separation to docking
m_target = 3000;        % [kg] Target mass

%% Model
% Space Tug propellant mass
m_prop_st(3) = (m_dry_st+m_target)*(exp(dv_margin*dv_llo2geo/c_main)-1);
m_prop_st(2) = (m_dry_st+m_prop_st(3))*(exp(dv_margin*dv_rdvz/c_main)-1);
m_prop_st(1) = (m_dry_st+m_prop_st(3)+m_prop_st(2))*(exp(dv_margin*dv_llo2geo/c_main)-1);
M_prop_st = sum(m_prop_st);
m_st = m_dry_st+M_prop_st;  % [kg] Wet Mass of the Space Tug


%% Scenario 1: LLM 2 launch-landing
 % 1° landing
    m_fin_1land = m_dry_llm;
    m_in_1land = m_fin_1land*exp(dv_margin*dv_land/(isp_main*g0));
    m_prop_1land = m_in_1land - m_fin_1land;
 % 1° lift-off
    m_fin_1launch = m_in_1land + m_st;
    m_in_1launch = m_fin_1launch*exp(dv_margin*dv_launch/(isp_main*g0));
    m_prop_1launch = m_in_1launch - m_fin_1launch;
 % 2° landing
    m_fin_2land = m_dry_llm + m_dry_st + m_target;
    m_in_2land = m_fin_2land*exp(dv_margin*dv_land/(isp_main*g0));
    m_prop_2land = m_in_2land - m_fin_2land;
 % 2° lift-off
    m_fin_2launch = m_in_2land - m_dry_st - m_target;
    m_in_2launch = m_fin_2launch*exp(dv_margin*dv_launch/(isp_main*g0));
    m_prop_2launch = m_in_2launch - m_fin_2launch;
    
m_prop_phase1 = m_prop_1land + m_prop_1launch;
m_prop_phase2 = m_prop_2land + m_prop_2launch;

fprintf("\n____________\n1° SCENARIO:\n")
fprintf("%35s\t%15s\t%22s\t%10s\n","First Launch","Landing","Second Launch","Landing")
fprintf("%10s\t%18.2f\t%15.2f\t%18.2f\t%14.2f\n","Iniz. Mass kg",m_in_1launch, m_in_1land, m_in_2launch, m_in_2land)
fprintf("%10s\t%18.2f\t%15.2f\t%18.2f\t%14.2f\n","Final Mass kg",m_fin_1launch,m_fin_1land,m_fin_2launch,m_fin_2land)
fprintf("%10s\t%18.2f\t%15.2f\t%18.2f\t%14.2f\n\n","Propellant kg",m_prop_1launch,m_prop_1land,m_prop_2launch,m_prop_2land)

m_prop_tot_scenario1 = m_prop_phase1 + m_prop_phase2;

fprintf("First phase propellant: %10.2f kg\n",m_prop_phase1)
fprintf("Second phase propellant: %9.2f kg\n",m_prop_phase2)
fprintf("Total propellant: %16.2f kg\n",m_prop_tot_scenario1)

%% Scenario 2: LLM Station-keeping
 % Landing
    m_fin_land = m_dry_llm+m_dry_st+m_target;
    m_in_land = m_fin_land*exp(dv_margin*dv_land/c_main);
    m_prop_land = m_in_land-m_fin_land;         % propellant mass w/o considering boil-off
 % Station-Keeping
    m_fin_sk = m_in_land;
    m_in_sk = m_fin_sk*exp(dv_margin*dv_sk*mission_duration/c_main);
    m_prop_sk = m_in_sk-m_fin_sk;
    % With boil-off for landing and station-keeping propellant
    m_hyd_sk = (m_prop_land+m_prop_sk)*1/(mix_ratio+1) .* (1+boil_off_hyd*mission_duration);
    m_ox_sk =  (m_prop_land+m_prop_sk)*mix_ratio/(mix_ratio+1) .* (1+boil_off_ox*mission_duration);
    m_prop_bo = m_ox_sk+m_hyd_sk-(m_prop_land+m_prop_sk);  % propellant mass considering boil-off

 % Launch
    m_fin_launch = m_in_sk+m_st;
    m_prop_launch = m_fin_launch*(exp(dv_margin*dv_launch/c_main)-1);
    % Boil-off
    m_fin_launch_bo = m_in_sk+m_st+m_prop_bo;
    m_in_launch_bo = m_fin_launch_bo*exp(dv_margin*dv_launch/c_main);
    m_prop_launch_bo = m_in_launch_bo-m_fin_launch_bo;

 % LLM total mass
    m_prop_llm = [m_prop_land*ones(size(mission_duration)); m_prop_sk; m_prop_bo; m_prop_launch_bo];
    m_prop_tot_scenario2 = sum(m_prop_llm);
    m_prop_tot_zbo = m_prop_land+m_prop_sk+m_prop_launch;

idx = 6;
fprintf("\n____________\n2° SCENARIO:\n")
fprintf("%30s\t%25s\t%12s\n","Launch","Station-keeping","Landing")
fprintf("%10s\t%18.2f\t%18.2f\t%17.2f\n","Iniz. Mass kg", m_in_launch_bo(idx), m_in_sk(idx), m_in_land)
fprintf("%10s\t%18.2f\t%18.2f\t%17.2f\n","Final Mass kg", m_fin_launch_bo(idx), m_fin_sk, m_fin_land)
fprintf("%10s\t%18.2f\t%18.2f\t%17.2f\n\n","Propellant kg", m_prop_launch_bo(idx), m_prop_sk(idx), m_prop_land)

fprintf("Propellant loss due to boil-off: %16.2f kg\n",m_prop_bo(idx))
fprintf("Total propellant: %16.2f kg\n",m_prop_tot_scenario2(idx))


%% Results
plot(mission_duration, m_prop_tot_scenario1*ones(size(m_prop_tot_scenario2)),'-r', ...
     mission_duration, m_prop_tot_scenario2,'-b', ...
     mission_duration, m_prop_tot_zbo,'--b')
grid on
xlabel('Mission Duration [days]')
ylabel('LLM Propellant Mass [kg]')
legend('Scenario 1', 'Scenario 2 - RBO','Scenario 2 - ZBO')
xlim([mission_duration(1) mission_duration(end)])
ylim([15e3 45e3])