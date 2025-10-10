%% Values for testing
% Target
a_cmd = [5.3550; 0.2; 1.9];
R = in2body(a_cmd);
q = rotm2quat(R);   
eul = force2eul(a_cmd);

% Inertia
I = 1e5*diag([1.8382, 1.8626, 1.0379])

% Control
kr = [0.5 0.5 0.8]';
ka = 0.2*kr;

% RCS Model
km = 4;
tau = 1;
Uon = 0.4;
Uoff = 0.1;
min_on = 0.03;
min_off = 0.1;
dwell_param = [Uon; Uoff; min_on; min_off];
larm = 2*[4; 4; 4];
Fmax_rcs = 500;
tau_max = larm.*Fmax_rcs
db = 0*[-100 100];


Ts = 0.01; fc = 6;
alpha = exp(-2*pi*fc*Ts);
numd = (1-alpha);
dend = [1 -alpha];

mdl = 'test'; Ts = 0.01;
open_system(mdl);
relays = find_system(mdl,'BlockType','Relay');
for i = 1:numel(relays)
    set_param(relays{i}, ...
        'SampleTime',num2str(Ts), ...
        'ZeroCross','off');
end

function R = in2body(F_cmd)
    
    R1 = @(bank)  [cos(bank) sin(bank) 0; -sin(bank) cos(bank) 0; 0 0 1];
    R2 = @(pitch) [sin(pitch) 0 -cos(pitch); 0 1 0; cos(pitch) 0 sin(pitch)];
    bank = atan2(F_cmd(2,:),F_cmd(1,:));
    F_cmd_rot = R1(bank)*F_cmd;
    pitch = atan2(F_cmd_rot(3,:),F_cmd_rot(1,:));
    R = R2(pitch)*R1(bank);

end


function eul = force2eul(F)
    
    R1 = @(bank)  [cos(bank) sin(bank) 0; -sin(bank) cos(bank) 0; 0 0 1];
    R2 = @(pitch) [sin(pitch) 0 -cos(pitch); 0 1 0; cos(pitch) 0 sin(pitch)];
    bank = atan2(F(2,:),F(1,:));
    F_cmd_rot = R1(bank)*F;
    pitch = atan2(F_cmd_rot(3,:),F_cmd_rot(1,:));
    R = R2(pitch)*R1(bank);
    eul = rotm2eul(R,'ZYX');
    eulx = eul(3); eulz = eul(1);
    eul(1) = eulx; eul(3) = eulz;
    
end

