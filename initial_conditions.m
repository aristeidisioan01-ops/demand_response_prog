%Simulation Parameters for No Demand Response

%% Wind Speed (uwind) - Step Generator Parameters**
% Time (hours)
uwind_time = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 ,17 ,18, 19, 20, 21, 22, 23];

% Amplitude (wind speed m/s)
uwind_amplitude = [8, 3.75, 2.25, 1.75, 1.25, 1, 0.75, 0.5, 1, 1.5, 2, 3, 8, 10, 13, 14.5, 15.5, 16.5, 14, 10.5, 9.5, 8, 6, 5];

%% Irradiation Data**
irradiation_time = [0, 3, 4, 5, 6, 9, 12, 13, 15, 16, 18, 24];
irradiation_amplitude = [0, 5.48, 377, 600, 800, 850, 800, 500, 100, 50, 0, 0];


%% *Temperature Data**
temperature_time = [0, 6, 16, 24];
temperature_amplitude = [18.22, 24.5, 20.8, 20.8];


%% Line Parameters (20kV & 400V)**
% Resistance (Ω/km) * Length (km)
R_20kv = 0.576;  % High Voltage Line Resistance
R_400v = 0.410;  % Low Voltage Line Resistance
km_20kv = 113358/12668;   %Length 20kV line in km
km_400v = 0.2;   %Length 400V line in km

% Inductance (H/km) 
fnom=50;
X_20kv_ohm = 0.397; %L ohms/km
L_H_20kv = X_20kv_ohm/(2*pi*fnom); %L H/km
X_400v_ohm = 0.071; %L ohms/km
L_H_400v = X_400v_ohm/(2*pi*fnom); %L H/km 

assignin('base', 'R_20kv', R_20kv);
assignin('base', 'R_400v', R_400v);
assignin('base', 'L_20kv', L_H_20kv);
assignin('base', 'L_400v', L_H_400v);
assignin('base', 'km_20kv', km_20kv);
assignin('base', 'km_400v', km_400v);

%% Transformer Parameters (Three-Phase, Two-Winding)**
% Rated Power & Frequency
S_rated_trans = 40e3; % 40 kVA
nominal_trans_param = [S_rated_trans, fnom];
assignin('base', 'nominal_trans_param', nominal_trans_param);

% Primary Winding (20kV Side)
winding_1_param = [20e3 , 0.002 , 0.08]; % [V1PH_PH_RMS, R1_pu, L1_pu]
assignin('base', 'winding_1_param', winding_1_param);

% Secondary Winding (400V Side)
winding_2_param = [400 , 0.002 , 0.08]; % [V2PH_PH_RMS, R2_pu, L2_pu]
assignin('base', 'winding_2_param', winding_2_param);

% Magnetization Parameters
Rm_pu = 500; % Magnetization Resistance in pu
Lm_pu = 500; % Magnetization Inductance in pu

% Saturation Curve (Per-unit)
saturation_curve = [0,0; 0.0024,1.2; 1.0,1.52];


assignin('base', 'Rm_pu', Rm_pu);
assignin('base', 'Lm_pu', Lm_pu);
assignin('base', 'saturation_curve', saturation_curve);

%% Supply parame
Ssupply = 7e6;
XRratio = 6;
Vph_ph_rms_nom_supply = 20e3;
assignin('base', 'Ssupply', Ssupply);
assignin('base', 'XRratio', XRratio);
assignin('base', 'Vph_ph_rms_nom_supply', Vph_ph_rms_nom_supply);

%% NO  demand response load W
no_drload_time = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23];
no_drload_amplitude = [5.96, 5.36, 5.01, 5.07, 5.66, 7.03, 8.82, 9.72, 9.83, 9.95, 10.13, 10.25 ,10.25, 10.19, 10.55, 11.74, 13.47, 14.90, 15.73, 15.62, 14.30, 11.98, 9.48, 7.39]*1000;
no_drload_amplitude_24_1 = no_drload_amplitude';

%% Time
stop_time = max(no_drload_time)+1;

%% Final Message
disp('All simulation parameters have been initialized and stored in the workspace.');
