clc; clear;
% close all;
%% Energy Band Analysis, Electric Field Analysis

% tspan = logspace(-12, 12, 5e3)';
% tspan = [1e-12;1e-11;1e-10;1e-9;1e-8;1e-7;1e-6;1e-5;1e-4;1e-3;1e-2;1e-1;.2;.5;2.^(0:26)'];
% tspan_sym = ...
%     [1E-12; 3E-12;          7E-12;          1E-11;      1.8E-11;    3.4E-11;    6.6E-11;
%     1E-10;  2.07517E-10;    4.22552E-10;    8.52622E-10;    1E-9;     1.86014E-9;    3.58042E-9;
%     7.02098E-9;    1E-8;    1.68811E-8;    3.06433E-8;    5.81678E-8;    1E-7;    2E-7;
%     4.20196E-7;    8.60587E-7;    1E-6;    1.88078E-6;    3.64235E-6;    7.16548E-6;
%     1E-5;    1.70463E-5;    3.11388E-5;    5.93238E-5;    1E-4;    2.1274E-4;    4.38221E-4;
%     8.89181E-4;    1E-3;    0.0019;    0.00371;    0.00731;    0.01;    0.01722;    0.03165;
%     0.06051;    0.1;    0.2;    0.43089;    0.5;    0.96178;    1;    1.92357;    2;    3.84714;
%     4;    7.69427;    8;    15.38854;    16;    30.77708;    32;    61.55417;    64;    123.10834;
%     128;    246.21668;    256;    492.43335;    512;    984.8667;    1024;    1969.7334;    2048;
%     3939.4668;    4096;    7878.9336;    8192;    15757.867;    16384;    31515.735;    32768;
%     63031.469;    65536;    126062.94;    131072;    252125.88;    262144;    504251.75;    524288;
%     1.0085E6;    1.04858E6;    2.01701E6;    2.09715E6;    4.03401E6;    4.1943E6;    8.06803E6;
%     8.38861E6;    1.61361E7;    1.67772E7;    2.67772E7;    3.35544E7;    4.35544E7;    5.35544E7;
%     6.35544E7;    6.71089E7;];
% tspan = tspan_sym;
% tspan_sym = ...
%     [1e-12;1e-11;1e-10;1e-9;1e-8;1e-7;1e-6;1e-5;1e-4;1e-3;1e-2;1e-1;.2;.5;2.^(0:26)'];
% tspan_sym = logspace(-12, 12, 5e3)';

temperature = 298;

Para = InitializeParameters(temperature);


% nt_ctn = 1.1e19; V_ch = -1.2;% 1us program, 2.4-2.5V
% nt_ctn = 0.01e19; V_ch = .2;% 0V
% nt_ctn = 0.2e19; V_ch = .05;% 0.46V
% nt_ctn = 0.405e19; V_ch = -0.1;% 1V
% nt_ctn = 0.643e19; V_ch = -0.25;% 1.5V
% nt_ctn = 0.851e19; V_ch = -.4;% 2V
% nt_ctn = 1.081e19; V_ch = -.55; % 2.5V


% nt_ctn = 0.405e19; V_ch = .3;% 1V, to1=2.0
% nt_ctn = 0.405e19; V_ch = -0.1;% 1V, to1=2.5 % scale variation setting
% nt_ctn = 0.405e19; V_ch = -0.1;% 1V, to1=3.0
% nt_ctn = 0.405e19; V_ch = -0.1;% 1V, to1=3.5
% nt_ctn = 0.405e19; V_ch = -0.1;% 1V, to1=4.0


% nt_ctn = 2e19; V_ch = -0.1;% 1V, to1=4.0

% [Et_trap(:,1), Et_trap(:,2), BE_Conc] = TrapSet_BETox(Para.Et_res, Para.E_bandgap_SiON);

n = 1; %
sigma = .05;
dt = 6e17/sqrt(2*pi)/sigma;
mu = 1.5;

[Et_trap(:,1), Et_trap(:,2), BE_Conc] = TrapSet_BETox(Para.Et_res, Para.E_bandgap_SiON, n, dt, mu, sigma);
nt_tn1 = BE_Conc;


Vg = 0;
Vch = 0;
CTN_slice = 20;
Rad = setRadius(20e-7, 8e-7, 2e-7, 2e-7, 0.9e-7, 6e-7, 8e-7, 2e-7, CTN_slice); % DRL Spec


% Vol = pi*((Rad.r_tn1)^2-Rad.r_tox1^2); % finding Vt in simple uniform distribution
% Vol2 = pi*((Rad.r_ctn)^2-Rad.r_tox^2);
% C1 = 2*pi*Para.e_n/log(Rad.r_tn1/(Rad.r_tox1+0.5*Rad.t_tn1)); %Cylindrical Capacitance: 2pi/e/log(rout/rin)
% C2 = 2*pi*Para.e_ox/log(Rad.r_tox/Rad.r_tn1);
% C3 = 2*pi*Para.e_n/log(Rad.r_ctn/Rad.r_tox);
% C32 = 2*pi*Para.e_n/log(Rad.r_ctn/(Rad.r_tox + 0.5*Rad.t_ctn));
% C4 = 2*pi*Para.e_ox/log(Rad.r_box/Rad.r_ctn);
% Ctg = (1/C1 + 1/C2 + 1/C3 + 1/C4)^(-1);
% CNG = (1/C32 + 1/C4)^(-1);

% [~, ~, ~, ~, ~,~, Vth_ini, Vth_det ] = setFieldParameters(Rad, Para, V_ch, 0, nt_tn1, nt_ctn) % finding Vth by field parameters

% Vth_ini = nt_tn1 *Para.q*Vol/Ctg + nt_ctn * Para.q*Vol2/CNG % Vt calculation with capacitive constants
% Vth_det = nt_tn1 *Para.q*Vol/Ctg
% psi_surface = psi_s(Rad, Para, Vth_ini)

nt_ctn = zeros(CTN_slice,1);
for i = 1:5
    nt_ctn(i) = 4e19;
end

% x = setFieldParameters(Rad, Para, Vch, Vg, nt_tn1, nt_ctn, CTN_slice)

[E_c0, E_v0, E_field0] = GAA_analysis_Band(Rad, Para, Vch, Vg, nt_tn1, nt_ctn, CTN_slice); % 1us program, Vch = -1.7V
E_field1 = abs(E_field0);
Plot_BandDiagram(Rad, E_c0, E_v0, E_field1);

% figure(2);
% plot(Et_trap(:,1), Et_trap(:,2));

% [~,Vt_ratio1] = Vt_Vtratio_SemilogxPlot(Rad, Para, tspan, BETOX_mesh, nt)

% tic
% BETOX_mesh = 15;
% Vt_ratio = Vt_Vtratio_solve(Rad, Para, tspan_sym, BETOX_mesh, nt_tn1, nt_ctn, V_ch, 0.5, Et_trap, BE_Conc);
% % Vt_ratio = Vth_det*Vt_ratio;

%% Stretched Exponential fitting
% toc
% figure(3);
% tspan = tspan;
% semilogx(tspan, Vt_ratio, 'o'); hold on;
% % semilogx(tspan_sym, Vt_ratio, 'o'); hold on;
% % [y_fit, tau_fit, beta_fit] = SE_fit(tspan_sym,Vt_ratio);
% % tau_fit
% % beta_fit
% [tau_compact, beta_compact] = compact_tau(Rad, Para, nt_tn1, 0, nt_ctn, mu, sigma)
% y_compact = nt_tn1*exp(-(tspan./tau_compact).^beta_compact);
% % y_fit = Vt_ratio(1)*exp(-(tspan./.50783).^.2);
% % semilogx(tspan, y_fit); hold on;
% semilogx(tspan, y_compact, LineWidth=2); hold on;

% % plot(fitresult, tspan, Vt_ratio); hold on;

% % Vt_tau_ratio = 1/BETOX_mesh*(1-exp(-1)):1/BETOX_mesh:1;
% % semilogx(tau0, Vt_tau_ratio, 'o'); hold on;

%% nu-TC relation under eTBT
% 
% [nu_Vt, TC_Vt, e_Vt] = nu_TC_relation(Rad, Para, tspan, BETOX_mesh, nt_tn1, nt_ctn, V_ch, 0.5, Et_trap, BE_Conc);
% figure(4)
% subplot(3,1,1);
% plot(nu_Vt(:,1), nu_Vt(:,2)); hold on;
% subplot(3,1,2);
% plot(TC_Vt(:,1), TC_Vt(:,2)); hold on;
% subplot(3,1,3);
% plot(e_Vt(:,1), e_Vt(:,2)); hold on;
