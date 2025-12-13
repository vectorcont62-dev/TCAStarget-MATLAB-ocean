clc; clear;
% close all;
%% Parameter initialize

% tspan = logspace(-12, 12, 5e3)';

temperature = 298;

Para = InitializeParameters(temperature);

Vg = 0;
Vch = 0;
CTN_slice = 20;
Et_res = 100;

[Et, trap, Conc] = TrapSet(3e19, Et_res, Para.E_bandgap_SiN);

Rad = setRadius(20e-7, 8e-7, 2e-7, 2e-7, 0.9e-7, 6e-7, 8e-7, 2e-7, CTN_slice); % DRL Spec BETOX
% Rad = setRadius(20e-7, 8e-7, 5e-7, .01e-7, .01e-7, 6e-7, 8e-7, 2e-7, CTN_slice); % DRL Spec ONO

cell_info_initial = cell(Et_res,3+CTN_slice);
for i = 1:Et_res
    for j = 1:3+CTN_slice
        cell_info_initial{i,j} = zeros(5,1); %radius, trap depth, current conc, nu, TC
        if(j==1)
            cell_info_initial{i,j}(1) = Rad.r_si+0.5*Rad.t_o1;
        elseif(j==2)
            cell_info_initial{i,j}(1) = Rad.r_o1+0.5*Rad.t_n1;
        elseif(j==3)
            cell_info_initial{i,j}(1) = Rad.r_n1+0.5*Rad.t_o2;
        else
            cell_info_initial{i,j}(1) = Rad.r_o2+(j-3.5)*Rad.t_ctn/CTN_slice;
        end
        cell_info_initial{i,j}(2) = (i-0.5)*Para.E_bandgap_SiN/Et_res;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CTN Programmed Poisson Solution & Vth
nt_tn1 = 0;
nt_ctn = zeros(CTN_slice,1);
for i = 1:5
    nt_ctn(i) = 3e19;
end
nt_ctn(1) = 1e19;

Vg = -.25;

for j = 1:3+CTN_slice
    if(j>=4)
        [Et_temp, trap_temp, ~] = TrapSet(nt_ctn(j-3), Et_res, Para.E_bandgap_SiN);
        for i = 1:Et_res
            cell_info_initial{i,j}(3) = trap_temp(i);
        end
    end
end

[~,Vth] = setFieldParameters(Rad, Para, Vch, Vg, nt_tn1, nt_ctn, CTN_slice)

[E_c0, E_v0, E_field0] = GAA_analysis_Band(Rad, Para, Vch, Vg, nt_tn1, nt_ctn, CTN_slice); % 1us program, Vch = -1.7V
E_field1 = abs(E_field0);
Plot_BandDiagram(Rad, E_c0, E_v0, E_field1);

N = zeros(Et_res,1);
for i = 1:Et_res
    N(i) = cell_info_initial{i,8}(3);
end
figure(2);
plot(Et, N)
%% nu-TC relation under eTBT
for j = 1:3+CTN_slice
    if(j>=4)
        for i = 1:Et_res
            
        end
    end
end


[TC_ch, TC_WL] = calc_TC_quantized(Rad, Para, E_c0, E_v0, E_field1, cell_info_initial{40,2})
[nu_ch, nu_WL] = nu_value(Rad, Para, E_c0, E_v0, E_field1, cell_info_initial{40,2})