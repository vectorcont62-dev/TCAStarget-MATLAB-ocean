function [Ec, Ev, Efield] = GAA_analysis_Band(r, p, Vch, Vg, nt_n1, nt_ctn, N)

x = setFieldParameters(r, p, Vch, Vg, nt_n1, nt_ctn, N);
Ko1 = x(1);
Kn1 = x(2);
Ko2 = x(3);
KCTN = x(4:N+3);
KBOX = x(N+4);
V_c = zeros(length(r.r),1);
Ec = zeros(length(r.r),1);
Ev = zeros(length(r.r),1);
Efield = zeros(length(r.r),1);
q = p.q;

V_instant_CTN = zeros(N,1);

for i=1:length(r.r)
    if (r.r(i)<r.r_0)
        V_instant_ch = 0;
        V_c(i) = V_instant_ch;
        Ec(i) = V_c(i);
        Ev(i) = Ec(i) - p.E_bandgap_Si;
        Efield(i) = 0;

    elseif(r.r(i)>=r.r_0 && r.r(i)<r.r_tox1)
        V_instant_o1 = Ko1/p.e_ox*log(r.r(i)/r.r_0);
        V_c(i) = V_instant_o1; % instant potential utilized: not yet used 
        Ec(i) = p.PHI_B_SiO2-V_c(i);
        Ev(i) = Ec(i) - p.E_bandgap_SiO2;
        Efield(i) = - Ko1/p.e_ox/r.r(i);    

    elseif(r.r(i)>=r.r_tox1 && r.r(i)<r.r_tn1)
        V_instant_n1 = Ko1/p.e_ox*log(r.r_tox1/r.r_0) ...
            + Kn1/p.e_SiON*log(r.r(i)/r.r_tox1) + q*nt_n1/4/p.e_SiON*(r.r(i)^2-r.r_tox1^2);
        V_c(i) = V_instant_n1; % instant potential utilized: not yet used 
        Ec(i) = p.PHI_B_SiON-V_c(i);
        Ev(i) = Ec(i) - p.E_bandgap_SiON;
        Efield(i) = - Kn1/p.e_SiON/r.r(i) - q*nt_n1/2/p.e_SiON*r.r(i);

    elseif(r.r(i)>=r.r_tn1 && r.r(i)<r.r_tox)
        V_instant_o2 = Ko1/p.e_ox*log(r.r_tox1/r.r_0) ...
            + Kn1/p.e_SiON*log(r.r_tn1/r.r_tox1) + q*nt_n1/4/p.e_SiON*(r.r_tn1^2-r.r_tox1^2)...
            + Ko2/p.e_ox*log(r.r(i)/r.r_tn1);
        V_c(i) = V_instant_o2; % instant potential utilized: not yet used 
        Ec(i) = p.PHI_B_SiO2-V_c(i);
        Ev(i) = Ec(i) - p.E_bandgap_SiO2;
        Efield(i) = - Ko2/p.e_ox/r.r(i);

    elseif (r.r(i)>=r.r_tox && r.r(i)<r.r_CTN(1))
        V_instant_CTN(1) = V_instant_o2 + KCTN(1)/p.e_n*log(r.r(i)/r.r_o2) + q*nt_ctn(1)/4/p.e_n*(r.r(i)^2-r.r_o2^2);
        V_c(i) = V_instant_CTN(1);
        Ec(i) = p.PHI_B_SiN-V_c(i);
        Ev(i) = Ec(i) - p.E_bandgap_SiN;
        Efield(i) = - KCTN(1)/p.e_n/r.r(i) - p.q*nt_ctn(1)/2/p.e_n*r.r(i);

    elseif (r.r(i)>=r.r_CTN(1) && r.r(i)<r.r_CTN(N))
        slice_num = SliceDetection(r, r.r(i), N);
        V_instant_CTN(slice_num) = V_instant_CTN(slice_num-1) + KCTN(slice_num)/p.e_n*log(r.r(i)/r.r_CTN(slice_num-1)) + q*nt_ctn(slice_num)/4/p.e_n*(r.r(i)^2-r.r_CTN(slice_num-1)^2);
        V_c(i) = V_instant_CTN(slice_num);
        Ec(i) = p.PHI_B_SiN-V_c(i);
        Ev(i) = Ec(i) - p.E_bandgap_SiN;
        Efield(i) = - KCTN(slice_num)/p.e_n/r.r(i) - p.q*nt_ctn(slice_num)/2/p.e_n*r.r(i);

    elseif (r.r(i)>= r.r_CTN(N) && r.r(i)<r.r_box)
        V_instant_BOX = V_instant_CTN(N) + KBOX/p.e_ox*log(r.r(i)/r.r_CTN(N));
        V_c(i) = V_instant_BOX;
        Ec(i) = p.PHI_B_SiO2 - V_c(i);
        Ev(i) = Ec(i) - p.E_bandgap_SiO2;
        Efield(i) = - KBOX/p.e_ox/r.r(i);

    else
        V_c(i) = Vg;
        Ec(i) = -Vg;
        Ev(i) = Ec(i) - 0;
        Efield(i) = 0;
    end

end

% V_diff = diff(V_c);
% Efield = -[V_diff(1,1),V_diff']'/r.dr;
end