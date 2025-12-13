function [nu_ch, nu_WL] = nu_value(r, p, E_c0, E_v0, E_field1, cell_info)

[idx, radius_quantized, error] = getRadiusIndex(cell_info(1), r);

TrapVolume = 4.64492e-12*1e-18;
W = 2*pi/p.h_bar*cell_info(2)^2*TrapVolume;
Etrapheight_ch = E_c0(idx)-cell_info(2);
Etrapheight_ch = max(0,p.q*Etrapheight_ch);
Etrapheight_WL = E_c0(idx)-E_c0(r.n_num)-cell_info(2);
Etrapheight_WL = max(0,p.q*Etrapheight_WL);

D_ch = sqrt(2*p.m_e^3*.40)/pi^2/p.h_bar^3*sqrt(Etrapheight_ch);
D_WL = sqrt(2*p.m_e^3*.40)/pi^2/p.h_bar^3*sqrt(Etrapheight_WL);

nu_ch=D_ch*W;
nu_WL=D_WL*W;

end