function TC=TC_value(x1, x2, Phi_1, Phi_2, Field, effectivemass, p)

m = effectivemass*p.m_e; % [kg]
hbar = p.h_bar; %[Js]
q = p.q; %[C]
Feq = 100*Field; %Field: [V/cm] -> Feq: [V/m]


H = max(Phi_1, Phi_2)*p.q; %Phi_h: [eV] -> H: [J]
H = max(0,H);
L = min(Phi_1, Phi_2)*p.q; %Phi_l: [eV] -> H: [J]
L = max(0,L);

if(Feq<1)
    TC_exp = -2*sqrt(2*m)/(hbar*q)*sqrt(H^1.5)*(abs(x1-x2)/100);
    TC = exp(TC_exp);
else
    TC_exp = -4*sqrt(2*m)/(3*hbar*q)*(H^1.5-L^1.5)/Feq;
    TC = exp(TC_exp);
end