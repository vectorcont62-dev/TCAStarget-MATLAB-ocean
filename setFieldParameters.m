function [x, Vth] = setFieldParameters(r, p, V_ch, V_G, nt_tn1, nt_ctn, N)
q = p.q;

% setFieldParameter solves 1. coefficients of Poisson Equation and 2. Vth value

%% Solving Ko1, Kn1, Ko2, KCTN(1-N), KBOX coeff
A = zeros(N+4, N+4); % CTN slice with N + Ko1 Kn1 Ko2 KBOX
b = zeros(N+4,1);

e_o1 = p.e_ox;
e_n1 = p.e_SiON;
e_o2 = p.e_ox;
e_CTN = p.e_n;
e_BOX = p.e_ox;
r_si = r.r_si;
r_o1 = r.r_o1;
r_n1 = r.r_n1;
r_o2 = r.r_o2;
r_CTN = r.r_CTN;
r_BOX = r.r_box;

for i=1:(N+3)
     for j=1:(N+4) % simultaneous equation connection of ... 0 0 1 -1 0 0 ... except last row
          if (i==j)
               A(i,j)=1;
          elseif (i==(j-1))
               A(i,j)=-1;
          else 
               A(i,j)=0;
          end
     end
end
for k=1:(N+4) % simulataneous equation last row
     if (k==1) % Ko1 coeff
          A(N+4,k) = 1/e_o1*log(r_o1/r_si);
     elseif (k==2) % Kn1 coeff
          A(N+4,k) = 1/e_n1*log(r_n1/r_o1);
     elseif (k==3) % Ko2 coeff
          A(N+4,k) = 1/e_o2*log(r_o2/r_n1);
     elseif (k==N+4) % KBOX coeff
          A(N+4,k) = 1/e_BOX*log(r_BOX/r_CTN(N));
     elseif (k==4) % KCTN,1 coeff
          A(N+4,k) = 1/e_CTN*log(r_CTN(1)/r_o2);
     else % KCTN,2-N coeff
          A(N+4,k) = 1/e_CTN*log(r_CTN(k-3)/r_CTN(k-4));
     end
end

for l=1:(N+3) % b solution except last row
     
     if (l==1) % Ko1 coeff
          b(l) = 1/2*q*nt_tn1*r_o1^2;
     elseif (l==2) % Kn1 coeff
          b(l) = -1/2*q*nt_tn1*r_n1^2;
     elseif (l==3) % Ko2 coeff
          b(l) = 1/2*q*nt_ctn(1)*r_o2^2;
     elseif (l==N+3) % KBOX coeff
          b(l) = -1/2*q*nt_ctn(N)*r_CTN(N)^2;
     else % KCTN,2-N coeff
          b(l) = 1/2*q*nt_ctn(l-2)*r_CTN(l-3)^2 - 1/2*q*nt_ctn(l-3)*r_CTN(l-3)^2;
     end
end

b(N+4) = V_G-V_ch...
     - (q*nt_tn1/4/e_n1*(r_n1^2-r_o1^2)+q*nt_ctn(1)/4/e_CTN*(r_CTN(1)^2-r_o2^2));
for ii = 2:N
     b(N+4) = b(N+4)-q*nt_ctn(ii)/4/e_CTN*(r_CTN(ii)^2-r_CTN(ii-1)^2);
end

x = A\b;  % Solution for Ax = b >> x = A^-1 b

%% Vth calculation: Vth = -det[Ab0]

b0 = b;
b0(N+4) = - (q*nt_tn1/4/e_n1*(r_n1^2-r_o1^2)+q*nt_ctn(1)/4/e_CTN*(r_CTN(1)^2-r_o2^2));
for iii = 2:N
     b0(N+4) = b0(N+4)-q*nt_ctn(iii)/4/e_CTN*(r_CTN(iii)^2-r_CTN(iii-1)^2);
end
Ab0 = A;
Ab0(:,1) = b0;
Vth = -det(Ab0);

end