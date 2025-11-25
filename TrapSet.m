function [Et, trap, Conc] = TrapSet(Nt, resEt, Eg)

Et = linspace(0,Eg,resEt)';
trap = zeros(resEt,1);
mu = [0.85, 1.8];
sigma = [0.25, 0.1];
Conc = [1/3, 2/3];
Dt = Nt*Conc./sqrt(2*pi*sigma);


for i = 1:length(mu)
    for ii = 1:resEt
        trap(ii) = trap(ii)+Dt(i)*exp(-((ii*Eg/resEt-mu(i))^2/2/sigma(i)^2))*Eg/resEt;
    end
end