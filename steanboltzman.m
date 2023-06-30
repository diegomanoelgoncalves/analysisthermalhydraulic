function [ sb ] = stefanboltzman(T_fo,T_ci)
sigma=(5.6697*10^(-8));
emiss_fuel=0.15;
emiss_clad=0.20;

sb=sigma*(1/(1/emiss_fuel)+(1/emiss_clad))*((T_fo^4-T_ci^4)/(T_fo-T_ci));
end

