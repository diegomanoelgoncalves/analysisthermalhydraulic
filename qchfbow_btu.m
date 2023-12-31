function [ qchf_bow ] = qchfbow_btu(de,hfg,h_in,p,g,xaeq)
%The Bowring Critical Heat Flux Mode
%   Detailed explanation goes here
%   Critical heat flux [W/m]
%   hfg - entalpy[J/kg]
%   p - pressure[MPa]
%   g - mass rate[kg/m�s]
%   xaeq - quality[-]
%   de - diameter hydraulic[m]
%   z - axial length[m]
%   Compatibilization units

pr = 0.001*p 
n = 2.0-0.5*pr 
if (pr<=1);
    F1 = (pr^18.942*exp(20.89*(1-pr))+0.917)/1.917 
    F2 = 1.309*F1/(pr^1.316*exp(2.444*(1-pr))+0.309) 
    F3 = (pr^17.023*exp(16.658*(1-pr))+0.667)/1.667 
    F4 = F3*pr^1.649 
else
    F1 = pr^(-0.368)*exp(0.648*(1-pr))
    F2 = F1/(pr^(-0.448)*exp(0.245*(1-pr))) 
    F3 = pr^0.219 
    F4 = F3*pr^1.649 
end 
    B = (de*12)*(g/1e6)/4.0 
    C = 104.4*F3*(de*12)*(g/1e6)/(1+0.347*F4*((g/1e6)^n)) 
    A = 2.317*hfg*B*F1/(1+3.092*F2*(g/1e6)*((de*12)^0.5)) 
    qchf_bow = 3.15*(A-B*hfg*xaeq)/C*1e6;
end

