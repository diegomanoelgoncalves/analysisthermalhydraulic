function [qchfbow] = qchfbow(de,hfg,p,G,xaeq)
%The Bowring Critical Heat Flux Mode
%   Detailed explanation goes here
%   Critical heat flux [W/m]
%   hfg - entalpy[J/kg]
%   p - pressure[MPa]
%   g - mass rate[kg/m²s]
%   xaeq - quality[-]
%   de - diameter hydraulic[m]
%   z - axial length[m]
%   Compatibilization units
    pr = 0.145*(p/1e6);%[MPa]
    na = 2.0-0.5*pr;

    if (pr<=1.00) ;
    F1 =[(pr^18.492)*exp(20.89*(1-pr))+0.917]/1.917;
    F2 = F1/[(pr^1.316)*exp(2.444*(1-pr))+0.309]/1.309;
    F3 = [(pr^17.023)*(exp(16.658*(1-pr)))+0.667]/1.667;
    F4 = F3*(pr^1.649);

    else (pr>1.00);
    F1 = (pr^-0.368)*(exp(0.648*(1-pr)));
    F2 = F1/[(pr^-0.448)*exp(0.648*(1-pr))];
    F3 = (pr^0.219);
    F4 = F3*(pr^1.649);  
    end
    B=(0.25)*(de*G);
    A=(2.317*(hfg*B)*F1)/(1+0.0143*F2*de^(0.5)*G);
    C=(0.077*F3*de*G)/(1+0.347*F4*(G/1356)^na);
    qchfbow = [A-(B*hfg*xaeq)]/C;%[W/m^2]
end

