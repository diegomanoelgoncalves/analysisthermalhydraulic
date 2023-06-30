function [qchfw3] = qchfw3(De,h,P,G,Xe,z)
%Critical Heat flux - For uniformly heated channels, the critical heat flux is given by
%   Detailed explanation goes here
%hf - entalpy[kJ/kg]
%h - entalpy in [kJ/kg]
%P - pressure[MPa]
%G - mass_rate[kg/m²s]
%Xe - quality[-]
%De - heated diameter equivalent[m]
hf=XSteam('hL_p',P/1e5)
% P - P/1e6;%[Mpa]
% hf - %[kj/kg]

K1=(2.022-0.06238*P) + (0.1722-0.01427*P)*exp((18.177-0.5987*P)*Xe);
K2=((0.1484-1.596*Xe + 0.1729*Xe*abs(Xe))*(2.236*G) + 3271)*(1.157-0.869*Xe);
K3=(0.2664+0.8357*exp(-124.1*De));
K4=(0.8258+0.0003413*(hf-h));

c=(4.23*10^6)*(((1-Xe)^(7.9))/(G^(1.72)));

qchfw3=K1*K2*K3*K4 %[kW/m²-W/m2]
end

