function [qchfw3u] = qchfw3u(De,hf,h,P,G,Xe,z)
%Critical Heat flux - For uniformly heated channels, the critical heat flux is given by
%   Detailed explanation goes here
%hf - entalpy[kJ/kg]
%h - entalpy in [kJ/kg]
%P - pressure[MPa]
%G - mass_rate[kg/m�s]
%Xe - quality[-]
%De - heated diameter equivalent[m]
% P - P/1e6;%[Mpa]
% hf - %[kj/kg]

P=P/1e6;
f1=(2.022-0.06238*P) + (0.1722-0.01427*P)*exp((18.177-0.5987*P)*Xe);
f2=((0.1484-1.596*Xe + 0.1729*Xe*abs(Xe))*(2.236*G) + 3271);
f3=(1.157-0.869*Xe)*(0.2664+0.8357*exp(-124.1*De));
f4=(0.8258+0.0003413*(hf-h));

qchfw3u=f1*f2*f3*f4*1e3; %[W/m�]
end

