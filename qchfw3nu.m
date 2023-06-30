function [qchfw3nu] = qchfw3nu(De,hf,h,Pr,G,Xe,z,qdot)
%Critical Heat flux - For uniformly heated channels, the critical heat flux is given by
%   Detailed explanation goes here
%hf - entalpy[kJ/kg]
%h - entalpy in [kJ/kg]
%P - pressure[MPa]
%G - mass_rate[kg/m²s]
%Xe - quality[-]
%De - heated diameter equivalent[m]
%qdot - Heat Flux[W/m]
% P - P/1e6;%[Mpa]
% hf - %[kj/kg]
P=Pr/1e6;
f1=(2.022-0.06238*P) + (0.1722-0.01427*P)*exp((18.177-0.5987*P)*Xe);
f2=((0.1484-1.596*Xe + 0.1729*Xe*abs(Xe))*(2.236*G) + 3271);
f3=(1.157-0.869*Xe)*(0.2664+0.8357*exp(-124.1*De));
f4=(0.8258+0.0003413*(hf-h));
c=(4.23*1e6)*(((1-Xe)^(7.9))/(G^(1.72)));
qchfw3u=f1*f2*f3*f4 ;%[kW/m²]

fun = @(z) qdot*exp(-c*(z(end)-z));
q1 = integral(fun,0,z(end));
q2 =  qdot*(1-exp(-c*(z(end))));
f=(c*q1)/q2;

qchfw3nu=(qchfw3u./f)*1e3; %[W/m²]
end

