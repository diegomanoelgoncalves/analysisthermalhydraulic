function [ K_clad,Resist_clad ] = Kzirc(Tc,thickness_clad,radius_fuel)
%Calculation Zirz K_clad [W/m.K]
%   Tc-Temperature clad [ºK]
%   thickness_clad[m]
%   radius_fuel[m]
K_clad=[7.51 + 2.09e-2*Tc-1.45e-5*(Tc^2)+7.67e-9*(Tc^3)];
Resist_clad=[ thickness_clad/(K_clad*radius_fuel)];
end

