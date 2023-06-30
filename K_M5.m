function [lambda1] = K_M5(tc)
%Calculation Zirz K_clad [W/m.K]
%   Tc-Temperature clad [ºK]
%   thickness_clad[m]
%   radius_fuel[m]
% Reference : Development of M5 Cladding
% Material Correlations in the TRANSURANUS Code
% lambda1 - FRAPCON (Luscher et al. 2011)
% lambda2 - RELAP5-NRC-approved AREVA report (Mitchel et al. 2000)
    if tc>2133;
    lambda1=36
    else tc<=2133;
    %lambda1=15.0636*exp(0.000461843*tc);
    lambda1 = 8.6383*exp(0.0007*tc)
    end
end