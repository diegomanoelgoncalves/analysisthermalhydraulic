function [ K_fuel ] = K_fuel_UThO2( X_UO2,Temperature_fuel)
%Determination of thermal conductivity of Fuel
%input parameters
%XUO2 - %w.t UO2
%X_ThO2 - %w.t ThO2
% T0 - Tempereture input
%Conductivity K_UO2
A2=11.8;
B2=0.0238;
k_uo2=((1/(A2+B2*T0))+((8.775*10^(-13))*T0^3))*10^2;%[W/mK]

%Conductivity K_ThO2
A1=0.0213;
B1=1.597*10^(-4);

k_tho2=(1/(A1+B1*T0));%[W/mK]
K_fuel=k_uo2*X_UO2+(1-X_UO2)*k_tho2;

end

