function [ K_fuel ] = K_fuel_UO2( T0 )
%Determination of thermal conductivity of Fuel
%input parameters
%X UO2 - %w.t UO2
%X_ThO2 - %w.t ThO2
A=11.8;
B=0.0238;
C=8.775e-13;
kf= [(1/(A+B*T0))+C*(T0^3) ];
K_fuel=[kf*100];
end

