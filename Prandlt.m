function [ Prandlt ] = Prandlt(vd,C_p,K_fuel)
%Prandlt  Summary of this function goes here
%   Detailed explanation goes here
%vd-viscosity dynamic [(Pa.s)]
%Cp - Specific isobaric heat capacity 
Prandlt=(vd*C_p)/(K_fuel)
end

