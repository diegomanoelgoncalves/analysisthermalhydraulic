function [ Nusselt ] = Nusselt(Re,Pr)
%Calculation Nusselt number
%   Detailed explanation goes here
Nusselt = 0.023*(Re^(0.8))*(Pr^(0.4))
end

