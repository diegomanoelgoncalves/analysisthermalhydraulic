function [ H_gap ] = Hgap_He(T0,P0)
%Special reference: Gas gap conductance in contact heat transfer (Wahid.S 2001)

a=0.425-2.3*(10^(-4))*T0
R=8.314 %[J/mole K]
CP=20.967;
CV=12.863;
gama=CP/CV;
H_gap=(P0/4)*((2*R/pi)^(0.5))*((gama+1)/(gama-1))*(1/T0^(0.5))*(a/(2-a))
end

