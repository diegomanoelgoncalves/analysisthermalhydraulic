function [kf] = condutvUO2A(T0)
A=11.8;
B=0.0238;
C=8.775E-13;
kf= [1/(A+B*T0)+C*(T0^3) ];
kf=[kf*100];
end