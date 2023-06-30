function [kf] = condutvUO2C(T0,BU)
kf=[0.053+0.0016*BU+(2.2E-4-5.0E-8*BU)*(T0+273) ];
kf=[1/kf];
end