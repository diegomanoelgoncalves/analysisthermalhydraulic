function [kf] = condutvUO2B(T0)
kf= [ 1.53E-13*((T0+273)^4)+38.24*log(402.4+T0) ];
end