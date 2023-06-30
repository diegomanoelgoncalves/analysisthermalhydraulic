function [ Reynolds ] = Reynolds(d,v,Dhe,vd)
%Reynolds Summary of this function goes here
%Detailed explanation goes here
%density - d [g/cm³]
%*velocity - v[m/s]
%Diameter_hydraulic - Dhe[m]
%(viscosity_dyn) - (vd)(m2/s);
Reynolds=d*v*Dhe/(vd)
end

