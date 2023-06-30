function [ lambda] = Kgap_He(T,P)
%Special reference: Thermophysical Properties of Materials 
%For Nuclear Engineering: A Tutorial and Collection of Data. IAEA-THPH, IAEA, Vienna, 2008. ISBN 978–92–0–106508–7.
% Detailed explanation goes here
% To -temperature in kelvin [K]
% index 1 - for pressure 0.1MPa
% index 2 - for pressure 40MPa
% lambda - conducitvity [W.m/K]

cond(1,:)=[179.0 212.6 244.1 274.0 302.3 328.8 353.3 376.2 397.8 418.8 439.4];%[W.m/K]
cond(2,:)=[180.2 213.5 244.9 274.7 302.8 329.3 353.8 376.6 398.2 419.1 439.7];%[W.m/K]
Te=[100 200 300 400 500 600 700 800 900 1000 1100];
Te=Te+273;%[degree Celsius to Kelvin]
lambda1(1)=interp1(Te,cond(1,:),T);%[W.m/K]
lambda1(2)=interp1(Te,cond(2,:),T);%[W.m/K]

p0=linspace(0.1,40,length(lambda1));
p=[0.1 40];
lambda=interp1(p,lambda1,P);
end

