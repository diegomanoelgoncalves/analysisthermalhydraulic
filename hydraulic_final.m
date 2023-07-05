%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program Analysis hydraulic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Parameters
% Consult for more information
%%*https://www.mathworks.com/matlabcentral/
%%fileexchange/9817-x-steam-thermodynamic-properties-of-water-and-steam
%%http://www.iapws.org/relguide/IF97-Rev.html
%% Import txt files power density result from pwd 
file0 = fopen('hydra_boc_nuscale.txt','w');
PDZ1 = importdata('qboc_nuscale.txt');
Z1 =  importdata('z_axis_nuscale.txt');
%Power_density linear [W/cm]
power_density=PDZ1.data;
z = Z1.data;
g=9.81;%[cte-gravity m/s²]
% Count nodes
nodes_z=length(power_density);
nodes_r=length(power_density);
length_active=(z(end)-z(1));
% Axial Length Element
for i=1:(length(z)-1);
    dz(i)=z(i+1)-z(i);
end
dz(end+1)=dz(1);%[cm]
z_plus = cumsum(dz);%[cm]
Z_axis=(abs(z(1))+abs(z(end)));%[cm]
heigth_active=Z_axis/100;%[m]
extrapolated_distance=heigth_active-max(2*z/100);%[m]
%% Calculation Down Pressure, entalpy, quality, MNDBR
% linear density linear power %[W/m]
QL=100*power_density;
%linear power of the pin%[W]
Q1=QL.*(dz/100)';
% linear energy generation rate or linear power %[W/m³]
Q2R=QL./(pi*radius_fuel*radius_fuel);
% linear energy generation rate or linear power %[W/m²]
Q2A=(Q2R.*radius_fuel*radius_fuel)/(2*(radius_fuel));
%heat flux assembly[W/m²]
Q2L=(QL*heigth_active./(pi*radius_fuel^2));
% Heat Flux Volumetric %[W/m³]
Q3=(QL)./(pi*(radius_fuel+thickness_gap+thickness_clad)^2); 

% Axial Length
Z=+extrapolated_distance+(heigth_active);
% Heigth Active add extrapolated distance
H=heigth_active+extrapolated_distance;
HA=linspace(-heigth_active/2,heigth_active/2,nodes_z);
%% Characteristics of design

%Dimensions[m]
radius_clad=radius_fuel+thickness_clad+thickness_gap+thickness_oxidus;
radius_gap=(radius_fuel+thickness_clad);
radius_oxidus=radius_clad+thickness_oxidus;
diameter=2*(radius_fuel+thickness_gap+thickness_clad);%[m]

%Mass rate core[kg/s]
mass_core=(mass_rate)*number_assembly;%(mass_rate)*number_assembly
%Mass rate assembly[kg/s]
mass_a=(mass_rate*(1-core_bypass))/number_assembly;
%Area assembly[m^2]
A_fa=(channel_width^2)-(pi*0.25)*[number_tg*((2*radius_tube)^2)+number_pin*(diameter^2)];
%Flowrate assembly[kg/m^2s]
G_a=abs(mass_a/A_fa);
%Mass rate channel[kg/s]
mass_ch=(mass_a)/(number_assembly);%mass_a/number_rods^2
%Area channel[m^2]
A_ch=[(pitch^2)-(pi*0.25)*(diameter^2)];
%Subchannel flow rate[kg/m^2s]
G_ch=abs(mass_ch/(A_ch));

%average linear heat rate[W/m]
ALHR=(core_power*1e6)/(number_assembly*number_pin*(heigth_active));

%core mass flux
core_mass_flux=(mass_a)/(number_pin*(A_ch));

%% Dimensions calculation
Area_eq=pi*(eq_diameter^2)/4;
%Volume core
Volume_core=Area_eq*heigth_active;
%wetted perimeter
Pw=pi*(2*(radius_clad));
%equivalent_diameter or diameter hydraulic
De=(4*A_ch)/(Pw);

%% Fluid Parameters Extract XSteam parameters 
thermal_pw=100*number_assembly*number_pin*power_density*heigth_active;
% Viscosity dynamic
viscosity_dyn=XSteam('my_pT',Pressure_in,Tm_in-273);   
% Thermal Capacity
Cp=XSteam('Cp_pT',Pressure_in,Tm_in-273);%[kJ/kcal.K]
Cpin=XSteam('Cp_pT',Pressure_in,Tm_in-273);
Cpout=XSteam('Cp_pT',Pressure_in,Tm_out-273);
% Density
density=XSteam('rho_pT',Pressure_in,Tm_in-273);
% Temperature Coolant
T_coolant=XSteam('T_ph',Pressure_in,XSteam('h_pT',Pressure_in,Tm_in-273));
flow_rate=mass_rate./density;%[kg/s/kg/m³=m³/s]
velocity=mass_rate/(density*Area_eq);
%Viscosity dynamic
vd=viscosity_dyn;
viscosity_kyn=viscosity_dyn.*density;
% Fluid Flow constants
Reynolds=(density*velocity*De)/(vd);
Prandlt=(vd*Cp*10^3)/(K_fuel);
% Darcy Friction Factor
friction_factor=64/Reynolds;
% Correlation Dittus-Boelter 
Nusselt_DB = 0.023*(Reynolds^(0.8))*(Prandlt^(0.4));
% Correlation Gnielinski 
Nusselt_G = ((friction_factor/8)*(Reynolds-1000)*Prandlt)/(1+12.7*((friction_factor/8)^(0.5))*((Prandlt^(0.667)-1)));
% Convection fluid
H_fluid=(Nusselt_DB*K_clad)/De;
% Flow Rate
flow_rate=mass_rate./density;%[kg/s/kg/m³=m³/s]

%% Calcuculation Pressure in axial length
Pressure(1)=Pressure_in*1e6/10;
for i=1:length(z);
    % Down Pressure Friction%[Pa]
    dP_fric(i)=-(2*1e-2*z(i)/De)*density*(velocity^2)*0.341*(Reynolds^(-0.25));
    % Down Pressure Form%[Pa]
    dP_form(i)=-(1/2)*(11.63*Reynolds^(-0.25))*density*(velocity^2)*1e-2*z(i);
    % Down Pressure Gravity%[Pa]
    dP_grav(i)=-density*g*1e-2*z(i);
end
%Calculation DowPressure [kPa]
dP=dP_grav+dP_form+dP_fric;
%Calculation Pressure [MPa]
Pressure=-dP+Pressure;
%Calculation Entalpy in liquid [kJ/kg]
hsat_in=XSteam('h_pT',Pressure(1)*1e-5,(Tm_in-273));
%Calculation Entalpy out liquid [kJ/kg]
hsat_out=XSteam('h_pT',Pressure(1)*1e-5,(Tm_out-273));

for i=1:nodes_z;
C_p(i)=XSteam('Cp_pT',Pressure(i)*1e-5,(Tm_in-273));%[kJ/kcal.K]
%Calculation Entalpy liquid [kJ/kg]
hl(i)=XSteam('hL_p',Pressure(i)*1e-5);
%Calculation Entalpy vapour [kJ/kg]
hv(i)=XSteam('hV_p',Pressure(i)*1e-5);
%Calculation Entalpy liquid-vapour [kJ/kg]
hlv(i)=hv(i)-hl(i);
end

% Calcuculation for Entalpy
hf=hsat_in;
for i=1:nodes_z-1;
    hf(i+1)=(hf(i)+(z(i)*1e-2/mass_rate)*((QL(i)-QL(i+1))));
end

%Cálculation local statitc quality
Xe=([hf-hl]./[hv-hl]);
Xin=Xe(1);
Xout=Xe(end);

% Calcuculation Critical Heat Flux uniform and non uniform
for ii=1:length(hf);
    q_w3u(ii)=qchfw3u(De,hv(ii),hl(ii),Pressure(ii),G_a,(Xe(ii)),(z_plus(ii))/100);%G_a(end
    q_w3nu(ii)=qchfw3nu(De,hv(ii),hl(ii),Pressure(ii),G_a,(Xe(ii)),(z_plus(ii)/100),Q1(ii));%G_a(end
    q_bow(ii)=qchfbow(De,1e3*(hlv(ii)),Pressure(ii),G_a,(Xe(ii)));%problema h-entalpia do XSteam
end

% Calcuculation factor peaking
FN=Q2R/mean(Q2R);
% Nuclear Enthalpy Rise Hot Channel Factor 
FNH=((max(QL)*heigth_active)/number_pin)./(((QL)*heigth_active)/number_pin);
pp=1.0;%100% POWER
FNDH=linspace(1.65*(1+0.3*(1-pp)),1.65*(1+0.3*(1-pp)),length(z));

% Calcuculation Critical Heat Flux uniform and non uniform
DNBu=q_w3u'./(Q2A);
DNBnu=q_w3nu'./(Q2A);
DNBbw=q_bow'./(Q2A);

%% Figures analysis hidraulic

figure (1)
hold on
title('Fator de Pico')
plot(z,FN,'.')
xlabel('Axial [cm]')
ylabel('Fn [-]')

figure (3)
hold on
title('Calor específicio a pressão cte')
plot(z,C_p,'.')
xlabel('Axial [cm]')
ylabel('C_{p} [kJ/kg]')

figure (4)
hold on
title('Entalpia')
plot(z,hf,'.')
xlabel('Axial [cm]')
ylabel('h [kJ/kg]')

figure (5)
hold on
title('Pressão')
plot(z,Pressure/1e6,'.')
xlabel('Axial [cm]')
ylabel('Pressão [MPa]')

figure (6)
hold on
title('Título')
plot(z,Xe,'.')
xlabel('Axial [cm]') % left y-axis 
ylabel('Xe [%]') % right y-axis

figure (7)
hold on
title('DNB')
plot(z,(DNBu),'x')
legend('W3-uniform')
xlabel('Axial [cm]') % left y-axis 
ylabel('DNB [-]') % right y-axis
        
figure (8)
hold on
title('Densidade de potência superficial')
plot(z,Q2A*10^-6,'.')
xlabel('Axial [cm]') % left y-axis 
ylabel('Densidade potência superfical [MW/m²]') % right y-axis

figure (9)
hold on
title('W3')
plot(z,q_w3u'*(10^-6),'.',z,q_bow'*(10^-6),'.')
xlabel('Axial [cm]') % left y-axis 
ylabel('Q_{cr} [MW/m²]') % right y-axis
legend('W3-uniform','bowring')

disp('  ******************************  ');
disp('   Results Analysis Hydraulic   ');
disp('  ******************************  ');
disp('   -- Reactor Design--   ');
disp('  Pressure[Pa]    ')
disp(max(Pressure));
disp('  Pressure[Pa]    ')
disp(min(Pressure));
disp('  heigth active[cm]    ')
disp(max(Z_axis));
disp('  channel_width[m]   ')
disp(max(channel_width));
disp('  thickness gap[m]   ')
disp(max(thickness_gap));
disp('  thickness clad[m]   ')
disp(max(thickness_clad));
disp('  Pitch[m]   ')
disp(max(pitch));
disp('  radius Fuel[m]   ')
disp(max(radius_fuel));
disp('  Diameter[m]   ')
disp(max(diameter));
disp('  Diameter hydraulic[m]   ')
disp(max(De));
disp('  Velocity[m/s]   ')
disp(max(velocity));
disp('  Volume core[m³]   ')
disp(max(Volume_core));
disp('  wetted perimeter[m²]   ')
disp(max(Pw));
disp('  -- Heat --  ');
disp('  actual surface heat flux[W/m3]   ')
disp(max(Q2A));
disp('  linear linear power[W/m]   ')
disp(max(Q2L));
disp('  linear power of the pin[W]]   ')
disp(max(Q1));
disp('  Heat volumetric[MW/m3]   ')
disp(max(Q3*1e-6));
disp('  Area channel[m2]   ')
disp(max(A_ch));
disp('  Area subchannel[m2]   ')
disp(max(A_fa));
disp('  Mass subchannel[kg/s]   ')
disp(max(mass_a));
disp('  Mass assembly[kg/s]   ')
disp(max(mass_ch));
disp('  Mass rate[kg/s]   ')
disp(max(mass_rate));
disp('  Flowrate assembly[kgm/s]   ')
disp(max(G_a));
disp('  Flowrate subchannel[kgm/s]   ')
disp(max(G_ch));
disp('  extrapolated_distance[m]   ')
disp(max(Z-heigth_active));
disp(' Flow rate[kg/m²s]   ')
disp(flow_rate);
disp('  actual heat flux[MW/m2]   ')
disp(max(Q2A/1e6));
disp('  Thermal Capacity[kJ/kg]   ')
disp(max(q_w3u/1e6));
disp('  -- Departure Nucleate Boiling --  ');
disp('  DNBR-W3-uniform')
disp(min((DNBnu)));
disp('  DNBR-W3-nouniform  ')
disp(min((DNBu)));
disp('  DNBR-Bowring   ')
disp(min((DNBbw)));
disp('  Thermal Capacity[kJ/kg]   ')
disp(max(C_p));
disp('  Entalpy [kJ/kg]  ')
disp(max(hlv));
disp('  density [kg/m³]  ')
disp(max(density));
disp('  Viscosity dynamic[-]  ')
disp(max(viscosity_dyn));
disp('   -- Constants --   ');
disp('  Prandlt[-]  ')
disp(max(Prandlt));
disp('  Reynolds[-]  ')
disp(max(Reynolds));
disp('  Nusselt Dittus-Boelter Equation%[-]  ')
disp(max(Nusselt_DB));
disp('  Nusselt Gnielinski Equation%[-]  ')
disp(max(Nusselt_G));
disp('  Conv. fluid[W/mK]  ')
disp(max(H_fluid));
disp('  Max total thermal[W/m]   ')
disp(max(thermal_pw));

fprintf(file0,'\n ****************************** \n');
fprintf(file0,' \n Results Analysis Hydraulic \n ');
fprintf(file0,'\n ****************************** \n');
fprintf(file0,'\n\n -- Reactor Design-- \n\n');

fprintf(file0,'\n heigth active[m] %36s %6.6f \n ',max(Z_axis));
fprintf(file0,'\n channel_width[m] %36s %6.6f \n',max(channel_width));
fprintf(file0,'\n thickness gap[m] %36s %6.6f \n',max(thickness_gap));
fprintf(file0,'\n thickness clad[m] %36s %6.6f \n',max(thickness_clad));
fprintf(file0,'\n Pitch[m] %36s %6.6f \n',max(pitch));
fprintf(file0,'\n radius Fuel[m] %36s %6.6f \n',max(radius_fuel));
fprintf(file0,'\n Diameter[m] %36s %6.6f \n',max(diameter));
fprintf(file0,'\n Diameter hydraulic[m] %36s %6.6f \n',max(De));
fprintf(file0,'\n Velocity[m/s] %36s %6.6f \n',max(velocity));
fprintf(file0,'\n Volume core[m³] %36s %6.6f \n',max(Volume_core));
fprintf(file0,'\n wetted perimeter[m²] %36s %6.6f \n',max(Pw));
fprintf(file0,'\n\n -- Heat -- \n\n');
fprintf(file0,'\n actual surface heat flux[W/m3] %36s %6.6f \n',max(Q2A));
fprintf(file0,'\n linear linear power[W/m] %36s %6.6f \n',max(Q2L));
fprintf(file0,'\n linear power of the pin[W]] %36s %6.6f \n',max(Q1));
fprintf(file0,'\n Heat volumetric[MW/m3] %36s %6.6f \n',max(Q3*1e-6));
fprintf(file0,'\n Area channel[m2] %36s %6.6f \n',max(A_ch));
fprintf(file0,'\n Area subchannel[m2] %36s %6.6f \n',max(A_fa));
fprintf(file0,'\n Mass subchannel[kg/s] %36s %6.6f \n',max(mass_a));
fprintf(file0,'\n Mass assembly[kg/s] %36s %6.6f \n',max(mass_ch));
fprintf(file0,'\n Mass rate[kg/s] %36s %6.6f \n',max(mass_rate));
fprintf(file0,'\n Flowrate assembly[kg/s.m^2] %36s %6.6f \n',max(G_a));
fprintf(file0,'\n Flowrate subchannel[kg/s.m^2] %36s %6.6f \n',max(G_ch));
fprintf(file0,'\n extrapolated_distance[m] %36s %6.6f \n',max(Z-heigth_active));
fprintf(file0,'\n Flow rate[kg/m²s] %36s %6.6f \n',flow_rate);
fprintf(file0,'\n moderator entalpy in [kJ/kg] %36s %6.6f \n',max(hlv));
fprintf(file0,'\n actual heat flux[MW/m2] %36s %6.6f \n',max(Q2A/1e6));
fprintf(file0,'\n Thermal Capacity[kJ/kg] %36s %6.6f \n',max(q_w3u/1e6));
fprintf(file0,'\n Quality max [-] %36s %6.6f \n ',max(Xe));
fprintf(file0,'\n Quality min [-] %36s %6.6f \n ',min(Xe));
fprintf(file0,'\n Entalpy max [kJ/kg] %36s %6.6f \n ',max(hf));
fprintf(file0,'\n Entalpy min [kJ/kg] %36s %6.6f \n ',min(hf));
fprintf(file0,'\n Pressure max [Pa] %36s %6.6f \n ',max(Pressure));
fprintf(file0,'\n Pressure min [Pa] %36s %6.6f \n ',min(Pressure));
fprintf(file0,'\n\n -- Departure Nucleate Boiling -- \n\n');
fprintf(file0,'\n DNBR-W3-uniform %36s %6.6f \n',min((DNBnu)));
fprintf(file0,'\n DNBR-W3-nouniform [kJ/kg] %36s %6.6f \n',min((DNBu)));
fprintf(file0,'\n DNBR-Bowring %36s %6.6f \n',min((DNBbw)));
fprintf(file0,'\n Thermal Capacity[kJ/kg] %36s %6.6f \n',max(C_p));
fprintf(file0,'\n Entalpy sat.liquid [kJ/kg]%36s %6.6f \n',max(hl));
fprintf(file0,'\n density [kg/m³]%36s %6.6f \n',max(density));
fprintf(file0,'\n Viscosity dynamic[-]%36s %6.6f \n',max(viscosity_dyn));
fprintf(file0,'\n\n -- Constants -- \n\n');
fprintf(file0,'\n Prandlt[-]%36s %6.6f \n',max(Prandlt));
fprintf(file0,'\n Reynolds[-]%36s %6.6f \n',max(Reynolds));
fprintf(file0,'\n Nusselt Dittus-Boelter Equation%[-]%36s %6.6f \n',max(Nusselt_DB));
fprintf(file0,'\n Nusselt Gnielisnk Equation%[-]%36s %6.6f \n',max(Nusselt_G));
fprintf(file0,'\n Conv. fluid[W/mK]%36s %6.6f \n',max(H_fluid));
fprintf(file0,'\n Max total thermal[W/m] %36s %6.6f \n',max(thermal_pw));
fprintf(file0,'\n -- Results End --');