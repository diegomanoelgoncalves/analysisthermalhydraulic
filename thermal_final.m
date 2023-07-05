%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program Analysis Thermal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Parameters
%Use XSteam para calcular as propriedades termodinamicas
% Consult for more information
%*** Nomenclature of XSteam code ****
% First the wanted property then a _ then the wanted input properties. 
% Example. T_ph is temperature as a function of pressure and enthalpy. 
% For a list of valid functions se bellow or XSteam for MS Excel.

%% Import txt files power density result from pwd 
file0 = fopen('thermal_boc_sbu.txt','w');
PDZ1 = importdata('qboc_sbu.txt');
Z1 =  importdata('z_axis_nuscale.txt');
power_density=PDZ1.data;
z = Z1.data;
% gravitational acceleration 
g=9.810;%[m/s^2]
% Count nodes
nodes_z=length(power_density);
nodes_r=length(power_density);
% Axial Length Element
for i=1:(length(z)-1);
    dz(i)=z(i+1)-z(i);
end
dz(end+1)=dz(1);%[cm]
z_plus = cumsum(dz);%[cm]
Z_axis=(abs(z(1))+abs(z(end)));%[cm]
heigth_active=Z_axis/100;%[m]
extrapolated_distance=heigth_active-max(2*z/100);%[m]
HA=linspace(-heigth_active/2,heigth_active/2,nodes_z);
%% Characteristics of design
%Dimensions[m]
radius_clad=radius_fuel+thickness_clad+thickness_gap+thickness_oxidus;
radius_gap=(radius_fuel+thickness_clad);
radius_oxidus=radius_clad+thickness_oxidus;
diameter=2*(radius_fuel+thickness_gap+thickness_clad);%[m]

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
%average linear heat rate[kW/m^2]
ALHR=(core_power*1e6)/(number_assembly*number_pin*(heigth_active));
%core mass flux
core_mass_flux=(mass_a)/(number_pin*(A_ch));

%% Dimensions calculation
Area_eq=pi*(eq_diameter^2)/4;
%Volume core
Volume_core=Area_eq*heigth_active;
%wetted perimeter
Pw=pi*(2*(radius_clad));
%equivalent_diameter
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
flow_rate=mass_rate./density;%[kg/s/kg/m³]=[m³/s]
velocity=mass_rate/(density*Area_eq);%[kg/s]/[kg/m³]*[m^2]
%Viscosity dynamic
vd=viscosity_dyn;
viscosity_kyn=viscosity_dyn./density;
% Fluid Flow constants
Reynolds=(density*velocity*De)/(vd);
Prandlt=(1e3*vd*Cp)/(K_clad);
Nusselt = 0.023*(Reynolds^(0.8))*(Prandlt^(0.4));
H_fluid=(Nusselt*K_clad)/De;
%% Method Analytical
%The Equations of Heat Conduction - Nuclear Heat Transfer - El Waki 1971
f=linspace(0,radius_fuel,length(Q3)/10);
g=linspace(radius_fuel,radius_fuel+thickness_gap,length(Q3)/10);
c=linspace(radius_fuel+thickness_gap,radius_fuel+thickness_gap+thickness_clad,length(Q3)/10);
o=linspace(radius_fuel+thickness_gap+thickness_clad,radius_fuel+thickness_gap+thickness_clad+thickness_oxidus,length(Q3)/10);
m=linspace(radius_fuel+thickness_gap+thickness_clad,radius_fuel+thickness_gap+thickness_clad+thickness_oxidus,length(Q3)/10);

Tmod=((max(Q2R)*f(end)^2)/(2*(thickness_gap+thickness_clad+f(end))*H_fluid))+Tm_out;
%To=((max(Q2R)*f(end)^2)/((f(end)+thickness_oxidus)*(K_oxidus./thickness_oxidus)))+Tmod;
Tc=((max(Q2R)*f(end)^2)*(log((f(end)+thickness_clad)/(f(end)+thickness_gap)))/(2*K_clad))+Tmod;
Tg=((max(Q2R)*f(end)^2)/(2*(f(end)+thickness_gap)*H_gap))+Tc;
Tcl=((max(Q2R)*f(end)^2)/(4*K_fuel))+Tg;

for i=1:length(Q3)/10;
T_mod(i)=Tmod-((max(Q2R)*f(i)^2)/(2*(thickness_gap+thickness_clad+f(end))*(H_fluid)));
%T_o(i)=To-((max(Q2R)*f(i)^2)/(2*K_clad))*log((f(i)+thickness_oxidus)/(f(i)+thickness_clad));
T_c(i)=Tc-((max(Q2R)*(f(i)^2)*(log((f(end)+thickness_clad+thickness_gap)/(f(end)+thickness_gap)))/(2*K_clad)));
T_g(i)=Tg-((max(Q2R)*f(i)^2)/(2*(f(i)+thickness_gap)*H_gap));
T_cl(i)=Tcl-((max(Q2R)*f(i)^2)/(4*K_fuel));
end

figure(1)
hold on
title('Temperatura Radial')
plot(1e3*m,T_mod,'o',1e3*c,T_c,'o',1e3*g,T_g,'o',1e3*f,T_cl,'o')
legend('1-moderador','2-revestimento','3-gap','4-combsutível')
xlabel('Radial [cm]')
ylabel('Temperatura [ºK]')

%% Temperature axial Analytical Method El Waki
%Crossectional area
Ac=(pitch^2-pi*f(end)^2);
%Crossectional fuel rod area
Af=(pi*o(end)^2);
%Heigth extrapolated distance
He=(heigth_active+extrapolated_distance);
H=heigth_active;
z0=linspace(z(1)/100,z(end)/100,length(z));
%Others parameters
W=velocity*density*(Ac);%[m/s]*[kg/m^3]*[m^2]
R=radius_fuel+thickness_clad+thickness_gap;
%Circunferecional area
Cc=2*pi*(radius_fuel+thickness_clad+thickness_gap);
%Distance of max temperatures axis
zc=(He/pi)*atan(He/Cp*mass_rate*[0.5*(1/H_fluid*(f(end)+c(end)))]);
zg=(He/pi)*atan(He/Cp*mass_rate*[0.5*(log(f(end)+c(end)/f(end))/K_clad)+1/H_fluid*(f(end)+c(end))]);
zf=(He/pi)*atan(He/Cp*mass_rate*[0.25*(1/K_fuel)+0.5*(log(f(end)+c(end)/f(end))/K_clad)+1/H_fluid*(f(end)+c(end))]);

Rt=radius_fuel+thickness_clad+thickness_gap;
He=(z(end)-z(1))/length(z);
z0=(z)/length(z);
H=heigth_active/length(z);
z_axis=100*(z0+extrapolated_distance);

for i=1:length(z);
    %moderator
     Ta_1(i)=Tmod(end)+((max(Q2R)*He*Ac)/(mass_rate*pi*Cp))*(sin((pi*z0(i))/He)+sin((pi*H)/2*He));
     %cladding
     Ta_2(i)=Tmod(end)+((max(Q2R)*He*Ac)/(mass_rate*pi*Cp))*(sin((pi*z0(i))/He)+sin((pi*H)/2*He))+((max(Q2R)*Ac)/(H_fluid*Cc))*cos(pi*z0(i)/He);
     %gap
     Ta_3(i)=Tmod(end)+((max(Q2R)*He*Ac)/(mass_rate*pi*Cp))*(sin((pi*z0(i))/He)+sin((pi*H)/2*He))+((max(Q2R)*Ac)/(H_fluid*Cc))*cos(pi*z0(i)/He)+((max(Q2R)*Ac)/(2*K_clad))*(log((radius_fuel+thickness_clad+thickness_gap)/radius_fuel+thickness_clad))*cos(pi*z0(i)/He);
     %fuel
     Ta_4(i)=Tmod(end)+((max(Q2R)*He*Ac)/(mass_rate*pi*Cp))*(sin((pi*z0(i))/He)+sin((pi*H)/2*He))+((max(Q2R)*Ac)/(H_fluid*Cc))*cos(pi*z0(i)/He)+((max(Q2R)*Ac)/(H_gap*Cc))*cos(pi*z0(i)/He)+((max(Q2R)*Ac)/(2*K_clad))*(log((radius_fuel+thickness_clad)/radius_fuel))*cos(pi*z0(i)/He)+((max(Q2R)*Ac)/(4*K_fuel))*(log((radius_fuel+thickness_gap)/radius_fuel))*cos(pi*z0(i)/He);
     
end
z_axis=z_axis';

figure(2)
hold on
title('Temperatura Axial')
plot(z_axis,Ta_1,'o',z_axis,Ta_2,'o',z_axis,Ta_3,'o',z_axis,Ta_4,'o')
legend('1-moderador','2-revestimento','3-gap','4-combsutível')
xlabel('Axial[cm]')
ylabel('Temperatura[ºK]')

%% Results
disp('******************************');
disp('  Results Analysis Thermal   ');
disp('******************************');
disp('  -- Reactor Design--  ');
disp('   -- Constants --   ');
disp('  Prandlt[-]')
disp([max(Prandlt)])
disp('  Reynolds[-]')
disp([max(Reynolds)])
disp('  Nusselt%[-]')
disp([max(Nusselt)])
disp('  heigth active[cm]   ')
disp(max(heigth_active))
disp('  channel_width[m]  ')
disp(max(channel_width))
disp('  thickness gap[m]  ')
disp(max(thickness_gap))
disp('  thickness clad[m]  ')
disp(max(thickness_clad))
disp('  Pitch[m]  ')
disp(max(pitch))
disp('  radius Fuel[m]  ')
disp(max(radius_fuel))
disp('  Diameter[m]  ')
disp(max(diameter))
disp('  Diameter hydraulic[m]  ')
disp(max(De));
disp('  Velocity[m/s]  ')
disp(max(velocity))
disp('  Volume core[m³]  ')
disp(max(Volume_core))
disp('  wetted perimeter[m²]  ')
disp(max(Pw))
disp('  Area channel[m2]  ')
disp((A_ch));
disp('  Area subchannel[m2]  ')
disp((A_fa));
disp('  Mass subchannel[kg/s]  ')
disp((mass_a));
disp('  Mass assembly[kg/s]  ')
disp((mass_ch))
disp('  Mass rate[kg/s]  ')
disp((mass_rate));
disp('  Flowrate assembly[kgm/s]  ')
disp((G_a))
disp('  Flowrate subchannel[kgm/s]  ')
disp((G_ch))
disp('   -- Heat  --   ');
disp('  radial linear power[W/m]  ')
disp(max(Q2R))
disp('  actual surface heat flux[W/m3]  ')
disp(max(Q2A))
disp('  linear linear power[W/m]  ')
disp(max(Q2L))
disp('  linear power of the pin[W]]  ')
disp(max(Q1))
disp('  Heat volumetric[MW/m3]  ')
disp(max(Q3*1e-6))
disp('  Flow rate[kg/m²s]  ')
disp(flow_rate)
disp('  density [kg/m³] ')
disp((density))
disp('  Viscosity dynamic[-] ')
disp((viscosity_dyn))
disp('  Conv. fluid[W/m.K] ')
disp((H_fluid));
disp('  Cond. fuel[W/m.K]  ')
disp((K_fuel))
disp('  Cond. gap[W/m.K]  ')
disp((H_gap));
disp('  Cond. clad[W/m.K]')
disp((K_clad))
disp(' -- Temperature -- ');
disp('  Fuel center temperature[ºK]')
disp([max(Tcl)])
disp('  Pelet surface temperature[ºK]')
disp([max(Tg)])
disp('  inner cladding temperature[ºK]')
disp([max(Tc)])
%disp('  outer cladding temperature[ºK]')
%disp([max(To)])
disp('  coolant out[ºK]')
disp([max(T_mod)])
disp(' -- Results End -- ');
fprintf(file0,'\n ******************************\n');
fprintf(file0,'\n Results Analysis Thermal \n ');
fprintf(file0,'\n ******************************\n');
fprintf(file0,'\n -- Reactor Design-- \n');
fprintf(file0,'\n heigth active[cm] %36s %6f \n ',max(Z_axis));
fprintf(file0,'\n channel_width[m] %36s %6.6f \n',max(channel_width));
fprintf(file0,'\n thickness gap[m] %36s %6.6f \n',max(thickness_gap));
fprintf(file0,'\n thickness clad[m] %36s %6.6f \n',max(thickness_clad));
fprintf(file0,'\n mass rate[kg/s] %36s %6.6f \n',(mass_rate));
fprintf(file0,'\n pitch[m] %36s %6.6f \n',max(pitch));
fprintf(file0,'\n radius Fuel[m] %36s %6.6f \n',max(radius_fuel));
fprintf(file0,'\n diameter[m] %36s %6.6f \n',max(diameter));
fprintf(file0,'\n diameter hydraulic[m] %36s %6.6f \n',max(De));
fprintf(file0,'\n velocity[m/s] %36s %6.6f \n',max(velocity));
fprintf(file0,'\n volume core[m³] %36s %6.6f \n',max(Volume_core));
fprintf(file0,'\n wetted perimeter[m²] %36s %6.6f \n',max(Pw));
fprintf(file0,'\n area channel[m2] %36s %6.6f \n',(A_ch));
fprintf(file0,'\n area subchannel[m2] %36s %6.6f \n',(A_fa));
fprintf(file0,'\n mass subchannel[kg/s] %36s %6.6f \n',(mass_a));
fprintf(file0,'\n mass assembly[kg/s] %36s %6.6f \n',(mass_ch));
fprintf(file0,'\n flowrate assembly[kgm/s] %36s %6.6f \n',(G_a));
fprintf(file0,'\n flowrate subchannel[kgm/s] %36s %6.6f \n',(G_ch));
fprintf(file0,'\n extrapolated_distance[m] %36s %6.6f \n',max(z_plus-heigth_active));
fprintf(file0,'\n\n -- Heat  -- \n\n');
fprintf(file0,'\n radial linear power[W/m] %36s %6.6f \n',max(Q2R));
fprintf(file0,'\n actual surface heat flux[W/m3] %36s %6.6f \n',max(Q2A));
fprintf(file0,'\n linear linear power[W/m] %36s %6.6f \n',max(Q2L));
fprintf(file0,'\n linear power of the pin[W]] %36s %6.6f \n',max(Q1));
fprintf(file0,'\n Heat volumetric[MW/m3] %36s %6.6f \n',max(Q3*1e-6));
fprintf(file0,'\n Flow rate[kg/m²s] %36s %6.6f \n',flow_rate);
fprintf(file0,'\n density [kg/m³]%36s %6.6f \n',(density));
fprintf(file0,'\n Viscosity dynamic[-]%36s %6.6f \n',(viscosity_dyn));
fprintf(file0,'\n Conv. fluid[W/m.K]%36s %6.6f \n',(H_fluid));
fprintf(file0,'\n Cond. fuel[W/m.K] %36s %6.6f \n',(K_fuel));
fprintf(file0,'\n Cond. gap[W/m.K] %36s %6.6f \n',(H_gap));
fprintf(file0,'\n Cond. clad[W/m.K] %36s %6.6f \n',(K_clad));
fprintf(file0,'\n\n -- Max Temperature Radial -- \n\n');
fprintf(file0,'\n Fuel center temperature[ºK] %36s %6.6f  \n',[max(T_cl)]);
fprintf(file0,'\n Pelet surface temperature[ºK] %36s %6.6f  \n',[max(T_g)]);
fprintf(file0,'\n inner cladding temperature[ºK] %36s %6.6f  \n ',[max(T_c)]);
fprintf(file0,'\n outer cladding temperature[ºK] %36s %6.6f \n ',[min(T_c)]);
fprintf(file0,'\n coolant out[ºK] %36s %6.6f -%6.6f  \n',[max(T_mod)]);
fprintf(file0,'\n\n -- Max Temperature Axial -- \n\n');
fprintf(file0,'\n Fuel temperature[ºK] %36s %6.6f  \n',[max(Ta_4)]);
fprintf(file0,'\n Gap temperature[ºK] %36s %6.6f  \n',[max(Ta_3)]);
fprintf(file0,'\n Cladding temperature[ºK] %36s %6.6f  \n ',[max(Ta_2)]);
fprintf(file0,'\n Moderator temperature[ºK] %36s %6.6f \n ',[max(Ta_1)]);
fprintf(file0,'\n\n -- Constants -- \n\n');
fprintf(file0,'\n Prandlt[-]%36s%6.6f \n',max(Prandlt));
fprintf(file0,'\n Reynolds[-]%36s%6.6f \n',max(Reynolds))
fprintf(file0,'\n Nusselt%[-]%36s%6.6f \n',max(Nusselt));
fprintf(file0,'\n -- Results End --\n');
