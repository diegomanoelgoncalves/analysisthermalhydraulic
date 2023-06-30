%%%%%%%%%%%%%%%%%%%%%
%INPUT PARAMATERS 
%%%%%%%%%%%%%%%%%
%Reference:NuScale Final Safety Analysis -NuScale Core %Design 

clear all;
clc;

number_pin=264;%[-]
number_assembly=37;%[-]
number_tg=24;%[-]
number_rods=17;%[-]
radius_tube=0.004750;%[-]
thickness_gap=0.00008255;%[m]
thickness_clad=0.0006096;%[m]
thickness_oxidus=25*10^-6;%[m]
radius_fuel=0.00405765;%[m]
pitch=0.012598400;%[m]
mass_rate=587;%[kg/s]
Tm_in=531.26;%[ºK]
Tm_out=587.70;%[ºK]
Temperature_fuel=700;%[K]
eq_diameter=1.51;%[m]
Pressure_in=128.33;%[bar]
percent_down_mass=1.00;%[%]
core_length=17.800;%[m]
core_bypass=8/100;%8.5%[%]
K_fuel=K_fuel_UO2(700);%[W/mK]
K_clad=K_M5(Tm_in);%[W/mK]
H_gap=Hgap_He(531,1e5);%[W/cmK]
K_oxidus=K_oxidus(Tm_in);%[W/mK]
channel_width=0.21;%[m]
core_power=160;%[MW]
