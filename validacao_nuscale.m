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
Tm_in=531.2610;%[ºK]
Tm_out=598.279;%[ºK]
Temperature_fuel=700;%[K]
heigth_active=2.00;%[m]
eq_diameter=1.51;%[m]
Pressure_in=128.439;%[bar]
percent_down_mass=1.00;
core_length=17.800;%[m]
core_bypass=8.5/100;%[%]
K_fuel=K_fuel_UO2(700);%[W/mK]
K_clad=18;%[W/mK]
H_gap=2/10^-4;%[W/mK]
K_oxidus=0.835+(1.81*10^-4)*Temperature_fuel;%[W/mK]
channel_width=0.214173;%[m]
core_power=160;%[MW]