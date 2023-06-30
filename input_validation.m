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
radius_fuel=0.004058;%[m]
pitch=0.01259;%[m]
mass_rate=593.6766;%[kg/s]
Tm_in=531;%[�K]
Tm_out=587;%[�K]
Temperature_fuel=700;%[K]
heigth_active=2.00;%[m]
eq_diameter=1.51;%[m]
Pressure_in=128.439;%[bar]
percent_down_mass=1.00;%[%]
core_length=17.800;%[m]
core_bypass=8.5/100;%8.5%[%]
K_fuel=K_fuel_UO2(700);%[W/mK]
K_clad=K_M5(Tm_in);%[W/mK]
H_gap=Hgap_He(700,2e5);%[W/mK]
K_oxidus=K_oxidus(Tm_in);%[W/mK]
channel_width=0.19;%[m]
core_power=160;%[MW]