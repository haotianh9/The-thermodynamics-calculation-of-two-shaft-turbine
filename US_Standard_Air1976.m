function [p,T,rho]=US_Standard_Air1976(alt_geom)
% Compute the pressure, density, and temperature of standard air
% at a given altitude.
% Standard according to US Standard Air 1976
%
% alt_geom: geometric altitude above sea level, meter.
%           i.e., height above sea level;
% p   : air pressure, Pa
% rho : air density, kg/m^3
% T   : air temperature, degree K.

% ====== Part (1): setting up ======
% Convert geometric altitude to geopotential altitude
r_Earth=6.369e+3; % radius of earth, km
alt_geom=alt_geom*0.001; % from meter to km
alt_geop=r_Earth*alt_geom/(r_Earth+alt_geom); % km

% The atmosphere is divided into 8 layers
% The lowest geopotential altitude of each layer is called the base
% altitude
N_Layer=8;
alt_geop_base=[0.0 11.0 20.0 32.0 47.0 51.0 71.0 84.852]; % km

% Determine which layer it is in
layer=0;
for i1=1:N_Layer;
    if alt_geop>=alt_geop_base(i1) && alt_geop<alt_geop_base(i1+1);
        layer=i1;
        break;
    end;
end;
layer;
delta_geop=alt_geop-alt_geop_base(layer);


% ====== Part (2): calculating p/p_sl, T, rho/rho_sl ======
% In each layer, T is a linear function of the geopotential altitude
% T(alt_geop)=T_base(layer)+T_grad(layer)*(alt_geop - alt_geop_base)
T_base=[288.15 216.65 216.65 228.65 270.65 270.65 214.65 186.946];
T_grad=[-6.5 0.0 1.0 2.8 0.0 -2.8 -2.0 0.0];
T=T_base(layer)+T_grad(layer)*delta_geop; % K

% Pressure
p_base=[1.0 2.233611e-1 5.403295e-2 8.5666784e-3 1.0945601e-3 6.6063531e-4 3.9046834e-5 3.68501e-6];
%p_base=p_base*1.01325e+5; % Pa
GMR=3.4163195e+1;
if layer==2 || layer==5 || layer==8;
    p=p_base(layer)*exp(-GMR*delta_geop/T_base(layer));
else;
    p=p_base(layer)*(T_base(layer)/T)^(GMR/T_grad(layer));
end;

% Density
rho=p*T_base(1)/T;

% ====== Part (2): calculating real p, rho ======
% Back to real values
rho0=1.225; % rho_sl, kg/m^3
p0=1.01325e+5; % p_sl, Pa
p=p*p0;
rho=rho*rho0;
