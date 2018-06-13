clc;
clear all;
close all;
%% parameter 
cT0=3;%relative power extraction coefficient, kJ/kg
beta=0.01;%extraction of air from high pressure compressor to cabin and so on
delta1=0.08;%extraction of air from high pressure compressor to high pressure turbine
delta2=0.02;%extraction of air from high pressure compressor to low pressure turbine
% xigemai=0.97; %Total pressure recovery coefficient of inlet
xigemab=0.97; %Total pressure recovery coefficient of main combustor
xigemam=0.97; %Total pressure recovery coefficient of mixer
xigemabab=0.97; %Total pressure recovery coefficient of afterburner
xigemac=0.97; %Total pressure recovery coefficient of nozzle
xigemaII=0.98; %Total pressure recovery coefficient of the outer flow duct from fan exit to mixer entrance
ngLPC=0.870;% adiabatic efficiency of fan
ngHPC=0.879;% adiabatic efficiency of high pressure compressor
ngb=0.99;% adiabatic efficiency of combuster
ngHPT=0.90;% adiabatic efficiency of high pressure turbine
ngLPT=0.917;% adiabatic efficiency of low pressure turbine
ngbab=0.97;% adiabatic efficiency of aftbunner
ngmH=0.99;%mechanical efficiency of high pressure axis
ngmL=0.99;%mechanical efficiency of low pressure axis
k=1.4;%  adiabatic exponent of air
cp=1.005;%cp pf air kJ/(kg*K)
R=(k-1)/k*cp*1000;%R of air,J/(kg*K)
kg=1.3;%  adiabatic exponent of gasÈ¼Æø
cpg=1.244;%cp pf gasÈ¼Æø kJ/(kg*K)
Rg=(kg-1)/kg*cpg*1000;%R of agas,J/(kg*K)
Hu=42900;% low heating value of fuel,kJ/kg
%% design point
Ma0=2.0;%flying Mach number
H=11;%flying height,km
%% thermal cycle parameter
B=0.38; %by pass ratio
piLPC=3.95;%  pressure ratio of fan 
piHPC=5.216;%  pressure ratio of high pressure compressor
pic=piLPC*piHPC;% total pressure ratio
Tt4=1900;%temperature before turbine, K
Ttab=2050;%temperature of the exit of aftbunner
%% calculation 0 section
[p0,T0,rho0]=US_Standard_Air1976(H*1000);%atmosphere static temperature,static pressure, density
a0=sqrt(k*R*T0);%speed of sound
c0=a0*Ma0;%speed
pt0=p0*(1+(k-1)*Ma0^2/2)^(k/(k-1));%total pressure 
Tt0=T0*(1+(k-1)*Ma0^2/2);%total temperature
%% calculation 2 section
if Ma0<1
    xigemai=0.97;%Total pressure recovery coefficient of inlet
else
    xigemai=0.97*(1-0.075*(Ma0-1)^1.35);
end
pt2=xigemai*pt0;
Tt2=Tt0;
%% calculation 22 section
pt22=pt2*piLPC;
Tt22=Tt2*(1+(piLPC^((k-1)/k)-1)/ngLPC);
LLPC=cp*(Tt22-Tt2);%power consumed by fan per weight air
%% calculation 3 section
pt3=pt22*piHPC;
Tt3=Tt22*(1+(piHPC^((k-1)/k)-1)/ngHPC);
LHPC=cp*(Tt3-Tt22);
%% calculation 4 section
pt4=pt3*xigemab;
f=(cpg*Tt4-cp*Tt3)/(ngb*Hu-cpg*Tt4);%the mass ratio between fuel and air,=qmf/((qma3*(1-beta-delta1-delta2)))
%% calculation 4.5 section high pressure turbine
taum1=((1-beta-delta1-delta2)*(1+f)+cp*delta1*Tt3/(cpg*Tt4))/((1-beta-delta1-delta2)*(1+f)+delta1);%Tt4a/Tt4
Tt4a=Tt4*taum1;%4a:afterthe mixing with the cooling air delta1
pt4a=pt4;
Tt45=Tt4a*(1-LHPC/(((1-beta-delta1-delta2)*(1+f)+delta1)*ngmH*cpg*Tt4a));
piHPT=(1-(1-Tt45/Tt4a)/ngHPT)^(-kg/(kg-1));
pt45=pt4a/piHPT;
%% calculation 5 section low pressure turbine
taum2=((1-beta-delta1-delta2)*(1+f)+delta1+cp*delta2*Tt3/(cpg*Tt45))/((1-beta-delta1-delta2)*(1+f)+delta1+delta2);%Tt4c/Tt45
Tt4c=taum2*Tt45;
pt4c=pt45;
Tt5=Tt4c*(1-(LLPC+cT0/ngmL)*(1+B)/(((1-beta-delta1-delta2)*(1+f)+delta1+delta2)*ngmL*cpg*Tt4c));
piLPT=(1-(1-Tt5/Tt4c)/ngLPT)^(-kg/(kg-1));
pt5=pt4c/piLPT;
%% calculation 6 section mixer
Bm=B/((1-beta-delta1-delta2)*(1+f)+delta1+delta2);%the real bypass ratio of the entrance of mixer
cp6=(cpg+Bm*cp)/(1+Bm);
Tt6=Tt5*cpg/cp6*(1+Bm*cp*Tt22/(cpg*Tt5))/(1+Bm);
pt5II=xigemaII*pt22;
pm=(pt5+Bm*pt5II)/(1+Bm);
pt6=xigemam*pm;
%% calculation 7 section aftbunner
pt7=pt6*xigemabab;
Tt7=Ttab;
fab=(1+f*(1-beta-delta1-delta2)/(1+B-beta))*(cpg*Tt7-cpg*Tt6)/(ngbab*Hu-cpg*Tt7);
f0=(f*(1-beta-delta1-delta2)+(1+B-beta)*fab)/(1+B);
%% calculation 9 section nozzle
pt9=xigemac*pt7;
Tt9=Tt7;
p9=p0;
Ma9=sqrt(2/(kg-1)*((pt9/p9)^((kg-1)/kg)-1));
T9=Tt9*(1+(kg-1)/2*Ma9^2)^(-1);
a9=sqrt(kg*Rg*T9);
c9=a9*Ma9;
%% performance
Fs=(1+f0-beta/(1+B))*(c9+R*T9/c9*(1-p0/p9))-c0
sfc=3600*f0/Fs