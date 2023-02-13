clc, clear, close all;
%this is the V-n diagram projecr for the group of
%Benjamin Young, Jordan Decker, and Michael Blunt
%problem 1
%1 in the variable name denotes for number 1
%Given constants/related values
%m1 = 2300; %  mass in kg
W1 = 15;%m1 * 9.813; %Weight in lbsp
S1 = 8;%19.33; %planform area in m^2 
AR1 = 8; %Aspect ratio
Cla1 = 4.48; %lift curve slope in /rad
Clmaxes1 = [1.398 .6]; %cl max and clmax neg

nmax = 4.4; %n+ from FAR23 for aerobatic
nmin  = -0.4*nmax; %nmin from FAR23 for aerobatic
rho = 0.002377;% density in kg/m^2 at 10kft
Vc1 = 45*1.688; %Design Cruise Speed (EAS) in kts at 10000
b1 = sqrt(AR1*S1); %span in m
c1 = S1/b1; %chord 
Vd1 = 1.55*Vc1; %dive Speed
n2 = 0.75*nmax;%n2 for problem 1
Udec1 = 50; 
Uded1 = Udec1/2;
[npar, Vpar, Vspos, Vsneg,Va,Vg] = VNparabola(W1,S1,Clmaxes1, nmax, nmin,rho);
[Vgust,ngust] = VNgust(W1,S1,rho,c1,Cla1,Vd1,Udec1,Uded1);

VNplot1=VNplot(Vpar, npar, Va, Vc1, Vd1, Vg, Vspos, Vsneg, nmax, nmin, Vgust, ngust,n2);

%problem 2
m2 = 15/32.2; %  mass in kg
W2 = m2 * 32.2;
S2 = 8; %planform area in m^2 
Clmaxes2 = [1.4,0.8]*1.25; %Clmax+ and -
Cd0 = 0.0101;
AR2 = 8; %Aspect ratio
b2 = sqrt(AR2*S2); %span in m
c2 = S2/b2; %chord 
Vc2 = 45*1.688; %Design Cruise Speedin m/s
Vmax = 50*1.688; %maximum speed
Vd2 = Vmax;
oswald = Clmaxes2(1)^2/(Cd0*pi*AR2); %calculating an e0 for the Cla
%a0 = 4.48; % assuming ideal
Udec2 = 50;
Uded2 = Udec2/2;
Cla2 = 4.48;%a0/(sqrt(1+(a0/(pi*oswald*AR2)))+((2*pi)/(pi*oswald*AR2))); %Clalpha assuming ideal slope
[npar2, Vpar2, Vspos2, Vsneg2,Va2,Vg2] = VNparabola(W2,S2,Clmaxes2, nmax, nmin,rho);
[Vgust2,ngust2,mug2] = VNgust(W2,S2,rho,c2,Cla2,Vd2,Udec2,Uded2);

VNplot2=VNplot(Vpar2, npar2, Va2, Vc2, Vd2, Vg2, Vspos2, Vsneg2, nmax, nmin, Vgust2, ngust2,n2);

