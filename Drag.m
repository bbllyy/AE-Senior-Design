%Drag Buildup
clc, clear, close all
%Aircraft Characteristics
W = 13;
b = 8;
cw = 1;
Sw = b*cw;
AR = b^2/Sw;%6;
e = 0.85;
tc_max = .12;
bh = 2+(2/3);
ch = (2/3);
Sh = bh*ch;
ARh = bh^2/Sh;

bv = 1+(1/3);
cvr = 12/12;
cvt = 8.5/12;
cv = (cvr+cvt)/2;
Sv = ((cvr+cvt)/2)*bv;
ARv = bv^2/Sv;

%Control Surface sizes
ce = 0.3*ch;
cr = 0.3*cv;
ca = 0.3*cw;

%Wing + Htail characteristics
CLmax = 1.41;
CLa = 4.894; %lift slope per rad;
CLadeg = 0.85416; %Lift slope per degree
alpha0 = -0.0622;
alpha0_deg = -3.564;
CL0 = 0.3114;
CD0 = 0.01176;

lfuse = 22.9/12;
dfuse = 4.5/12;
Swet_fuse = 2*(lfuse*dfuse+dfuse*dfuse+lfuse*dfuse);
Sref = Sw+Sv+lfuse*dfuse;
%Flight characteristics
rho = 0.002377;%0.002377;
mu = 3.737E-7;
nu = mu/rho;
V = 40*1.688;
aoa_lo = 3;
aoa_td = 3;
aoa_cr = 0.5;
Vstall = sqrt((2*W)/(rho*Sw*CLmax));

[CL_lo,CD_lo] = DragBuildup(aoa_lo, 'PlaneInfo.mat');
[CL_td,CD_td] = DragBuildup(aoa_td, 'PlaneInfo.mat');
[CL_cr,CD_cr] = DragBuildup(aoa_cr, 'PlaneInfo.mat');
mur = 0.05;
T = 12; Tr = 7;
Vlo = 1.2*Vstall;
Dlo = CD_lo*0.5*rho*0.7*Vlo^2*Sw;
Llo = (CL_lo+0)*0.5*rho*0.7*Vlo^2*Sw;
Slo = (1.44*W^2)/(32.2*rho*Sref*CLmax*(T-(Dlo+mur*(W-Llo))))
slo = 1.44*W^2/(32.2*rho*Sref*CLmax*T);

Vtd = 1.3*Vstall;
Dtd = CD_td*0.5*rho*0.7*Vtd^2*Sw;
Ltd = CL_td*0.5*rho*0.7*Vtd^2*Sw;
Std = (1.69*W^2)/(32.2*rho*Sref*CLmax*(Tr+(Dtd+mur*(W-Ltd))))
LC = CL_cr*0.5*rho*V^2*Sw
DC = CD_cr*0.5*rho*V^2*Sw
Thust_req = DC/2*16 %Thrust required per motor in oz