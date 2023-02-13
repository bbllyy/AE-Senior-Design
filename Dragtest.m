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
ch = 0.66;
Sh = bh*ch;
ARh = bh^2/Sh;

bv = 1+(1/3);
cvr = 12/12;
cvt = 8/12;
cv = (cvr+cvt)/2;
Sv = ((cvr+cvt)/2)*bv;
ARv = bv^2/Sv;

%Wing + Htail characteristics
CLmax = 1.41;
CLa = 4.894; %lift slope per rad;
CLadeg = 0.085416; %Lift slope per degree
alpha0 = -0.0622;
alpha0_deg = -3.564;
CL0 = 0.3114;
CD0 = 0.01176;

lfuse = 22.9/12;
dfuse = 4.5/12;
Swet_fuse = 2*(lfuse*dfuse+dfuse*dfuse+lfuse*dfuse);

%Flight characteristics
rho = 0.002377;
mu = 3.737E-7;
nu = mu/rho;
V = 40*1.688;
%%
aoa = 0;
CL = CLadeg*(aoa-alpha0_deg);

% 
if tc_max >= 0.05
    K = 1.9767+0.5333*tc_max;
else 
    K = 2.0;
end
if tc_max >= 0.3
    L = 1.2;
else 
    L = 2.0;
end
c = [cw ch cv lfuse];
S = [Sw Sh Sv];

Swet = [(S.*K) Swet_fuse];
ld = lfuse/dfuse;
RE = (V.*c)./nu;
CLw = 1.05*CL;
Cf = 0.455./(log10(RE).^2.58);
Rls = 1.075;
Rwf = 1; %Setting as 1 because RE is so low this doesnt matter

CD0_foils = Rwf*Rls.*Cf(1:3).*(1+L*tc_max+100*tc_max^4).*(Swet(1:3)./Sw);
CD0_fuse = Cf(4)*(1+(60/ld^3)+0.0025*ld)*(Swet_fuse/Sw);
CDi = CLw^2/(pi*e*AR);
CDbase = 0.05;
CDgear = 0.005;

CD = sum(CD0_foils)+CD0_fuse+sum(CDi)+CDbase+CDgear