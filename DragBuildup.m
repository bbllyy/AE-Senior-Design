function [CL,CD] = DragBuildup(aoa, PlaneCharacteristics)
load(PlaneCharacteristics)
CL = CLadeg*(aoa-alpha0_deg);
cv = (cvr+cvt)/2;
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
CDbase = 0.04;
CDgear = 0.04;

CD = sum(CD0_foils)+CD0_fuse+sum(CDi)+CDbase+CDgear;
end