clc, clear, close all

bw = 8;
cw = 1;
Sw = bw*cw
ARw = bw^2/Sw;%6;

Vh = 0.786;
Vv = 0.062;
Slh = Vh*Sw*cw;
Sh = 0.3*Sw
Lh = Slh/Sh
lh = 3.5;
Sh2 = Slh/lh
lv = 3.5;
Slv = Vv*bw*Sw;
Sv = Slv/lv
ARh = 4;
bh = (ARh*Sh2)^0.5
Ch = Sh2/bh

arv = 2
bv = (arv*Sv)^0.5
Cv = Sv/bv
