 clc, clear, close all

g = 32.2;
foil = 4412;
%Can either specify AR, or do span/chord
b = 8;
c = 1;
S = b*c;
AR = b^2/S;%6;
e = 0.90;
rho = 0.002377;
sigma = 1; %rho/rho_sl
mu = 3.737E-7/32.2;
nu = mu/rho;
s_to = 60;
s_l = 100;
V = 40*1.688;
E_motor = 600*2;
W_bat = 2*2; %Weight for 2 batteries
bat_E = 16*11.1; %Battery energy in W*h
Ebat = bat_E*3600; %Battery energy in Joules
K_bat = Ebat/W_bat; %Battery K Value
Ws  = 0.03:0.01:5;
RC = 5;
q = 0.5*rho*V^2;
n = 1.15;
Wreal = 15; %With Fudge Factor
AR_real = b^2/S;
WS_real = Wreal/S;
W_land = Wreal-2.5;
if foil == 2412
    CLmax = 1.272;
    CD0 = 0.0075;
    CL_a = 4.037; %CL_alpha in /rad
    CL0 = 0.173;
    alpha0 = -2.25; %AOA at 0 lift in deg
elseif foil == 4412
    CLmax = 1.398;
    CD0 = 0.0101;
    CL_a = 4.4803; %CL_alpha in /rad
    CL0 = 0.35;
    alpha0 = -4.34; %AOA at 0 lift in deg
else 
    CLmax = 1.2;
    CD0 = 0.08;
end
%%  Constraint Analysis
Vs_real = sqrt((2*Wreal)/(rho*S*CLmax));
WTO_S_stall = CLmax*.5*rho*1.2*Vs_real^2;
Vs_land = sqrt((2*W_land)/(rho*S*CLmax));
etap = 0.75;%prop eff
etam = 0.85; %motor eff
alpha = 1;
CLmaxLD = sqrt(CD0*e*pi*AR);
K_clcd = 1/(pi*e*AR);
LDmax = 1/(2*sqrt(K_clcd*CD0));

%Takeoff

Vstall = ((2*Ws)./(sigma*rho*CLmax)).^0.5;
Vto = 1.2*Vstall;
PW_to = ((0.7*Vto).^3*746)./(2*550*alpha*etap*etam*g*s_to);
PW_land = (Vto.^3.*746)./(550*alpha*etap*etam*g*s_l);
Cl = sqrt((3*CD0)/K_clcd);
Vy = ((2*Ws)./(rho*Cl)).^0.5;
PW_ceiling = (Vy*746)./(0.866*550*alpha*etap*etam*g*(sigma^0.5));
Vminpower = ((2*Ws)./(rho*Cl)).^0.5;
PW_rc = (746./(550.*alpha.*etap.*etam)).*(RC+(Vminpower./(0.866.*LDmax.*sigma^0.5)));
TW = (1/alpha)*(((q*CD0)./Ws)+((K_clcd.*Ws)./(q)));
PW_maxV = (V.*TW*746)./(550*etap*etam);
TW_turn = (1/alpha)*(((q*CD0)./Ws)+(n^2*(K_clcd.*Ws)./q));
PW_turn = (V.*TW_turn*746)./(550*etap*etam);
stall_speed = sprintf('Clean Stall %2.2f ft/s',Vs_real);
figure
plot(Ws,PW_to)
hold on
plot(Ws,PW_land)
plot(Ws,PW_ceiling)
plot(Ws,PW_rc)
plot(Ws,PW_maxV)
plot(Ws,PW_turn)
scatter(WS_real,80,"or")
xline(WTO_S_stall,'-',{stall_speed},'LineWidth',1.2)
legend('Takeoff','Landing','Ceiling','Rate of Climb','Max V','Turn', 'Design Point')
hold off
grid on
title('Constraint Analysis')
ylim([0 300])
%xlim([0 1.5])
ylabel('P/W (W/lbs)')
xlabel('W/S (lbs/ft^2)')
%%
%Aircraft Caluclations
LDmax = ((pi*e*AR)/(4*CD0))^0.5;
LD32 = 0.25*((3*pi*e*AR)/CD0^(1/3))^(3/4);

Motor_V = 14;
Motor_I = 24;
Motor_W = Motor_V*Motor_I;
Flight_time_full = (bat_E/(Motor_W*2))*60*etam;
Range_full = 50*(Flight_time_full/60);
motor_V_80 = 10.9;
motor_I_80 = 21.5;
Motor_W_80 = motor_I_80*motor_V_80;
flight_time_80 = (bat_E/(Motor_W_80*2))*60*etam
Range_80 = 45*(flight_time_80/60)
%%
%Thrust per speed
density_wings = 35;
density_fiber = 1.5;
density_epoxy = 1100;
density_cf = 1750.0;
density_wood = 600;
rho_wing_tot = density_wings+density_fiber;
Vwing = 0.01835351;
Vtail = 0.00480141;
V_epoxy_fuse = 0.000709765;%24 fl oz
V_epoxy_tail = 0.000177441; %6 fl oz
V_epoxy_wing = 0.000354882; %12 fl oz
Vfuselage = 0.00131097;
V_spar = 0.00034822511;
V_ribs = 0.000737418;
Mass_wing = Vwing*rho_wing_tot;
Mass_tail = Vtail*rho_wing_tot;
Mass_fuse = Vfuselage*density_fiber;
Mass_wing_epoxy = V_epoxy_wing*density_epoxy;
Mass_fuse_epoxy = V_epoxy_fuse*density_epoxy;
mass_spar = V_spar*density_cf;
mass_wood = V_ribs*density_wood;
weight_wing = (Mass_wing+Mass_wing_epoxy)*2.2;
Mass_tail_epoxy = V_epoxy_tail*density_epoxy;
weight_tail = (Mass_tail+Mass_tail_epoxy)*2.2;
weight_fuse = (Mass_fuse+Mass_fuse_epoxy)*2.2;
weight_spar = mass_spar*2.2;
weight_wood = mass_wood*2.2;
weight_machine = 1;
weight = weight_spar+weight_fuse+weight_tail+weight_wing+weight_wood;

%%
T = 6.225;
W = 13;
Sto = (1.44*W^2)/(g*rho*CLmax*S*T);
