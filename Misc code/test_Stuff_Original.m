<<<<<<< HEAD
%This code is meant to help guide someone through the process of
%preliminary design. Each calculation should be checked again and verified
%after initial design is finished. What this code does is start with a
%basic constraint analysis using an input AR and Weight. The user chooses a
%W/S (T/W is neglected) from which other wing parameters are
%chosen/calculated. Once the run is complete, rerun the code with the newly
%calculated AR value and cycle until desirable values are recieved.

%It is best to have some ideas about the general size of your aircraft.
%Basic ideas include max wingspan, root chord and weight. 


clc
clear all
close all

%Initial User Inputs
AR = 6.56;                    %Aspect Ratio
W = 15;                       %Weight (lbs)
S = 6.57;                   %Wing area

%Flight Parameters 
P_SL = 2116.22;              %Standard Day Sea Level Pressure (lb/ft^2)
rho_SL = 0.002377;           %Standard Day Sea Level Density (lb/ft^2)
g0_SL = 32.174;              %Sea Level Gravity (lb/ft^2)
del_SL = 1.00;               %Ratio of P/Pstd
theta_SL = 1.00;             %Ratio of T/Tstd (standard day)
sigma = del_SL/theta_SL;  %del/sigma ratio of ratios

%Weight Specific Excess Power
h = (0:10:5000);
V = (0:1:250);
z_e = zeros(size(h,2),size(V,2));
for i = 1:size(h,2)
    z_e(i,:) = h(i) + (V.^2)/(2*g0_SL);
end

%Varying Values
WTO_S = (1:0.1:10);              %Take off weight/wing area ratio range (lbs/ft^2)
%Master Equation

%Aircraft Parameters
e = 0.92;                       %Oswald Efficiency Factor
CD_0 = 0.015;                   %Estimate for min drag coefficient
CD_R = 0;
CL_min = 0.35;                     %CL at alpha = 0
CL_max_set = 1.4;                 %Set CL_max
%CL_Vmax = 0.1;
K1 = 1/(e*pi*AR);               %K1 factoring in CD = K1*CL^2 + K2*CL + CD0
K2 = 0;                         %K2 factoring in CD = K1*CL^2 + K2*CL + CD0
Vstall = sqrt(2*W/(CL_max_set*rho_SL*S));                    %V stall in fps
%Vstall = 53;
WTO_S_stall = CL_max_set*.5*rho_SL*Vstall^2;
WTO_S_stall = [WTO_S_stall WTO_S_stall];
TSL_WTO_stall = [0 2];

%% Case 1: Constant Altitude/Speed Cruise (Ps=0)
%dh_dt = 0; %dV_dt = 0; %n = 1;
%Given h and V
beta = 0.9;
alpha = 0.75;
Vset = 50*1.69;                %Max Cruise Speed (conversion of knots to fps)
qset = 0.5*rho_SL*Vset.^2;     %Max Cruise Dynamic Pressure

WTO_S_minTW_c1 = (qset/beta)*sqrt((CD_0+CD_R)/K1);
TSL_WTO_min_c1 = (beta/alpha)*(2*sqrt((CD_0+CD_R)*K1)+K2);
TSL_WTO_c1 = (beta/alpha)*(K1*(beta/qset)*(WTO_S) + K2 + ((CD_0+CD_R)./((beta/qset)*(WTO_S))));

%% Case 2: Constant Max Speed
beta = 0.9;
alpha = 1.0;
Vset = 50*1.69;                %Max Cruise Speed (conversion of knots to fps)
qset = 0.5*rho_SL*Vset.^2;     %Max Cruise Dynamic Pressure

WTO_S_minTW_c15 = (qset/beta)*sqrt((CD_0+CD_R)/K1);
TSL_WTO_min_c15 = (beta/alpha)*(2*sqrt((CD_0+CD_R)*K1)+K2);
TSL_WTO_c15 = (beta/alpha)*(K1*(beta/qset)*(WTO_S) + K2 + ((CD_0+CD_R)./((beta/qset)*(WTO_S))));
%% Case 3: Constant Speed Climb (Ps = dh/dt)
%dV_dt = 0; %n = 1;
%Given h, dh_dt, and V
beta = 0.9;
alpha = 0.8;
Vset = 50*1.69;                %Max Cruise Speed (conversion of knots to fps)
qset = 0.5*rho_SL*Vset.^2;     %Max Cruise Dynamic Pressure

dh_dt = 10;             %Climb (ft/s)
WTO_S_minTW_c2 = (qset/beta)*sqrt((CD_0+CD_R)/K1);
TSL_WTO_min_c2 = (beta/alpha)*(2*sqrt((CD_0+CD_R)*K1)+K2+(1/Vset)*dh_dt);
TSL_WTO_c2 = (beta/alpha)*(K1*(beta/qset)*(WTO_S) + K2 + ((CD_0+CD_R)./((beta/qset)*(WTO_S))) + (1/Vset)*dh_dt);

%% Case 4: Constant Altitude/Speed Turn (Ps = 0)
%dh_dt = 0; %dV_dt = 0;
%Given h and V
beta = 0.9;
alpha = 0.8;
Vset = 50*1.69;                %Max Cruise Speed
qset = 0.5*rho_SL*Vset.^2;     %Max Cruise Dynamic Pressure

Rc = 200;                                   %Turn Radius (ft)
n = sqrt(1+((Vset^2)/(g0_SL*Rc))^2);     %Load factor based on turn radius
WTO_S_minTW_c3 = (qset/(n*beta))*sqrt((CD_0+CD_R)/K1);
TSL_WTO_min_c3 = (n*beta/alpha)*(2*sqrt((CD_0+CD_R)*K1)+K2+(1/Vset)*dh_dt);
TSL_WTO_c3 = (beta/alpha)*(K1*(n^2)*(beta/qset)*(WTO_S) + K2*n + ((CD_0+CD_R)./((beta/qset)*(WTO_S))));

%% Constraint Analysis Graph
figure
plot(WTO_S,TSL_WTO_c1, WTO_S,TSL_WTO_c15, WTO_S,TSL_WTO_c2, WTO_S,TSL_WTO_c3, WTO_S_stall,TSL_WTO_stall,'LineWidth',2) 
xlabel('W_T_O/S (lbs/ft^2)');
ylabel('T_S_L/W_T_O');
title('Constraint Analysis - All Cases');
legend({'Constant V/H Cruise: 100 knots','Constant Max Speed: 125 knots','Constant V Climb: 10 ft/s','Constant V/H Turn: R = 200 ft','Clean Stall: 25.5 knots'},'Fontsize',14)
axis([1 10 0 2]);
grid on

[Wing_Loading,~] = ginput(1);
%% Wing Sizing
b_W = (.5:.1:10);    %wing span variation in ft
S_W = W/Wing_Loading;   %Wing Area
c_bar_W = S_W./b_W;   %Mean wing chord in ft

plot(c_bar_W,b_W);
title('Select a Wingspan & Mean Chord Length')
ylabel('Span (ft)')
xlabel('Mean Chord (ft)')
[c_bar_W,b_W] = ginput(1);

AR_W = b_W^2/S_W;   %Aspect ratio from user input

cr_W = (0:.01:10);   %root chord variation in ft
ct_W = 2*c_bar_W-cr_W;   %tip chord in ft
lambda_W = ct_W./cr_W;  %taper ratio 

plot(cr_W,lambda_W)
title('Select a Root Chord and Taper Ratio')
ylabel('Taper Ratio')
xlabel('Root Chord (ft)')
ylim([0,1])
[cr_W,lambda_W] = ginput(1);    %root chord & taper ratio based off user input
ct_W = 2*c_bar_W-cr_W;   %tip chord based off user input

Vstall = sqrt(2*Wing_Loading/(rho_SL*CL_max_set));   %Vstall in fps
Vcruise = sqrt(2*Wing_Loading/(rho_SL*CL_min));    %Vcruise in fps

%% Wing Planform

%This section calculates some wing characteristics based on your choices
%above. This assumes that your wing planform is "traditional"and not swept. It will
%allow for the tail sizing to be done relatively well for your planform.
%Again this code is only meant for preliminary data only

%Straight
if (lambda_W == 1)
    LE_LAMBDA = 0;
%Leadin/Trailing/Rectangular    
elseif (1>lambda_W)&&(lambda_W >=.6)
    prompt = 'Is the Wing a [Select a Number] :(1)Leading Edge, (2)Trailing Edge, (3)Rectangular  ';
    x = input(prompt);
    if x == 1
        LE_LAMBDA = tan((cr_W-ct_W)/(b_W/2));
    elseif x == 2
        LE_LAMBDA = 0;
    elseif x == 3
        LE_LAMBDA = tan(((cr_W-ct_W)/2)/(b_W/2));
    end
%Elliptical    
elseif (0.4 <= lambda_W)&&(lambda_W < .6)
    LE_LAMBDA = tan(((cr_W-ct_W)/2)/(b_W/2));
%Delta    
elseif (lambda_W < 0.4)
    LE_LAMBDA = tan((cr_W-ct_W)/(b_W/2));
end

x_LE_W = b_W - (b_W*.7); %Distance from nose to root leading edge assumed to be 30% of the wingspan. Can be edited
x_mac = b_W/6*(1 + 2*lambda_W)/(1 + lambda_W)*tan(LE_LAMBDA);    %MAC in ft
x_ac_W = x_LE_W + x_mac + 0.4*c_bar_W;    %Location of the AC in ft
x_cg = x_LE_W + x_mac + 0.2*c_bar_W;   %Location of cg in ft
CL_alpha_W = 2*pi/(sqrt(1 + (2*pi/(pi*AR_W).^2)) + 2*pi/(pi*AR_W));
%% Horizontal Tail Sizing

%This section is more of a guide on how to do the calculations. The tail
%volume coefficient is set at 0.4 as a safe estimate. From this, the static
%margin can be used to find the optimal span and area for the horizontal.

V_h = 0.4;   %0.4 is a safe estimate but can be modified
AR_h = (0:0.1:5);
CL_alpha_h = 2*pi./(sqrt(1 + (2*pi./(pi*AR_h).^2)) + 2*pi./(pi.*AR_h));
NP = x_ac_W - V_h.*CL_alpha_h.*c_bar_W./CL_alpha_W;
sm = (NP - x_cg)./c_bar_W*100;

plot(AR_h,sm)
xlabel('Aspect Ratio for Horizontal')
ylabel('Static Margin Associated with AR')
title('Aspect Ratio vs Static Margin')
[AR_h,sm] = ginput(1);

S_h = (0:0.1:S_W);
l_h = V_h*S_W*c_bar_W./S_h;
b_h = sqrt(AR_h.*S_h);   

yyaxis left
plot(l_h,S_h)
xlim([0,3])
ylabel('Horizontal Area in ft^2')
yyaxis right
plot(l_h,b_h)
ylabel('Horizontal Span in ft')
xlabel('Distance Between 1/4 Chords in ft')
title('Select Area & Distance to Tail')

[l_h,S_h] = ginput(1);

b_h = sqrt(AR_h.*S_h);   

=======
%This code is meant to help guide someone through the process of
%preliminary design. Each calculation should be checked again and verified
%after initial design is finished. What this code does is start with a
%basic constraint analysis using an input AR and Weight. The user chooses a
%W/S (T/W is neglected) from which other wing parameters are
%chosen/calculated. Once the run is complete, rerun the code with the newly
%calculated AR value and cycle until desirable values are recieved.

%It is best to have some ideas about the general size of your aircraft.
%Basic ideas include max wingspan, root chord and weight. 


clc
clear all
close all

%Initial User Inputs
AR = 6.56;                    %Aspect Ratio
W = 15;                       %Weight (lbs)
S = 6.57;                   %Wing area

%Flight Parameters 
P_SL = 2116.22;              %Standard Day Sea Level Pressure (lb/ft^2)
rho_SL = 0.002377;           %Standard Day Sea Level Density (lb/ft^2)
g0_SL = 32.174;              %Sea Level Gravity (lb/ft^2)
del_SL = 1.00;               %Ratio of P/Pstd
theta_SL = 1.00;             %Ratio of T/Tstd (standard day)
sigma = del_SL/theta_SL;  %del/sigma ratio of ratios

%Weight Specific Excess Power
h = (0:10:5000);
V = (0:1:250);
z_e = zeros(size(h,2),size(V,2));
for i = 1:size(h,2)
    z_e(i,:) = h(i) + (V.^2)/(2*g0_SL);
end

%Varying Values
WTO_S = (1:0.1:10);              %Take off weight/wing area ratio range (lbs/ft^2)
%Master Equation

%Aircraft Parameters
e = 0.92;                       %Oswald Efficiency Factor
CD_0 = 0.015;                   %Estimate for min drag coefficient
CD_R = 0;
CL_min = 0.35;                     %CL at alpha = 0
CL_max_set = 1.4;                 %Set CL_max
%CL_Vmax = 0.1;
K1 = 1/(e*pi*AR);               %K1 factoring in CD = K1*CL^2 + K2*CL + CD0
K2 = 0;                         %K2 factoring in CD = K1*CL^2 + K2*CL + CD0
Vstall = sqrt(2*W/(CL_max_set*rho_SL*S));                    %V stall in fps
%Vstall = 53;
WTO_S_stall = CL_max_set*.5*rho_SL*Vstall^2;
WTO_S_stall = [WTO_S_stall WTO_S_stall];
TSL_WTO_stall = [0 2];

%% Case 1: Constant Altitude/Speed Cruise (Ps=0)
%dh_dt = 0; %dV_dt = 0; %n = 1;
%Given h and V
beta = 0.9;
alpha = 0.75;
Vset = 50*1.69;                %Max Cruise Speed (conversion of knots to fps)
qset = 0.5*rho_SL*Vset.^2;     %Max Cruise Dynamic Pressure

WTO_S_minTW_c1 = (qset/beta)*sqrt((CD_0+CD_R)/K1);
TSL_WTO_min_c1 = (beta/alpha)*(2*sqrt((CD_0+CD_R)*K1)+K2);
TSL_WTO_c1 = (beta/alpha)*(K1*(beta/qset)*(WTO_S) + K2 + ((CD_0+CD_R)./((beta/qset)*(WTO_S))));

%% Case 2: Constant Max Speed
beta = 0.9;
alpha = 1.0;
Vset = 50*1.69;                %Max Cruise Speed (conversion of knots to fps)
qset = 0.5*rho_SL*Vset.^2;     %Max Cruise Dynamic Pressure

WTO_S_minTW_c15 = (qset/beta)*sqrt((CD_0+CD_R)/K1);
TSL_WTO_min_c15 = (beta/alpha)*(2*sqrt((CD_0+CD_R)*K1)+K2);
TSL_WTO_c15 = (beta/alpha)*(K1*(beta/qset)*(WTO_S) + K2 + ((CD_0+CD_R)./((beta/qset)*(WTO_S))));
%% Case 3: Constant Speed Climb (Ps = dh/dt)
%dV_dt = 0; %n = 1;
%Given h, dh_dt, and V
beta = 0.9;
alpha = 0.8;
Vset = 50*1.69;                %Max Cruise Speed (conversion of knots to fps)
qset = 0.5*rho_SL*Vset.^2;     %Max Cruise Dynamic Pressure

dh_dt = 10;             %Climb (ft/s)
WTO_S_minTW_c2 = (qset/beta)*sqrt((CD_0+CD_R)/K1);
TSL_WTO_min_c2 = (beta/alpha)*(2*sqrt((CD_0+CD_R)*K1)+K2+(1/Vset)*dh_dt);
TSL_WTO_c2 = (beta/alpha)*(K1*(beta/qset)*(WTO_S) + K2 + ((CD_0+CD_R)./((beta/qset)*(WTO_S))) + (1/Vset)*dh_dt);

%% Case 4: Constant Altitude/Speed Turn (Ps = 0)
%dh_dt = 0; %dV_dt = 0;
%Given h and V
beta = 0.9;
alpha = 0.8;
Vset = 50*1.69;                %Max Cruise Speed
qset = 0.5*rho_SL*Vset.^2;     %Max Cruise Dynamic Pressure

Rc = 200;                                   %Turn Radius (ft)
n = sqrt(1+((Vset^2)/(g0_SL*Rc))^2);     %Load factor based on turn radius
WTO_S_minTW_c3 = (qset/(n*beta))*sqrt((CD_0+CD_R)/K1);
TSL_WTO_min_c3 = (n*beta/alpha)*(2*sqrt((CD_0+CD_R)*K1)+K2+(1/Vset)*dh_dt);
TSL_WTO_c3 = (beta/alpha)*(K1*(n^2)*(beta/qset)*(WTO_S) + K2*n + ((CD_0+CD_R)./((beta/qset)*(WTO_S))));

%% Constraint Analysis Graph
figure
plot(WTO_S,TSL_WTO_c1, WTO_S,TSL_WTO_c15, WTO_S,TSL_WTO_c2, WTO_S,TSL_WTO_c3, WTO_S_stall,TSL_WTO_stall,'LineWidth',2) 
xlabel('W_T_O/S (lbs/ft^2)');
ylabel('T_S_L/W_T_O');
title('Constraint Analysis - All Cases');
legend({'Constant V/H Cruise: 100 knots','Constant Max Speed: 125 knots','Constant V Climb: 10 ft/s','Constant V/H Turn: R = 200 ft','Clean Stall: 25.5 knots'},'Fontsize',14)
axis([1 10 0 2]);
grid on

[Wing_Loading,~] = ginput(1);
%% Wing Sizing
b_W = (.5:.1:10);    %wing span variation in ft
S_W = W/Wing_Loading;   %Wing Area
c_bar_W = S_W./b_W;   %Mean wing chord in ft

plot(c_bar_W,b_W);
title('Select a Wingspan & Mean Chord Length')
ylabel('Span (ft)')
xlabel('Mean Chord (ft)')
[c_bar_W,b_W] = ginput(1);

AR_W = b_W^2/S_W;   %Aspect ratio from user input

cr_W = (0:.01:10);   %root chord variation in ft
ct_W = 2*c_bar_W-cr_W;   %tip chord in ft
lambda_W = ct_W./cr_W;  %taper ratio 

plot(cr_W,lambda_W)
title('Select a Root Chord and Taper Ratio')
ylabel('Taper Ratio')
xlabel('Root Chord (ft)')
ylim([0,1])
[cr_W,lambda_W] = ginput(1);    %root chord & taper ratio based off user input
ct_W = 2*c_bar_W-cr_W;   %tip chord based off user input

Vstall = sqrt(2*Wing_Loading/(rho_SL*CL_max_set));   %Vstall in fps
Vcruise = sqrt(2*Wing_Loading/(rho_SL*CL_min));    %Vcruise in fps

%% Wing Planform

%This section calculates some wing characteristics based on your choices
%above. This assumes that your wing planform is "traditional"and not swept. It will
%allow for the tail sizing to be done relatively well for your planform.
%Again this code is only meant for preliminary data only

%Straight
if (lambda_W == 1)
    LE_LAMBDA = 0;
%Leadin/Trailing/Rectangular    
elseif (1>lambda_W)&&(lambda_W >=.6)
    prompt = 'Is the Wing a [Select a Number] :(1)Leading Edge, (2)Trailing Edge, (3)Rectangular  ';
    x = input(prompt);
    if x == 1
        LE_LAMBDA = tan((cr_W-ct_W)/(b_W/2));
    elseif x == 2
        LE_LAMBDA = 0;
    elseif x == 3
        LE_LAMBDA = tan(((cr_W-ct_W)/2)/(b_W/2));
    end
%Elliptical    
elseif (0.4 <= lambda_W)&&(lambda_W < .6)
    LE_LAMBDA = tan(((cr_W-ct_W)/2)/(b_W/2));
%Delta    
elseif (lambda_W < 0.4)
    LE_LAMBDA = tan((cr_W-ct_W)/(b_W/2));
end

x_LE_W = b_W - (b_W*.7); %Distance from nose to root leading edge assumed to be 30% of the wingspan. Can be edited
x_mac = b_W/6*(1 + 2*lambda_W)/(1 + lambda_W)*tan(LE_LAMBDA);    %MAC in ft
x_ac_W = x_LE_W + x_mac + 0.4*c_bar_W;    %Location of the AC in ft
x_cg = x_LE_W + x_mac + 0.2*c_bar_W;   %Location of cg in ft
CL_alpha_W = 2*pi/(sqrt(1 + (2*pi/(pi*AR_W).^2)) + 2*pi/(pi*AR_W));
%% Horizontal Tail Sizing

%This section is more of a guide on how to do the calculations. The tail
%volume coefficient is set at 0.4 as a safe estimate. From this, the static
%margin can be used to find the optimal span and area for the horizontal.

V_h = 0.4;   %0.4 is a safe estimate but can be modified
AR_h = (0:0.1:5);
CL_alpha_h = 2*pi./(sqrt(1 + (2*pi./(pi*AR_h).^2)) + 2*pi./(pi.*AR_h));
NP = x_ac_W - V_h.*CL_alpha_h.*c_bar_W./CL_alpha_W;
sm = (NP - x_cg)./c_bar_W*100;

plot(AR_h,sm)
xlabel('Aspect Ratio for Horizontal')
ylabel('Static Margin Associated with AR')
title('Aspect Ratio vs Static Margin')
[AR_h,sm] = ginput(1);

S_h = (0:0.1:S_W);
l_h = V_h*S_W*c_bar_W./S_h;
b_h = sqrt(AR_h.*S_h);   

yyaxis left
plot(l_h,S_h)
xlim([0,3])
ylabel('Horizontal Area in ft^2')
yyaxis right
plot(l_h,b_h)
ylabel('Horizontal Span in ft')
xlabel('Distance Between 1/4 Chords in ft')
title('Select Area & Distance to Tail')

[l_h,S_h] = ginput(1);

b_h = sqrt(AR_h.*S_h);   

>>>>>>> ff4d41e2726217651ffaa74f2f21343ad48b8382
