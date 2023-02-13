clc, clear, close all

m_bat = (729*2)/1000;
m_plane = ((15*453.6)-m_bat)/1000;
V_bat = 14.21; %Voltage
A_bat = 16; %Capacity in AH
eta_prop = 0.8;
n_prop = 6;
rho_bat = (V_bat*A_bat)/m_bat;
E_bat = rho_bat*m_bat;
W_motor = 261*2;
time_h = E_bat/W_motor
time_minute = time_h*60
% 
% time_h = ((rho_bat*m_bat)/(m_bat+m_plane))*eta_prop*((m_bat+m_plane)/n_prop)
% time_minute = time_h*60