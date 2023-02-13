function [Vgust,ngust,mug] = VNgust(W,S,rho,c,Cla,Vd,Udec,Uded)
const = 498;

mug = (2*(W/S))/(rho*c*32.2*Cla);
Kg = (0.88*mug)/(5.3+mug);
Vgust = linspace(0,Vd+100,500);
ngust = zeros(4,500);
ngust(1,:) =1+((Kg*Udec*Vgust*Cla)/(const*(W/S)));%gust line +c load
ngust(2,:) = 1-((Kg*Udec*Vgust*Cla)/(const*(W/S)));%gust line -c load
ngust(3,:) =1+((Kg*Uded*Vgust*Cla)/(const*(W/S)));%gust line +d load
ngust(4,:) = 1-((Kg*Uded*Vgust*Cla)/(const*(W/S)));%gust line -d load
end