function [npar, Vpar, Vspos, Vsneg,Va,Vg] = VNparabola(W,S,Clmaxes, nmax, nmin,rho)
Clmax = Clmaxes(1); Clmaxneg = Clmaxes(2);
Vspos = sqrt((2*W)/(rho*1.25*Clmax*S)); %Positive Stall speed 
Vsneg = sqrt((2*W)/(rho*1.25*Clmaxneg*S));%Negative stall speed 
Va = Vspos*sqrt(nmax); %manuvering speed 
Vg = sqrt((2*nmax*W)/(rho*Clmaxneg*S));%negative manuver speed 

Vpar = zeros(4,250);
npar = zeros(4,250);
Vpar(1,:) = linspace(0,Vspos,250);%Velocities for line OA. separate in order to plot this section --
Vpar(2,:) = linspace(Vspos,Va,250);%Velocities for line from Vstall+ to A
Vpar(3,:) = linspace(0,Vsneg,250);%Velocities for line OG. separate in order to plot this section --
Vpar(4,:)= linspace(Vsneg,Vg,250);%Velocities for line from Vstall- to G
npar(1,:) = (1.25*Clmax*rho*Vpar(1,:).^2*S)./(2*W); %loading from o to stall
npar(3,:) = -(1.25*Clmaxneg*rho*Vpar(3,:).^2*S)./(2*W);

%Here i iterate the loadings in order to ensure the line does not go past
%the maximum loading 
for i = 1:length(Vpar(1,:))
    npar(2,i) = (1.25*Clmax*rho*Vpar(2,i)^2*S)/(2*W);
    if npar(2,i)>= nmax
        npar(2,i) = nmax;
        Vpar(2,i) = Va;
    end
end
for i = 1:length(Vpar(4,:))
    npar(4,i) = -(1.25*Clmaxneg*rho*Vpar(4,i)^2*S)/(2*W);
    if npar(4,i)<= nmin
        npar(4,i) = nmin;
        Vpar(4,i) = Vg;
    end
end
