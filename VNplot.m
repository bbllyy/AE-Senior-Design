function VNplot=VNplot(Vpar, npar, Va, Vc, Vd, Vg, Vspos, Vsneg, nmax, nmin, Vgust, ngust,n2)

A = [Va nmax];
C = [Vc nmax];
F = [Vd n2];
G = [Vg nmin];
D = [Vd nmin];
Voa = Vpar(1,:); Vsa = Vpar(2,:); Vog = Vpar(3,:); Vsg = Vpar(4,:);
noa = npar(1,:); nsa = npar(2,:); nog = npar(3,:); nsg = npar(4,:);
ngcpos = ngust(1,:); ngcneg = ngust(2,:); ngdpos = ngust(3,:); ngdneg = ngust(4,:);
figure
VNplot = plot(Vsa, nsa,'k');
hold on
plot(Voa,noa,'--k');plot(Vog,nog,'--k');plot(Vsg, nsg,'k') %Plots for the parabolic curves
nproofp = 1.25*nmax; nproofn = 1.25*nmin;
nultp = 1.5*nmax; nultn = 1.5*nmin;

plot([A(1) C(1)], [A(2) C(2)],'k') %Line AC
plot([C(1) F(1)], [C(2) F(2)],'k') %Line CF
plot([F(1) D(1)], [F(2) D(2)],'k') %Line FD
plot([G(1) D(1)], [G(2) D(2)],'k') %Line GD
plot([0 Vd], [0 0],'k') %Middle line
plot([0 Vspos],[1 1],'--k');plot([0 Vsneg],[-1 -1],'--k') %lines at n=1
pgust1 = plot(Vgust,ngcpos, '--b','LineWidth', 0.1);plot(Vgust,ngcneg, '--b','LineWidth', 0.1) %plots for the Vc gust lines
pgust2 = plot(Vgust,ngdpos, '--r','LineWidth', 0.1);plot(Vgust,ngdneg, '--r','LineWidth', 0.1)%plots for the Vd gust lines
plot([Vspos Vspos], [1 0], 'k');plot([Vsneg Vsneg], [-1 0], 'k') %plots for the stall lines
xline(Va, '--k','Va'); xline(Vc, '--k','Vc') %lines to label Va and Vc
yline(nmax, '--k','Max Load');yline(nmin, '--k','min Load') %lines to label max and min load
yline(nproofp, '--k','proof Load');yline(nultp, '--k','ultimate Load') %lines to label 
yline(nproofn, '--k','proof Load');yline(nultn, '--k','ultimate Load') %lines to label max and min load
axis([0 (Vd+10) (nultn-2) (nultp+1)])
yticks([-6 -4 -2 -1 0 1 2 4 6 8])
title('V-n Diagram')
ylabel('Load Factor n'); xlabel('V(m/s)')
legend([pgust1 pgust2],{'\pmV_c gust-line', '\pmV_d gust-line'},'Location', 'southwest')