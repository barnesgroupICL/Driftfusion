function [C_donor, C_acceptor] = C_MottSchottky(sol)
% Calculates Mott Schottky capacitance as a function of Vapp(t)
% Assumes single layer only
par = sol.par;
Vapp = dfana.calcVapp(sol);

C_donor = ((par.e^2*par.epp0*par.epp(1)*par.ND(1))./(2*(par.Vbi-Vapp))).^0.5;
C_acceptor = ((par.e^2*par.epp0*par.epp(1)*par.NA(1))./(2*(par.Vbi-Vapp))).^0.5;

end