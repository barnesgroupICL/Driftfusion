function C_MS = C_MottSchottky(eppr, Nion, V, par)
epp0 = par.epp0;
q = par.e;

C_MS =((q^2*eppr*epp0*Nion)/(2*V));


end
