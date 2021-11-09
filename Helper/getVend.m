function Vend = getVend(sol)
% Obtains the applied potential VEND at the final time point of SOL

Vappt = dfana.calcVapp(sol);
Vend = Vappt(end);

end