function [sol_Voc, Voc] = findVocDirect(sol_ini, light_intensity, mobseti)
% Obtain approximate open circuit voltage directly using high Rs
sol_Voc = lightonRs(sol_ini, light_intensity, -1, mobseti, 1e6, 400);

Voct = dfana.calcVQFL(sol_Voc);
Voc = Voct(end);

end