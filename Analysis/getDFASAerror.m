function getDFASAerror(DFsol, ASAsol)

Vapp_DF = dfana.calcVapp(DFsol);
J_DF = dfana.calcJ(DFsol);
Jtot_DF = J_DF.tot(:, 1)';

Vapp_ASA = ASAsol(1:261, 1)';
J_ASA = ASAsol(1:261, 2)'*1e-4;

error = 100*(J_ASA-Jtot_DF)./J_ASA;

% figure(300)
% semilogy(Vapp_DF, Jtot_DF, Vapp_ASA, J_ASA, '--')
% xlabel('Voltage (V)')
% ylabel('Current (Acm-2)')

figure(301)
plot(Vapp_DF, error)
xlabel('Voltage (V)')
ylabel('Difference (%)')
ylim([-5, 5])
end