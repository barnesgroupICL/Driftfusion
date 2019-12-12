% An example script for comparing the current-voltage curves for two different
% architects of perovskite solar cell
initialise_df

% Create parameters objects for Spiro/MAPI/TiO2 and PEDOT/MAPI/PCBM devices
% par = pc('Input_files/spiro_mapi_tio2_simple.csv');
par = pc('Input_files/3_layer_methods_test1.csv');
% par = pc('Input_files/pn_heterojunction.csv');
% Find equilibrium solutions
soleq = equilibrate(par);
soleq.no_ion.par.gx1 = ASA_gentot_DFmesh;
% Perform current-voltage scans to 1.4 V
%JV = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
JV_dk = doJV(soleq.no_ion, 1e-9, 281, 0, 0, 0, 1.4, 1);
JV_ill = doJV(soleq.no_ion, 1e-3, 281, 1, 0, 0, 1.4, 2);

% plot the current voltage curve
dfplot.JV(JV_ill,2)
%set(gca, 'YScale', 'log')
% legend('tio2-dk-f', 'tio2-dk-r', 'tio2-1sun-f', 'tio2-1sun-r')
% xlim([0, 1.2])
% ylim([-30e-3, 10e-3])
%hold off

% plot the energy level diagram and carrier densities for the tio2 device at
% 1 V (t= 20s) during the illuminated forward scan
% dfplot.ELx(JV.kscan50em3.ill.f, 20)
% 
% % plot the currents as a function of position in the PCBM device at 0.5 V
% % (t = 10 s) during the illuminated forward scan
% dfplot.Jx(JV_pcbm.kscan50em3.ill.f, 10)

% Save the workspace- this is commented out as the filepath should lead to
% a folder on your computer. It is not recommended to store large files in
% the Github repository as they will be cached even when deleted leading to
% a large folder sizes
% save('/Users/Username/Data/temp.mat')
