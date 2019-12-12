initialise_df

% Create parameters objects for Spiro/MAPI/TiO2 and PEDOT/MAPI/PCBM devices
% par = pc('Input_files/spiro_mapi_tio2_simple.csv');
par = pc('Input_files/pn_junction_DA_compare.csv');

soleq = equilibrate_short(par);
dfplot.ELx(soleq);

%JV = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
JV = doJV_short(soleq, 1, 241, 1, 0, 0, 1.1, 3);