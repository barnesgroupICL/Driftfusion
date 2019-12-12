initialise_df

par_pn = pc('Input_files/pn_junction_DA_compare.csv');

soleq_pn = equilibrate(par_pn);
dfplot.ELx(soleq);

% JV = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
JV = doJV(soleq, 1, 241, 1, 0, 0, 1.1, 3);

