par_tio2 = pc('Input_files/spiro_mapi_tio2.csv');
soleq_tio2 = equilibrate(par_tio2);
sol_CV_100mVs_tio2_dark = doCV(soleq_tio2.ion, 0, 0, 1.0, 0, 100e-3, 1, 241);
compare_rec_flux_general(sol_CV_100mVs_tio2_dark, 1, 0.02, 1e6);