%par_tio2 = pc('Input_files/spiro_mapi_tio2.csv');

par_tio2 = pc('Input_files/p_n_heterojunction.csv');
soleq_tio2 = equilibrate(par_tio2);
sol_CV_100mVs_tio2_ill = doCV(soleq_tio2.el, 1, 0, 1.0, 0, 100e-3, 1, 241);
compare_rec_flux(sol_CV_100mVs_tio2_dark, 0.02, 1e6, 1);