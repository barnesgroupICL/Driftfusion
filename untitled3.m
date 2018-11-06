% soleq_tio2 = equilibrate(par.tio2.base);
% soleq_pedot = equilibrate(par.pedot.base);
% 
% JV_tio2.JV0p5mVs = doJV(soleq_tio2.i_sr, 0.5e-3, 100, 1, 1, 0, 1.4, 2);
% JV_tio2.JV5mVs = doJV(soleq_tio2.i_sr, 5e-3, 100, 1, 1, 0, 1.4, 2);
% JV_tio2.JV50mVs = doJV(soleq_tio2.i_sr, 50e-3, 100, 1, 1, 0, 1.4, 2);
% JV_tio2.JV500mVs = doJV(soleq_tio2.i_sr, 500e-3, 100, 1, 1, 0, 1.4, 2);
% JV_tio2.JV5000mVs = doJV(soleq_tio2.i_sr, 5000e-3, 100, 1, 1, 0, 1.4, 2);
% 
% JV_pedot.JV0p5mVs = doJV(soleq_pedot.i_sr, 0.5e-3, 100, 1, 1, 0, 1.4, 2);
% JV_pedot.JV5mVs = doJV(soleq_pedot.i_sr, 5e-3, 100, 1, 1, 0, 1.4, 2);
% JV_pedot.JV50mVs = doJV(soleq_pedot.i_sr, 50e-3, 100, 1, 1, 0, 1.4, 2);
% JV_pedot.JV500mVs = doJV(soleq_pedot.i_sr, 500e-3, 100, 1, 1, 0, 1.4, 2);
% JV_pedot.JV5000mVs = doJV(soleq_pedot.i_sr, 5000e-3, 100, 1, 1, 0, 1.4, 2);
% 
% soleq_pedot_lowVbi = equilibrate(par.pedot.lowVbi);
% JV_pedot_lowVbi.JV50mVs = doJV(soleq_pedot_lowVbi.i_sr, 50e-3, 100, 1, 1, 0, 1.4, 2);
% JV_pedot_lowVbi.JV5000mVs = doJV(soleq_pedot_lowVbi.i_sr, 5000e-3, 100, 1, 1, 0, 1.4, 2);


JV_tio2.JV500mVs = doJV(soleq_tio2.i_sr, 500e-3, 100, 1, 1, 0, 1.4, 2);
JV_pedot.JV500mVs = doJV(soleq_pedot.i_sr, 500e-3, 100, 1, 1, 0, 1.4, 2);