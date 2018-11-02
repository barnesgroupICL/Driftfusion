soleq_tio2 = equilibrate(par.tio2.base);
soleq_pedot = equilibrate(par.pedot.base);
JV_tio2.JV50mVs = doJV(soleq_tio2.i_sr, 50e-3, 100, 1, 1, 0, 1.4, 2);
JV_tio2.JV5000mVs = doJV(soleq_tio2.i_sr, 5000e-3, 100, 1, 1, 0, 1.4, 2);
JV_pedot.JV50mVs = doJV(soleq_pedot.i_sr, 50e-3, 100, 1, 1, 0, 1.4, 2);
JV_pedot.JV5000mVs = doJV(soleq_pedot.i_sr, 5000e-3, 100, 1, 1, 0, 1.4, 2);