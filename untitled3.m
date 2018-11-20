%soleq_tio2 = equilibrate(par.tio2.base);
%soleq_tio2_Vbi1p1V = equilibrate(par.tio2.Vbi1p1V);
%soleq_pedot = equilibrate(par.pedot.base);
% soleq_toy_lowepp = equilibrate(par.toy.lowepp);
% soleq_toy_hiepp = equilibrate(par.toy.hiepp);
% 
% JV_tio2.JV0p5mVs = doJV(soleq_tio2.i_sr, 0.5e-3, 100, 1, 1, 0, 1.4, 2);
% JV_tio2.JV5mVs = doJV(soleq_tio2.i_sr, 5e-3, 100, 1, 1, 0, 1.4, 2);
% JV_tio2.JV50mVs = doJV(soleq_tio2.i_sr, 50e-3, 100, 1, 1, 0, 1.4, 2);
% JV_tio2.JV500mVs = doJV(soleq_tio2.i_sr, 500e-3, 100, 1, 1, 0, 1.4, 2);
% JV_tio2.JV5000mVs = doJV(soleq_tio2.i_sr, 5000e-3, 100, 1, 1, 0, 1.4, 2);
% 
% JV_tio2_Vbi1p1V.JV0p5mVs = doJV(soleq_tio2.i_sr, 0.5e-3, 100, 1, 1, 0, 1.4, 2);
% JV_tio2_Vbi1p1V.JV5mVs = doJV(soleq_tio2.i_sr, 5e-3, 100, 1, 1, 0, 1.4, 2);
% JV_tio2_Vbi1p1V.JV50mVs = doJV(soleq_tio2.i_sr, 50e-3, 100, 1, 1, 0, 1.4, 2);
% JV_tio2_Vbi1p1V.JV500mVs = doJV(soleq_tio2.i_sr, 500e-3, 100, 1, 1, 0, 1.4, 2);
% JV_tio2_Vbi1p1V.JV5000mVs = doJV(soleq_tio2.i_sr, 5000e-3, 100, 1, 1, 0, 1.4, 2);
% 
% JV_pedot.JV0p5mVs = doJV(soleq_pedot.i_sr, 0.5e-3, 100, 1, 1, 0, 1.4, 2);
% JV_pedot.JV5mVs = doJV(soleq_pedot.i_sr, 5e-3, 100, 1, 1, 0, 1.4, 2);
% JV_pedot.JV50mVs = doJV(soleq_pedot.i_sr, 50e-3, 100, 1, 1, 0, 1.4, 2);
% JV_pedot.JV500mVs = doJV(soleq_pedot.i_sr, 500e-3, 100, 1, 1, 0, 1.4, 2);
% JV_pedot.JV5000mVs = doJV(soleq_pedot.i_sr, 5000e-3, 100, 1, 1, 0, 1.4, 2);

% JV_toy_lowepp.JV0p5mVs = doJV(soleq_toy_lowepp.i_sr, 0.5e-3, 100, 1, 1, 0, 1.4, 2);
% JV_toy_lowepp.JV5mVs = doJV(soleq_toy_lowepp.i_sr, 5e-3, 100, 1, 1, 0, 1.4, 2);
% JV_toy_lowepp.JV50mVs = doJV(soleq_toy_lowepp.i_sr, 50e-3, 100, 1, 1, 0, 1.4, 2);
% JV_toy_lowepp.JV500mVs = doJV(soleq_toy_lowepp.i_sr, 500e-3, 100, 1, 1, 0, 1.4, 2);
% JV_toy_lowepp.JV5000mVs = doJV(soleq_toy_lowepp.i_sr, 5000e-3, 100, 1, 1, 0, 1.4, 2);

% JV_toy_hiepp.JV0p5mVs = doJV(soleq_toy_hiepp.i_sr, 0.5e-3, 100, 1, 1, 0, 1.4, 2);
% JV_toy_hiepp.JV5mVs = doJV(soleq_toy_hiepp.i_sr, 5e-3, 100, 1, 1, 0, 1.4, 2);
% JV_toy_hiepp.JV50mVs = doJV(soleq_toy_hiepp.i_sr, 50e-3, 100, 1, 1, 0, 1.4, 2);
% JV_toy_hiepp.JV500mVs = doJV(soleq_toy_hiepp.i_sr, 500e-3, 100, 1, 1, 0, 1.4, 2);
% JV_toy_hiepp.JV5000mVs = doJV(soleq_toy_hiepp.i_sr, 5000e-3, 100, 1, 1, 0, 1.4, 2);

% jumpto0p8V_tio2 = jumptoV(soleq_tio2.i_sr, 0.8);
% jumpto0p8V_tio2_Vbi1p1V = jumptoV(soleq_tio2_Vbi1p1V.i_sr, 0.8);
% jumpto0p8V_pedot = jumptoV(soleq_pedot.i_sr, 0.8);

%jumpto0p8V_toy_lowepp = jumptoV(soleq_toy_lowepp.i_sr, 0.8);
%jumpto0p8V_toy_hiepp = jumptoV(soleq_toy_hiepp.i_sr, 0.8);
% 
% Fjump_tio2 = moviemake(jumpto0p8V_tio2);
% Fjump_tio2_Vbi1p1V = moviemake(jumpto0p8V_tio2_Vbi1p1V);
% Fjump_pedot = moviemake(jumpto0p8V_pedot);
% 
% Fjump_toy_lowepp = moviemake(jumpto0p8V_toy_lowepp);
% Fjump_toy_hiepp = moviemake(jumpto0p8V_toy_hiepp);


JV_tio2_fixedion = doJV(soleq_tio2.i_sr, 5e-3, 100, 1, 0, 0, 1.0, 2);
JV_pedot_fixedion = doJV(soleq_pedot.i_sr, 5e-3, 100, 1, 0, 0, 1.0, 2);
