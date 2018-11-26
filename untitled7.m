% soleq_toy_ntyp0p3eV = equilibrate(par.toy.ntyp0p3eV);
% soleq_toy_ptyp0p3eV = equilibrate(par.toy.ptyp0p3eV);
% 
% sdpsol_toy_ntype0p2eV = sdp(soleq_toy_ntypewf_nobulksrh.i_sr, logspace(-6, 2, 9), 0.7, 10e-6, 1e-3, 1);
% sdpsol_toy_ntype0p3eV = sdp(soleq_toy_ntyp0p3eV.i_sr, logspace(-6, 2, 9), 0.7, 10e-6, 1e-3, 1);
% sdpsol_toy_ptype0p2eV = sdp(soleq_toy_ptypewf_nobulksrh.i_sr, logspace(-6, 2, 9), 0.7, 10e-6, 1e-3, 1);
% sdpsol_toy_ptype0p3eV = sdp(soleq_toy_ptyp0p3eV.i_sr, logspace(-6, 2, 9), 0.7, 10e-6, 1e-3, 1);


sdpanal(sdpsol_toy_ntype0p3eV, 6e-6);
figure(1)
hold on
sdpanal(sdpsol_toy_ntype, 6e-6);

sdpana(sdpsol_toy_ptype, 6e-6);

sdpana(sdpsol_toy_ptype0p3eV, 6e-6);
