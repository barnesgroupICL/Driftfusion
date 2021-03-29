% File to plot the array of parameters builted in  build_device.m
% run the code before to the end of build_device.m

% Constant properties
figure; plot(xmesh*1e7,dev.mucat); ylabel('mucat [cm^2/Vs]'); xlabel('x [nm]');
figure; plot(xmesh*1e7,dev.muani); ylabel('muani [cm^2/Vs]'); xlabel('x [nm]');
figure; plot(xmesh*1e7,dev.epp); ylabel('epp '); xlabel('x [nm]');
figure; plot(xmesh*1e7,dev.NTSD); ylabel('NTSD [cm^{-3}]'); xlabel('x [nm]');
figure; plot(xmesh*1e7,dev.NTSA); ylabel('NTSA [cm^{-3}]'); xlabel('x [nm]');
 
% Linearly graded properties
figure; plot(xmesh*1e7,dev.EA,'b',xmesh*1e7,dev.IP,'r',xmesh*1e7,dev.E0,'k');
ylabel('EA, IP, E0 [eV]'); xlabel('x [nm]'); legend('EA','IP','E0');

% Linearly graded properties
figure; plot(xmesh*1e7,dev.NA,'b',xmesh*1e7,dev.IP,'r',xmesh*1e7,dev.E0,'k');
ylabel('EA, IP, E0 [eV]'); xlabel('x [nm]'); legend('EA','IP','E0');

% Logarithmically graded properties
figure; plot(xmesh*1e7,dev.Nc,'r-',xmesh*1e7,dev.Nv,'b-'); 
ylabel('Nc, Nv [cm^{-3}]'); xlabel('x [nm]');legend('Nc','Nv');
figure; plot(xmesh*1e7,dev.n0,'r-',xmesh*1e7,dev.p0,'b-'); 
ylabel('n0, p0 [cm^{-3}]'); xlabel('x [nm]');legend('n0','p0');
figure; plot(xmesh*1e7,dev.Nani,'r-',xmesh*1e7,dev.Ncat,'b-'); 
ylabel('Nani, Ncat [cm^{-3}]'); xlabel('x [nm]');legend('Nani','Ncat');
figure; plot(xmesh*1e7,dev.ni); ylabel('ni (log_graded) [cm^{-3}]'); xlabel('x [nm]');
figure; plot(xmesh*1e7,dev.DOSani,'r-',xmesh*1e7,dev.DOScat,'b-'); 
ylabel('par.amax (DOSani), par.cmax (log_graded) [cm^{-3}]'); xlabel('x [nm]');legend('Nc','Nv');

% Properties that are zeroed in the interfaces
figure; plot(xmesh*1e7,dev.NA,'r-',xmesh*1e7,dev.NA,'b-'); 
ylabel('NA, ND [cm^{-3}]'); xlabel('x [nm]');legend('NA','ND');
figure; plot(xmesh*1e7,dev.g0,'b',xmesh*1e7,dev.B,'r');
ylabel('g0, B [cm^{-3}/s]'); xlabel('x [nm]'); legend('g0','B');

%  Gradient properties
figure; plot(xmesh*1e7,dev.gradEA,'r--',xmesh*1e7,dev.gradIP,'b-'); 
ylabel('gradEA, gradIP (lin graded) [cm^{-3}/cm]'); xlabel('x [nm]');legend('gradEA','gradIP');
figure; plot(xmesh*1e7,dev.gradNc,'r--',xmesh*1e7,dev.gradNv,'b-'); 
ylabel('gradNc, gradNv (log graded) [cm^{-3}/cm]'); xlabel('x [nm]');legend('gradNc','gradNv');

% Surface recombination velocity equivalence schemes
figure; plot(xmesh*1e7,dev.mue,'b--',xmesh*1e7,dev.muh,'r--');
ylabel('mue, muh [cm^{2}/Vs]'); xlabel('x [nm]');legend('mue','muh');
figure; plot(xmesh*1e7,dev.taun,'b--',xmesh*1e7,dev.taup,'r--');
ylabel('taun, taup [s]'); xlabel('x [nm]');legend('taun','taup');

% Trap Density
figure; plot(xmesh*1e7,dev.nt,'b--',xmesh*1e7,dev.pt,'r--');
ylabel('nt, pt [cm^{-3}]'); xlabel('x [nm]');legend('nt','pt');

figure; plot(xmesh*1e7,dev.ni_srh);
ylabel('ni_{srh} [cm^{-3}]'); xlabel('x [nm]');

figure; plot(xmesh*1e7,dev.int_switch);
ylabel('int_{switch} [cm^{-3}]'); xlabel('x [nm]');


