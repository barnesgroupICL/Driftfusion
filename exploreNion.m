function parexsol = exploreNion(par_base)

% par is the base parameter set
tic
%tau_ptype = logspace(-13, -6, 8);
N_ion = logspace(16, 20, 5);
E0_ptype = (-5.0:0.1:-4.4);
parexsol.Voc_f = zeros(length(N_ion), length(E0_ptype));
parexsol.Voc_r = zeros(length(N_ion), length(E0_ptype));

j = 1;

for i = 1:length(N_ion)
    par = par_base;
    par.Ana = 0;
    par.NI = N_ion(i);
    
    for j = 1:length(E0_ptype)
        
        runN = (i-1)*length(E0_ptype) + j;
        disp(['Run no. ', num2str(runN), ', taun = ', num2str(N_ion(i)), ', E0 = ', num2str(E0_ptype(j))]);
    
        par.E0(1) = E0_ptype(j);
        soleq = equilibrate(par);
        JV = doJV(soleq.i_sr, 50e-3, 100, 1, 1e-10, 0, 1.5, 2);

        parexsol.Voc_f(i, j) = JV.Voc_f;
        parexsol.Voc_r(i, j) = JV.Voc_r;
        
    end
    
end

parexsol.N_ion = N_ion;
parexsol.E0_ptype = E0_ptype;

figure(3000)
surf(E0_ptype-par.IP(1), N_ion, parexsol.Voc_f)
ylabel('N_ion [cm-3]')
xlabel('Fermi level offset [eV]')
zlabel('Voc F scan [V]')

figure(3001)
surf(E0_ptype-par.IP(1), N_ion, parexsol.Voc_r)
ylabel('N_ion [cm-3]')
xlabel('Fermi level offset [eV]')
zlabel('Voc R scan [V]')

toc
end
        
        
        
