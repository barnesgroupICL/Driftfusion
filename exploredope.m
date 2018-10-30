function parexsol = exploredope(par_base)

% par is the base parameter set
tic
tau_ptype = logspace(-13, -6, 8);
E0_ptype = (-5.0:0.1:-4.4);
parexsol.Voc_f = zeros(length(tau_ptype), length(E0_ptype));
parexsol.Voc_r = zeros(length(tau_ptype), length(E0_ptype));

j = 1;

for i = 1:length(tau_ptype)
    par = par_base;
    par.Ana = 0;
    par.taun(1) = tau_ptype(i);
    par.taup(1) = tau_ptype(i);
    
    for j = 1:length(E0_ptype)
        
        runN = (i-1)*length(E0_ptype) + j;
        disp(['Run no. ', num2str(runN), ', taun = ', num2str(tau_ptype(i)), ', E0 = ', num2str(E0_ptype(j))]);
    
        par.E0(1) = E0_ptype(j);
        soleq = equilibrate(par);
        JV = doJV(soleq.i_sr, 50e-3, 100, 1, 1e-10, 0, 1.5, 2);

        parexsol.Voc_f(i, j) = JV.Voc_f;
        parexsol.Voc_r(i, j) = JV.Voc_r;
        
    end
    
end

parexsol.tau_ptype = tau_ptype;
parexsol.E0_ptype = E0_ptype;

figure(3000)
surf(E0_ptype-par.IP(1), tau_ptype, parexsol.Voc_f)
ylabel('tau_ptype [s]')
xlabel('Fermi level offset [eV]')
zlabel('Voc F scan [V]')

figure(3001)
surf(E0_ptype-par.IP(1), tau_ptype, parexsol.Voc_r)
ylabel('tau_ptype [s]')
xlabel('Fermi level offset [eV]')
zlabel('Voc R scan [V]')

toc
end
        
        
        
