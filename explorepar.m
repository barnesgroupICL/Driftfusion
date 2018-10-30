function parexsol = explorepar(par_base)

% par_base is the base parameter set
tic
disp('Starting parameter exploration');
% disp(['Parameter 1: ', parname1]);
% disp(['Parameter 2: ', parname2]);
tau_ptype = logspace(-13, -6, 8);
E0_ptype = (-5.0:0.1:-4.4);
j = 1;

pararr1 = tau_ptype;
pararr2 = E0_ptype;

parfor i = 1:length(pararr1)
 
    par = par_base;
    par.Ana = 0;
    par.taun(1) = tau_ptype(i)
    par.taup(1) = tau_ptype(i);
    
    Voc_f = zeros(1, length(pararr2));
    Voc_r = zeros(1, length(pararr2));    
    Jsc_f = zeros(1, length(pararr2));
    Jsc_r = zeros(1, length(pararr2));
    mpp_f = zeros(1, length(pararr2));
    mpp_r = zeros(1, length(pararr2));
    FF_f = zeros(1, length(pararr2));
    FF_r = zeros(1, length(pararr2));
    Voc_stable = zeros(1, length(pararr2));
    PLint = zeros(1, length(pararr2));

    for j = 1:length(pararr2)
        
        runN = (i-1)*length(pararr2) + j;
        disp(['Run no. ', num2str(runN), ', taun = ', num2str(pararr1(i)), ', E0 = ', num2str(pararr2(j))]);
        par.E0(1) = E0_ptype(j);
        
        soleq = equilibrate(par);
        JV = doJV(soleq.i_sr, 50e-3, 100, 1, 1e-10, 0, 1.5, 2);
        
        Voc_f(j) = JV.stats.Voc_f;
        Voc_r(j) = JV.stats.Voc_r;
        Jsc_f(j) = JV.stats.Jsc_f;
        Jsc_r(j) = JV.stats.Jsc_r;
        mpp_f(j) = JV.stats.mpp_f;
        mpp_r(j) = JV.stats.mpp_r;
        FF_f(j) = JV.stats.FF_f;
        FF_r(j) = JV.stats.FF_r;
        
        % For PL
        [sol_Voc, Voc] = findVoc(soleq.i_sr, 1e-6, Voc_f(j), (Voc_f(j)+0.1))
        Voc_stable(j) = Voc;
        PLint(j) = sol_Voc.PLint(end);
            
    end
    
    A(i,:) = Voc_f;
    B(i,:) = Voc_r;
    C(i,:) = Jsc_f;
    D(i,:) = Jsc_r;
    E(i,:) = mpp_f;
    F(i,:) = mpp_r;
    G(i,:) = FF_f;
    H(i,:) = FF_r;
    J(i,:) = Voc_stable;
    K(i,:) = PLint;
    
end

parexsol.stats.Voc_f = A;
parexsol.stats.Voc_r = B;
parexsol.stats.Jsc_f = C;
parexsol.stats.Jsc_r = D;
parexsol.stats.mpp_f = E;
parexsol.stats.mpp_r = F;
parexsol.stats.FF_f = G;
parexsol.stats.FF_r = H;
parexsol.stats.Voc_stable = J;
parexsol.stats.PLint = K;
parexsol.tau_ptype = tau_ptype;
parexsol.E0_ptype = E0_ptype;

% 
% figure(3000)
% surf(pararr2-par.IP(1), pararr1, parexsol.Voc_f)
% ylabel('pararr1 [cm-3]')
% xlabel('Fermi level offset [eV]')
% zlabel('Voc F scan [V]')
% 
% figure(3001)
% surf(pararr2-par.IP(1), pararr1, parexsol.Voc_r)
% ylabel('pararr1 [cm-3]')
% xlabel('Fermi level offset [eV]')
% zlabel('Voc R scan [V]')

toc
end
        
        
        
