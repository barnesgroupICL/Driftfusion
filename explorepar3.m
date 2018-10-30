function parexsol = explorepar3(par_base, parnames, parvalues)

% par_base is the base parameter set
% parnames is a cell array with the parameter names in - check these
% carefully to avoid heartache later
% parvalues is matrix with the parameter value ranges e.g.

% Current version for 2 parameters only- will be extended in the future
tic
disp('Starting parameter exploration');
% disp(['Parameter 1: ', parnames(1)]);
% disp(['Parameter 2: ', parnames(2)]);
for k = 1:length(parnames)

disp(['Parameter ', num2str(k),':' parnames(k)]);


end


parval = cell2mat(parvalues(k));
parval2 = cell2mat(parvalues(2));
name = char(parnames(k));
str2 = char(parnames(2));

j = 1;

parfor i = 1:length(parval1)
 
    par = par_base;
    par.Ana = 0;
    par = exploreparhelper(par, str1, parval1(i));
    par.taup(1) = par.taun(1)
    
    Voc_f = zeros(1, length(parval2));
    Voc_r = zeros(1, length(parval2));    
    Jsc_f = zeros(1, length(parval2));
    Jsc_r = zeros(1, length(parval2));
    mpp_f = zeros(1, length(parval2));
    mpp_r = zeros(1, length(parval2));
    FF_f = zeros(1, length(parval2));
    FF_r = zeros(1, length(parval2));
    Voc_stable = zeros(1, length(parval2));
    PLint = zeros(1, length(parval2));

    for j = 1:length(parval2)
        
        runN = (i-1)*length(parval2) + j;
        disp(['Run no. ', num2str(runN), ', taun = ', num2str(parval1(i)), ', E0 = ', num2str(parval2(j))]);
  
        par = exploreparhelper(par, str2, parval2(j));
        
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
parexsol.parnames = parnames;
parexsol.parvalues = parvalues;
parexsol.parval1 = parval1;
parexsol.parval2 = parval2;

% 
% figure(3000)
% surf(parval2-par.IP(1), parval1, parexsol.Voc_f)
% ylabel('parval1 [cm-3]')
% xlabel('Fermi level offset [eV]')
% zlabel('Voc F scan [V]')
% 
% figure(3001)
% surf(parval2-par.IP(1), parval1, parexsol.Voc_r)
% ylabel('parval1 [cm-3]')
% xlabel('Fermi level offset [eV]')
% zlabel('Voc R scan [V]')

end

toc
end
        
        
        
