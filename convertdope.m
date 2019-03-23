function Ef0 = convertdope(Ndope,Nband,type,par,Eband)

% type is a string - eith 'NA' or 'ND'
kT = par.kB*par.T;
    
if strcmp('ND', type) == 1

    Ef0 = Eband+kT*log(Ndope/Nband);
    
elseif strcmp('NA', type) == 1
    
    Ef0 = Eband-kT*log(Ndope/Nband);
end

end