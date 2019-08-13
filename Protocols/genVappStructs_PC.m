function VappStruct = genVappStructs_PC(solini, Vapp_arr)

% solstruct should be a stablised solution

par = solini.par;    
pini = par;           % Save initial params


Vapp_arr = [solini.par.Vapp, Vapp_arr];    % include intial potential in array

% General settings
par.Ana = 1;
par.figson = 1;
par.pulseon = 0;
par.int1 = 0;

for i = 1:length(Vapp_arr)-1
    
    disp([mfilename ' - applied voltage ' num2str(Vapp_arr(i+1))])
    name = matlab.lang.makeValidName([inputname(1) '_Vapp_' num2str(Vapp_arr(i+1))]);

    par.mobseti = 1e4;
    par.JV = 1;
    par.Vapp = solini.par.Vapp;
    par.Vstart = Vapp_arr(1);
    par.Vend = Vapp_arr(i+1);
    par.tmesh_type = 1;
    par.tmax = 1e-2;
    par.t0 = 0;%par.tmax/1e6;
    
    sol = df(solini, par);
    
    % store new applied potential
    sol.par.Vapp = par.Vend;
    par.Vapp = par.Vend;
    %Switch on ion mobility and stabilise
    %par.mui = 1e-6;
    par.JV = 0;
    par.tmesh_type = 2;
    par.tpoints = 100;
    par.tmax = 1;
    par.t0 = par.tmax/1e8;
    
    sol = df(sol, par);
   
    verifyStabilization(sol.sol, sol.t, 0.7);
    
    sol.par.mobseti = 1;
    
    % if there's only one solution then duplicate sol structure
    if length(Vapp_arr)-1 == 1
        
        VappStruct = sol;
    
    % if there's multiple solutions, store in a master struct
    else
        
        VappStruct{1, i} = sol;
        VappStruct{2, i} = name;
        
    end
    
end

end

