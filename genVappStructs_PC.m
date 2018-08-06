function VappStruct = genVappStructs_PC(solini, Vapp_arr)

% solstruct should be a stablised solution

p = solini.p;    
pini = p;           % Save initial params

Vapp_arr = [solini.p.Vapp, Vapp_arr];    % include intial potential in array

% General settings
p.Ana = 1;
p.figson = 1;
p.pulseon = 0;
p.Int = 0;

for i = 1:length(Vapp_arr)-1
    
    disp([mfilename ' - applied voltage ' num2str(Vapp_arr(i+1))])
    name = matlab.lang.makeValidName([inputname(1) '_Vapp_' num2str(Vapp_arr(i+1))]);
    % run initial sweep with ion mobility off
    p.mui = 1e-6;
    p.JV = 1;
    p.Vapp = solini.p.Vapp;
    p.Vstart = Vapp_arr(1);
    p.Vend = Vapp_arr(i+1);
    p.tmesh_type = 1;
    p.tmax = 1e-2;
    p.t0 = 0;%p.tmax/1e6;
    
    sol = pindrift(solini, p);
    
    % store new applied potential
    sol.p.Vapp = p.Vend;
    p.Vapp = p.Vend;
    %Switch on ion mobility and stabilise
    %p.mui = 1e-6;
    p.JV = 0;
    p.tmesh_type = 2;
    p.tpoints = 200;
    p.tmax = 1;
    p.t0 = p.tmax/1e8;
    
    sol = pindrift(sol, p);
   
    verifyStabilization(sol.sol, sol.t, 0.7);
       
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

