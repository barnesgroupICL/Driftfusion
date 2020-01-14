function VappSol = genVappStructs(solini, Vapp_arr, mobseti)
% Obtain steady-state solutions for applied voltages defined by VAPP_ARR
% solstruct should be a stablised solution
% If MOBSETI = 1 ions are accelerated to obtain steady state
% P Calado 2019, Imperial College London
% Adapted from original code by I. Gelmetti

% Store parameters
par = solini.par;

Vapp_arr = [solini.par.Vapp, Vapp_arr];    % include intial potential in array in case

% Store initial solution
sol = solini;

for i = 1:length(Vapp_arr)-1
    
    disp([mfilename ' - applied voltage ' num2str(Vapp_arr(i+1))])
    name = matlab.lang.makeValidName([inputname(1) '_Vapp_' num2str(Vapp_arr(i+1))]);
    
    if mobseti
        % Take ratio of electron and ion mobilities in the active layer
        rat_anion = par.mue(par.active_layer)/par.muani(par.active_layer);
        rat_cation = par.mue(par.active_layer)/par.mucat(par.active_layer);
        
        % If the ratio is infinity (ion mobility set to zero) then set the ratio to
        % zero instead
        if isnan(rat_anion) || isinf(rat_anion)
            rat_anion = 0;
        end
        
        if isnan(rat_cation) || isinf(rat_cation)
            rat_cation = 0;
        end
        par.mobseti = 1;           % Ions are accelerated to reach equilibrium
        par.K_anion = rat_anion;
        par.K_cation = rat_cation;
    end
    
    par.tmesh_type = 1;
    par.t0 = 0;
    par.tmax = 1e-2;
    par.V_fun_type = 'sweep';
    par.V_fun_arg(1) = Vapp_arr(i);
    par.V_fun_arg(2) = Vapp_arr(i+1);
    par.V_fun_arg(3) = par.tmax;
    
    sol = df(sol, par);
    
    par.V_fun_type = 'constant';
    par.V_fun_arg(1) = Vapp_arr(i+1);

    sol = df(sol, par);
        
    % Check that the solution is steady-state
    all_stable = verifyStabilization(sol.u, sol.t, 0.7);
    j = 0;
    while any(all_stable) == 0
        disp(['increasing equilibration time, tmax = ', num2str(par.tmax*10^j)]);
        par.tmax = par.tmax*10;
        par.t0 = par.tmax/1e6;
        
        sol = df(sol, par);
        
        all_stable = verifyStabilization(sol.u, sol.t, 0.7);
    end
    
    % reset ion coeffs before storing
    if mobseti
        sol.par.K_anion = 1;
        sol.par.K_cation = 1;
    end
    % if there's only one solution then duplicate sol structure
    if length(Vapp_arr)-1 == 1
        VappSol = sol;
        % if there's multiple solutions, store in a master struct
    else
        VappSol{1, i} = sol;
        VappSol{2, i} = name;
    end
end

end

