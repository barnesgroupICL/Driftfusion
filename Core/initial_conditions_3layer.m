function [uic] = initial_conditions_3layer(par, figson)
% A function to obtain initial conditions for a 3 layer device with mobile
% ionic charge
xnm = par.xx*1e7;
dev = par.dev;
pcum1 = par.pcum0+1;

% obtain electrostatic potential from depletion approximation based on
% doping densites and ionic densities
[rhomat, Fmat, Phimat] = depletion_approx_modelX_numeric(par, 1, 0, 0.1, 2, [0],0);

V0 = Phimat(1,:);

% Assume QFL are constant at equilibrium to obtain carrier densities
n = dev.Nc.*exp((par.PhiA - (dev.EA - V0))/(par.kB*par.T));
p = dev.Nv.*exp(((dev.IP - V0) - par.PhiA)/(par.kB*par.T));

% To obtain ionic charge densities
% Active layer midpoint
p_midact = round((pcum1(3)+pcum1(4))/2);

c = ones(1, length(par.xx)).*dev.Ncat;
a = ones(1, length(par.xx)).*dev.Nani;

% Point limits for the active layer
p_active_l = pcum1(par.active_layer(1)-1)+1;
p_active_r = pcum1(par.active_layer(end)+2);

V0_active = ones(1, length(par.xx))*V0(p_midact);
V0_active(p_active_l:p_active_r) = V0(p_active_l:p_active_r);

%% Get appropriate ionic densities of states (Beta) that ensure the total integrated charge is the same as that
% of the uniform integrated charge
Ncat_int = trapz(par.xx(p_active_l:p_active_r), dev.Ncat(p_active_l:p_active_r));
Nani_int = trapz(par.xx(p_active_l:p_active_r), dev.Nani(p_active_l:p_active_r));

Beta_c = Ncat_int/(trapz(par.xx(p_active_l:p_active_r), ...
    exp((V0_active(p_midact) - V0_active(p_active_l:p_active_r))/(par.kB*par.T))));

Beat_a = Nani_int/(trapz(par.xx(p_active_l:p_active_r), ...
    exp((V0_active(p_active_l:p_active_r) - V0_active(p_midact))/(par.kB*par.T))));

% c(p_active_l:p_active_r) = Beta_c.*exp((V0_active(p_midact) - V0_active(p_active_l:p_active_r))/(par.kB*par.T));
% a(p_active_l:p_active_r) = Beat_a.*exp((V0_active(p_active_l:p_active_r) - V0_active(p_midact))/(par.kB*par.T));

c(p_active_l:p_active_r) = dev.Ncat(p_active_l:p_active_r).*exp((V0_active(p_midact) - V0_active(p_active_l:p_active_r))/(par.kB*par.T));
a(p_active_l:p_active_r) = dev.Nani(p_active_l:p_active_r).*exp((V0_active(p_active_l:p_active_r) - V0_active(p_midact))/(par.kB*par.T));

c = dev.Ncat.*exp((V0_active(p_midact) - V0_active)/(par.kB*par.T));
a = dev.Nani.*exp((V0_active - V0_active(p_midact))/(par.kB*par.T));

% Normalise to uniform integrated density in the active layer
Ncat_int = trapz(par.xx(p_active_l:p_active_r), dev.Ncat(p_active_l:p_active_r));
Nani_int = trapz(par.xx(p_active_l:p_active_r), dev.Nani(p_active_l:p_active_r));
c_int = trapz(par.xx(p_active_l:p_active_r), c(p_active_l:p_active_r));
a_int = trapz(par.xx(p_active_l:p_active_r), a(p_active_l:p_active_r));

% % We only want the ionic charge in the active layer to be normalised
% 
% 
% % Check
% c_norm_int = trapz(par.xx(p_active_l:p_active_r), c_norm(p_active_l:p_active_r));
% a_norm_int = trapz(par.xx(p_active_l:p_active_r), a_norm(p_active_l:p_active_r));

% Subtract electron and hole densities from the space charge density to
% obtain ionic charge- will be step function
rho = rhomat(1,:)/par.e;
rhoel = -n + p + dev.ND - dev.NA;
rhoc = c - a;

uic(:,1) = n;
uic(:,2) = p;
uic(:,3) = c;
uic(:,4) = V0;

if par.N_ionic_species == 2
    uic(:,5) = a;
end

if figson
    figure(501)
    semilogy(xnm, n, xnm, p)
    xlabel('Position [nm]')
    ylabel('Electronic carrier density [cm-3]')
    legend('n','p')
    
    figure(502)
    plot(xnm, rho, xnm, rhoc, xnm, rhoel)
    xlabel('Position [nm]')
    ylabel('Space charge density [cm-3]')
    
    figure(503)
    plot(xnm, a, xnm, c);%, c_norm, xnm, a_norm)
    xlabel('Position [nm]')
    ylabel('Ionic carrier density [cm-3]')
    legend('a','c');%,'c_{norm}','a_{norm}')
    
    figure(504)
    plot(xnm, -V0)
    xlabel('Position [nm]')
    ylabel('Potential [V]')
end
end

