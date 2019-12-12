function [tau, Cap, rhomat, Fmat, Phimat, Q, J, tout] = depletion_approx_modelX1_fungen(par, tarr, figson)
% Philip Calado 2019
% Department of Physics, Imperial College London
% Solves the depletion approximation for a device based on the parameters
% given in PAR to obtain the charge density RHO, electric field F, and electrostatic potential Phi
% PAR is a parameters object created by PC in DRIFTFUSION. It can be
% created using the command:
% par = pc('input_files/3 layer depletion test.csv'), where the file path must reference as .csv file
% This code asumes a three layer model, with each layer separated by a
% junction. Consequently the indices 1,3, & 5 reference the HTL, perovskite, and
% HTL respectively
% s = the scan rate
% VSTART = the scan starting applied voltage.
% VEND = the scan final applied voltage.
% TPOINTS = the number of desired time points
% TARR = An array containing the times that the user wishes to plot

%% MODEL
% This model asumes that the space charge in the contact regions is electronic,
% while that in the perovskite layer is ionic. Step functions are used for
% the charge in the depletion region leading to quadratic changes in the
% potential as with the original depletion approximation.
% The labels 1,2, & 3 refer to the different layers for example V2 is the
% potential drop in the perovskite at BOTH INTERFACES. This saves
% repetition and maintains consistency with the labels for materials
% properties e.g. epp2
% The model separates the charge into two components: space charge due to
% the junction capacitance and a second contribution due to the applied potential.
% This results in two components to the depletion width within the
% transport layers:
% In the HTL, w1 is the width due to the junction capacitance and w1app is the
% componenet due to the applied potential. Similarly in the ETL the
% componenets are w3 and wapp3 (in keeping with labelling standard)
% The rate of change of Q, dQdt is asumed to be determined by the drift
% current of ions acros the perovskite and this is solved for numerically
% using ODE15S
% Currently V(t) is set to be linear but any arbitrary function can be
% input
%% Version history
% v6 Qapp is asumed to be on the electrodes - allows field at the
% boundaires. Symmetirc device only.
% v7 zero field boundary conditions- all charge is contained within the
% system. Symmetric device only.
% v9 adds the ability to model asymmetric contacts by introducing 2
% variables for the additional depletion width due to the applied voltage,
% wapp1 and wapp3- this is solved for by using expresions F(0) = F(d) = 0 and
% Phi(0) = 0 & Phi(d) = V
% This is achieved in solve_wapp1_wapp3_model9
% vX1 adds function generator compatibility - takes parameters from par

%% Spatial mesh
x = par.xx;
xnm = x*1e7;
d = x(end);
pmid = par.pcum0(par.active_layer) +...
    round((par.pcum0(par.active_layer+1)-par.pcum0(par.active_layer))/2);


%% Constants
epp0 = 8.8541878128e-14;
q = par.e;

% Doping densities in each layer
N1 = par.NA(1);
N2 = par.Ncat(par.active_layer);
N3 = par.ND(end);

epp1 = par.epp(1)*epp0;
epp2 = par.epp(par.active_layer)*epp0;
epp3 = par.epp(end)*epp0;

Vbi = par.Vbi;

% Q0 is the initial charge in the space charge layer widths excluding Qapp1.
% This is solved for analytically at steady state by setting dQdt = 0;
Q0 = (2^(1/2)*N1*N2*N3*epp1*epp2*epp3*q*((Vbi*(N1*N2*epp1*epp2 + 2*N1*N3*epp1*epp3 + N2*N3*epp2*epp3))...
    /(N1*N2*N3*epp1*epp2*epp3*q))^(1/2))/(N1*N2*epp1*epp2 + 2*N1*N3*epp1*epp3 + N2*N3*epp2*epp3);

V1 = (Q0^2)/(2*q*N1*epp1);
V2 = (Q0^2)/(2*q*N2*epp1);
V3 = (Q0^2)/(2*q*N3*epp1);

% Equilibrium depletion widths
w1 = Q0/(q*N1);
w2 = Q0/(q*N2);
w3 = Q0/(q*N3);

d0 = 0;
d1 = par.dcum0(2);% + par.dcum0(3))/2;
d2 = par.dcum0(5);% + par.dcum0(4))/2;
d3 = par.dcum0(end);

%Voltage drops
Vint1 = (q*N1*w1.^2)./(2*epp1);
Vint2 = (q*N2*w2.^2)./(2*epp2);
Vint3 = (q*N3*w3.^2)./(2*epp3);

% Ionic resistance - asumed to be constant in this version
switch par.N_ionic_species
    case 1
        R = par.dcell(par.active_layer)./(par.e.*par.Ncat(par.active_layer).*par.mucat(par.active_layer));
    case 2
        R = par.dcell(par.active_layer)./(par.e.*(par.Ncat(par.active_layer).*par.mucat(par.active_layer) + par.Nani(par.active_layer).*par.muani(par.active_layer)));
end
% Capacitances
Cap.C1 = epp1/w1;
Cap.C2 = epp2/w2;
Cap.C3 = epp3/w3;

Cap.Cionic = Cap.C2;
% interface capacitances- called Cion in Moia et al. 2019
Cap.Cint1 = 1/(1/Cap.C1 + 1/Cap.C2);
Cap.Cint2 = 1/(1/Cap.C2 + 1/Cap.C3);
Cap.Cion = 1/(1/Cap.Cint1 + 1/Cap.Cint2);
Cap.CT_a = Cap.Cion;

% Perovskite PP capacitance
Cap.Cg_dash = epp2/par.d_active;
Cap.CT_b = 1/(1/Cap.Cion + 1/Cap.Cg_dash);
Cap.CT_c = Cap.Cion+Cap.Cg_dash;
Cap.CT_G = 1/((1/C1) + (1/((C2/2) + Cg_dash)) + (1/C3));

tau.RC_a = R*Cap.Cion;
tau.RC_b = R*Cap.CT_b;
tau.RC_c = R*Cap.CT_c;
tau.RC_g = R*Cap.CT_g;

Vbi = par.Vbi;

Q0 = (2^(1/2)*N1*N2*N3*epp1*epp2*epp3*q*((Vbi*(N1*N2*epp1*epp2 + 2*N1*N3*epp1*epp3 + N2*N3*epp2*epp3))...
    /(N1*N2*N3*epp1*epp2*epp3*q))^(1/2))/(N1*N2*epp1*epp2 + 2*N1*N3*epp1*epp3 + N2*N3*epp2*epp3);

%% Time array
t = meshgen_t(par);

%% Applied voltage function V(t)
Vapp_fun = fun_gen(par.V_fun_type);
Vapp = 0;

%% Solve the ODE for Q based on the ionic current through the bulk 'resistor' dQdt = Phi_per(Q(t))/Rion.
% V peroskite is simply the x-dependent terms for the Phi function (below) for
% d1+w2 < x < d2-w2
    function dQdt = odefun(t, Q)
        Vapp = Vapp_fun(par.V_fun_arg, t);
        V = Vbi - Vapp;
        % Second term on rhs is the Potential across the perovskite:
        dQdt = (1/R)*(((N1*epp1*Q + N3*epp3*Q)*(2*Q + N2*d1*q - N2*d2*q))/...
            (N2*epp2*q*(N1*epp1 + N3*epp3)) - ((N3*epp1*epp3*((N1*q*(2*N1*V*epp1*epp2^2 +...
            2*N3*V*epp2^2*epp3 - 2*N1*d1*epp1*epp2*Q + 2*N1*d2*epp1*epp2*Q - 2*N3*d1*epp2*epp3*Q +...
            2*N3*d2*epp2*epp3*Q - (2*N1*epp1*epp2*Q^2)/(N2*q) - (2*N3*epp2*epp3*Q^2)/(N2*q) +...
            N1*N3*d1^2*epp1*epp3*q + N1*N3*d2^2*epp1*epp3*q -...
            2*N1*N3*d1*d2*epp1*epp3*q))/(N3*epp1*epp3))^(1/2) + N1*N3*d1*epp1*epp3*q -...
            N1*N3*d2*epp1*epp3*q)*(2*Q + N2*d1*q - N2*d2*q))/(N2*epp2^2*q*(N1*epp1 + N3*epp3)));
    end

%% Call solver - default settings currently work well
options = odeset('MaxStep', 1e-5*abs(t(end)-t(1)));
[tout, Q] = ode15s(@odefun,t,Q0,options);
Q = real(Q');

rhomat = zeros(length(t), length(x));
Fmat = zeros(length(t), length(x));
Phimat = zeros(length(t), length(x));

%% Define peicewise functions for charge density rho, E-field F, and Potential Phi
% I couldn't find a convenient route to convert from the symbolic piecewise
% functions in depletion_approx_model8_analytical

    function rhosol = getrho
        rhosol = zeros(1, length(x));
        for ii = 1:length(x)
            if x(ii) <= d1 - (w1 + wapp1)
                rhosol(1,ii) = 0;
            elseif x(ii) <= d1 && d1 - (w1 + wapp1) < x(ii)
                rhosol(1,ii) = -q*N1;
            elseif x(ii) <= d1 + w2 && d1 < x(ii)
                rhosol(1,ii) = q*N2;
            elseif d1 + w2 < x(ii) && x(ii) <= d2 - w2
                rhosol(1,ii) = 0;
            elseif x(ii) <= d2 && d2 - w2 < x(ii)
                rhosol(ii) = -q*N2;
            elseif x(ii) <= d2 + (w3 + wapp3) && d2 < x(ii)
                rhosol(ii) = q*N3;
            elseif d2 + (w3 + wapp3) < x(ii) && x(ii) <= d3
                rhosol(ii) = 0;
            end
        end
    end

    function Fsol = getF
        Fsol = zeros(1, length(x));
        for ii = 1:length(x)
            if x(ii) <= d1 - (w1 + wapp1)
                Fsol(1,ii) = 0;
            elseif x(ii) <= d1 && d1 - (w1 + wapp1) < x(ii)
                Fsol(1,ii) = (-q*N1)*(x(ii)-(d1-(w1 + wapp1)))/epp1;
            elseif x(ii) <= d1 + w2 && d1 < x(ii)
                Fsol(1,ii) = (q*N2*(x(ii)-d1))/epp2 + (-q*N1)*(w1 + wapp1)/epp2;
            elseif d1 + w2 < x(ii) && x(ii) <= d2 - w2
                Fsol(1,ii) = (q*N2*w2/epp2)+((-q*N1)*(w1 + wapp1)/epp2);
            elseif x(ii) <= d2 && d2 - w2 < x(ii)
                Fsol(ii) = (-q*N2*(x(ii)-(d2-w2))/epp2)+((q*N2*w2)/epp2)+((-q*N1)*(w1 + wapp1)/epp2);
            elseif x(ii) <= d2 + (w3 + wapp3) && d2 < x(ii)
                Fsol(ii) = ((q*N3)*(x(ii)-d2)/epp3)+((-q*N2*w2)+(q*N2*w2)+(-q*N1)*(w1 + wapp1))/epp3;
            elseif d2 + (w3 + wapp3) < x(ii) && x(ii) <= d3
                Fsol(ii) = ((q*N3)*(w3 + wapp3)/epp3)+((-q*N2*w2)+(q*N2*w2)+(-q*N1)*(w1 + wapp1))/epp3;
            end
        end
    end

    function Phisol = getPhi
        %% Subfunction define the potential
        % This is currently written in terms of depletion widths and Qapp1
        % since this was the easiest approach in the first instance. These
        % could be substituted for the appropriate expresions to obtain
        % the expresions in terms of Q only.
        Phisol = zeros(1, length(x));
        
        for ii = 1:length(x)
            if x(ii) <= d1 - (w1 + wapp1)
                Phisol(1,ii) = 0;
            elseif x(ii) <= d1 && d1 - (w1 + wapp1) < x(ii)
                Phisol(1,ii) = -0.5.*(-q.*N1).*(x(ii)-(d1-(w1 + wapp1))).^2/epp1;
            elseif x(ii) <= d1 + w2 && d1 < x(ii)
                Phisol(1,ii) = -(0.5.*(q.*N2.*(x(ii)-d1).^2)/epp2 + (x(ii)-d1).*(-q.*N1).*(w1 + wapp1)/epp2) + (-0.5.*(-q.*N1).*(w1 + wapp1).^2/epp1);
            elseif d1 + w2 < x(ii) && x(ii) <= d2 - w2
                Phisol(1,ii) = -(((q.*N2.*w2)/epp2 + (-q.*N1).*(w1 + wapp1)/epp2).*(x(ii)-(d1+w2))) +...
                    -(0.5.*(q.*N2.*w2.^2)/epp2 + w2.*(-q.*N1).*(w1 + wapp1)/epp2) + (-0.5.*(-q.*N1).*(w1 + wapp1).^2/epp1);
            elseif x(ii) <= d2 && d2 - w2 < x(ii)
                Phisol(ii) = -((0.5.*-q.*N2.*(x(ii)-(d2-w2)).^2/epp2)+(((q.*N2.*w2)/epp2)+((-q.*N1).*(w1 + wapp1)/epp2)).*(x(ii)-(d2-w2))) +...
                    -(((q.*N2.*w2)/epp2 + (-q.*N1).*(w1 + wapp1)/epp2).*((d2-w2)-(d1+w2))) +...
                    -(0.5.*(q.*N2.*w2.^2)/epp2 + w2.*(-q.*N1).*(w1 + wapp1)/epp2) + (-0.5.*(-q.*N1).*(w1 + wapp1).^2/epp1);
            elseif x(ii) <= d2 + (w3 + wapp3) && d2 < x(ii)
                Phisol(ii) = -(0.5.*(q.*N3).*(x(ii)-d2).^2/epp3 + (((-q.*N2.*w2)+(q.*N2.*w2)+(-q.*N1).*(w1 + wapp1))/epp3).*(x(ii)-d2)) +...
                    -((0.5.*-q.*N2.*w2.^2/epp2)+(((q.*N2.*w2)/epp2)+((-q.*N1).*(w1 + wapp1)/epp2)).*w2) +...
                    -(((q.*N2.*w2)/epp2 + (-q.*N1).*(w1 + wapp1)/epp2).*((d2-w2)-(d1+w2))) +...
                    -(0.5.*(q.*N2.*w2.^2)/epp2 + w2.*(-q.*N1).*(w1 + wapp1)/epp2) + (-0.5.*(-q.*N1).*(w1 + wapp1).^2/epp1);
            elseif d2 + (w3 + wapp3) < x(ii) && x(ii) <= d3
                Phisol(ii) = -(((q.*N3).*(w3 + wapp3)/epp3  + ((-q.*N2.*w2)+(q.*N2.*w2)+(-q.*N1).*(w1 + wapp1)/epp3)).*(x(ii)-(d2 + (w3 + wapp3))))+...
                    -(0.5.*(q.*N3).*(w3 + wapp3).^2/epp3 + (((-q.*N2.*w2)+(q.*N2.*w2)+(-q.*N1).*(w1 + wapp1))/epp3).*(w3 + wapp3)) +...
                    -((0.5.*-q.*N2).*w2.^2/epp2 + (((q.*N2.*w2)/epp2)+((-q.*N1).*(w1 + wapp1)/epp2)).*w2) +...
                    -(((q.*N2.*w2)/epp2 + (-q.*N1).*(w1 + wapp1)/epp2).*((d2-w2)-(d1+w2))) +...
                    -(0.5.*(q.*N2.*w2.^2)/epp2 + w2.*(-q.*N1).*(w1 + wapp1)/epp2) + (-0.5.*(-q*N1)*(w1 + wapp1).^2/epp1);
            end
        end
    end

for i = 1:length(t)
    Vapp = Vapp_fun(par.V_fun_arg, t(i));
    V = Vbi - Vapp;
        
    w1 = Q(i)/(q*N1);
    w2 = Q(i)/(q*N2);
    w3 = Q(i)/(q*N3);
    wapp1 = (N3*epp1*epp3*((N1*q*(2*N1*V*epp1*epp2^2 + 2*N3*V*epp2^2*epp3 -...
        2*N1*N2*epp1*epp2*q*w2^2 - 2*N2*N3*epp2*epp3*q*w2^2 + N1*N3*d1^2*epp1*epp3*q +...
        N1*N3*d2^2*epp1*epp3*q - 2*N1*N2*d1*epp1*epp2*q*w2 + 2*N1*N2*d2*epp1*epp2*q*w2 -...
        2*N2*N3*d1*epp2*epp3*q*w2 + 2*N2*N3*d2*epp2*epp3*q*w2 - 2*N1*N3*d1*d2*epp1*epp3*q))/...
        (N3*epp1*epp3))^(1/2))/(N1*epp2*q*(N1*epp1 + N3*epp3)) -...
        (N1^2*epp1*epp2*w1 - N1*N3*d1*epp1*epp3 + N1*N3*d2*epp1*epp3 +...
        N1*N3*epp2*epp3*w1)/(N1*epp2*(N1*epp1 + N3*epp3));
    wapp3 = (epp1*epp3*((N1*q*(2*N1*V*epp1*epp2^2 + 2*N3*V*epp2^2*epp3 -...
        2*N1*N2*epp1*epp2*q*w2^2 - 2*N2*N3*epp2*epp3*q*w2^2 + N1*N3*d1^2*epp1*epp3*q +...
        N1*N3*d2^2*epp1*epp3*q - 2*N1*N2*d1*epp1*epp2*q*w2 + 2*N1*N2*d2*epp1*epp2*q*w2 -...
        2*N2*N3*d1*epp2*epp3*q*w2 + 2*N2*N3*d2*epp2*epp3*q*w2 - 2*N1*N3*d1*d2*epp1*epp3*q))/...
        (N3*epp1*epp3))^(1/2))/(epp2*q*(N1*epp1 + N3*epp3)) - ...
        (N1*d2*epp1*epp3 - N1*d1*epp1*epp3 + N1*epp1*epp2*w3 + N3*epp2*epp3*w3)...
        /(epp2*(N1*epp1 + N3*epp3));
    
    rhomat(i,:) = getrho;
    Fmat(i,:) = getF;
    Phimat(i,:) = getPhi;
end

% V2 = Vbi - Vapp - V1;
for i = 1:length(t)
    dVdx(i,:) = gradient(Phimat(i,:), x);
end

% Particle Current
Jpart = gradient(Q, tout);

% Displacement current;
Jdisp = par.e*par.epp(par.active_layer)*par.epp0.*(gradient(dVdx(:,pmid), tout));

J = Jpart + Jdisp';

if figson
    
    for i =1:length(tarr)
        
        try
                p2 = find(tout <= tarr(i));
                p2 = p2(end);
            
            % Inidividual figures
                figure(10)
                plot(xnm, rhomat(p2, :));
                xlabel('Position [nm]')
                ylabel('rho [cm-3]')
                xlim([xnm(1), xnm(end)])
                hold on
                grid on
            
                figure(11)
                plot(xnm, Fmat(p2, :));
                xlabel('Position [nm]')
                ylabel('E field [Vcm-1]')
                xlim([xnm(1), xnm(end)])
                hold on
                grid on
            
                figure(12)
                plot(xnm, Phimat(p2, :));
                xlabel('Position [nm]')
                ylabel('Potential [V]')
                xlim([xnm(1), xnm(end)])
                hold on
                grid on
            
            % Three panel figure
%             figure(22)
%             subplot(3, 1, 1)
%             plot(xnm, rhomat(p2, :)/q);
%             xlabel('Position [nm]')
%             ylabel('rho [cm-3]')
%             xlim([xnm(1), xnm(end)])
%             hold on
%             grid on
%             
%             subplot(3, 1, 2)
%             plot(xnm, Fmat(p2, :));
%             xlabel('Position [nm]')
%             ylabel('E field [Vcm-1]')
%             xlim([xnm(1), xnm(end)])
%             hold on
%             grid on
%             
%             subplot(3, 1, 3)
%             plot(xnm, Phimat(p2, :));
%             xlabel('Position [nm]')
%             ylabel('Potential [V]')
%             xlim([xnm(1), xnm(end)])
%             hold on
%             grid on
        catch
            warning('could not plot desired potential')
        end
    end
    figure(10)
    hold off
    figure(11)
    hold off
    figure(12)
    hold off
%     figure(22)
%     subplot(3, 1, 1)
%     hold off
%     subplot(3, 1, 2)
%     hold off
%     subplot(3, 1, 3)
%     hold off
    
    figure(17)
    plot(tout, Q);
    xlabel('Time [s]')
    ylabel('Charge [C cm-2]')
    hold off
    
    figure(18)
    plot(tout, J);
    xlabel('Time [s]')
    ylabel('Current density [J cm-2]')
    hold off
    
end

end
