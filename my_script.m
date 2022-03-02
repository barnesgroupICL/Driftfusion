%% Conductivity profiles
% So systematically you could look at the following.

% 1)Electrode workfunctions
% 2)Ion density 
% 3)1 ion, opposite charge (i.e. mobile anions)
% 4)2 ions
% 5)Different ion densities

% Ideally you would set some of these parameter explorations up as loops and extract peak conductivity then plot 
% on a contour plot with x = Ion density, y = Electrode workfunctions, z = peak conductivity for example.

%% Initialize driftfusion
initialise_df

%% Add parameter file to path 
% Filepath Mac
par_alox = pc('Input_files/alox.csv');
no_of_diff_ion_conc=abs(log10((par_alox.Ncat(1,3)/1e18)));
epoints=round((par_alox.Phi_left-par_alox.Phi_right)/-(0.1));%number of different electrode values
valuestore_n=zeros(epoints,no_of_diff_ion_conc);%create the matrix for storing n values
valuestore_p=zeros(epoints,no_of_diff_ion_conc);%create the matrix for storing p values
phi_left_electrode = par_alox.Phi_left;
row=1; %intialize 
column=1;
%% while
while par_alox.Ncat(1,3)>1e17
%% Equilibrium solutions 
 
 for electrode_change= par_alox.Phi_left:0.1:par_alox.Phi_right %loop to run for different electrode workfunction

 soleq_alox = equilibrate(par_alox);
 
%% Plot equilibrium energy level diagram

% dfplot.ELnpx(soleq_alox.ion) commented out since many schematics will be
% formed. Similiar with other plots



%% Current-voltage scan
% JVsol = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
% JVsol = doJV(soleq_sio2.ion, 100e-3, 201, 1, 0, 0, 1, 1);
k_scan = 0.1;
Vmax = 1.2;
Vmin = -1.2;

% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV = doCV(soleq_alox.ion, 0, 0, Vmax, Vmin, k_scan, 1, 241);
%% Plot Vapp vs time
% dfplot.Vappt(sol_CV)

%% Plot JV scan
%dfplot.JtotVapp(sol_CV, 0);
%set(gca,'YScale','log')

%% Plot anion and cation densities
%dfplot.acx(sol_CV, 1/k_scan*[0:Vmax/3:Vmax]);

%% Plot electron and hole profiles
%dfplot.npx(sol_CV, 1/k_scan*[0:Vmax/3:Vmax]);

%% Plot space charge density
%dfplot.rhox(sol_CV, 1/k_scan*[0:Vmax/3:Vmax]);

%% Calculate conductivity
[sigma_n, sigma_p] = dfana.calc_conductivity(sol_CV);

%% Debye length Calculation
par = par_alox;

e = par.e;
V_T = par.kB*par.T;                     % Thermal votlage
epp_pvsk = e*par.epp0*par.epp(3);       % Perovskite absolute dielectric constant
N0 = par.Ncat(3);                   
N0_courtier = 1.6e19;                   % cm-3

L_D = sqrt((epp_pvsk*V_T)/(e*N0));      % Deby width [cm]
L_D_courtier = sqrt((epp_pvsk*V_T)/(e*N0_courtier));

N_Debye = 3;                            % Number of Debye lengths to average electron density over
%%
x_perov_left = sol_CV.par.dcum0(3);
x_perov_right = sol_CV.par.dcum0(4);

x = sol_CV.x;
t = sol_CV.t;
Vappt = dfana.calcVapp(sol_CV);
%% Get point at which perovskite starts 
sigma_n_bar = mean(sigma_n(:, x > x_perov_left & x < x_perov_left + N_Debye*L_D), 2);
sigma_p_bar = mean(sigma_p(:, x > x_perov_left & x < x_perov_left + N_Debye*L_D), 2);
% 
% sigma_n_bar_entire = mean(sigma_n(:, x > x_perov_left & x < x_perov_right), 2);
% sigma_p_bar_entire = mean(sigma_p(:, x > x_perov_left & x < x_perov_right), 2);
% 
% sigma_n_bar_bulk = mean(sigma_n(:, x > x_perov_left + N_Debye*L_D & x < x_perov_right), 2);
% sigma_p_bar_bulk = mean(sigma_p(:, x > x_perov_left + N_Debye*L_D & x < x_perov_right), 2);

%% Find peak conductivity for applied bias
pp_Vmax = find(Vappt == max(Vappt));      %% pp = point position
pp_Vmin = find(Vappt == min(Vappt));      %% pp = point position

%% Max
sigma_n_bar_Vpeak = sigma_n_bar(pp_Vmax);
sigma_p_bar_Vpeak = sigma_p_bar(pp_Vmax);
%% Put value inside matrix
valuestore_n(row,column)= sigma_n_bar_Vpeak;  
valuestore_p(row,column)= sigma_p_bar_Vpeak;
row=row+1;
par_alox.Phi_left=par_alox.Phi_left-0.1;
 end
row=1; %reset back to row 1
par_alox.Phi_left=phi_left_electrode; %reset back to original eletrode value
column=column+1; %move to the next doping value
% %% Plot the outputs
% % figure(101)
% % Ntr = 6;            % Number of voltage transients
% % 
% % for i = 1:Ntr
% %     plot(sol_VTROTTR(i).t, DeltaVoc(i,:));
% %     hold on
% % end
% % 
% % xlabel('Time [s]')
% % ylabel('DeltaV [V]')
% % xlim([0, 1e-5])
% % hold off
% %%
% % %% Plot average conductivity
% % figure
% % plot(Vappt, sigma_n_bar, Vappt, sigma_p_bar)
% % axis([-1 1 0 inf])
% % xlabel('Voltage [V]')
% % ylabel('Average conductivity [Siemens]')
% % legend('Electron', 'Hole')
%  
% %% Plot average conductivity
% % figure(200)
% % semilogy(Vappt, sigma_n_bar, Vappt, sigma_p_bar)
% % xlabel('Voltage [V]')
% % ylabel('Average channel conductivity [Semilog]')
% % legend('Electron', 'Hole')
% % 
% % %% Plot average conductivity
% % % figure(201)
% % % plot(Vappt, sigma_n_bar, Vappt, sigma_p_bar)
% % % xlabel('Voltage [V]')
% % % ylabel('Average channel conductivity [Linear]')
% % % legend('Electron', 'Hole')
% % 
% % %% Plot average conductivity
% % figure(202)
% % semilogy(Vappt, sigma_n_bar_bulk, Vappt, sigma_p_bar_bulk)
% % xlabel('Voltage [V]')
% % ylabel('Average bulk conductivity [Semilog]')
% % legend('Electron', 'Hole')
% % 
% % %% Plot average conductivity
% % % figure(203)
% % % plot(Vappt, sigma_n_bar_bulk, Vappt, sigma_p_bar_bulk)
% % % xlabel('Voltage [V]')
% % % ylabel('Average bulk conductivity [Linear]')
% % % legend('Electron', 'Hole')
% % 
% % %% Plot average conductivity
% % figure(204)
% % semilogy(Vappt, sigma_n_bar_entire, Vappt, sigma_p_bar_entire)
% % xlabel('Voltage [V]')
% % ylabel('Average entire conductivity [Semilog]')
% % legend('Electron', 'Hole')
% % 
% % %% Plot average conductivity
% % % figure(205)
% % % plot(Vappt, sigma_n_bar_entire, Vappt, sigma_p_bar_entire)
% % % xlabel('Voltage [V]')
% % % ylabel('Average entire conductivity [Linear]')
% % % legend('Electron', 'Hole')
% % 
% % 

 par_alox.Ncat(1,3)= par_alox.Ncat(1,3)/10;
 end

 %%
figure
contour(valuestore_n)
contour(valuestore_p)
%%
% %% Make movie for anions and cations
% %makemovie(sol_CV, @dfplot.acx, 0, [0, 1.5e18], 'acx', true, true);