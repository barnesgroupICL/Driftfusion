%% Initialise Driftfusion
initialise_df

%% Build a base parameter set object
par_tio2_base = pc('Input_files/spiro_mapi_tio2.csv');

%% copy to a new object where we will chaneg the recombination parameters
par_tio2_low_rec = par_tio2_base;

%% Change the interfacial recombination parameters for interface 2 (layer no. 4)
par_tio2_low_rec.taun(4) = 1e-9;
par_tio2_low_rec.taup(4) = 1e-9;

%% Refresh the device- rebuilds the device structures par.dev and par.dev_ihalf etc.
par_tio2_low_rec = refresh_device(par_tio2_low_rec);

%% Get equilibrium solutions
soleq_tio2_base = equilibrate(par_tio2_base);
soleq_tio2_low_rec = equilibrate(par_tio2_low_rec);

%% DoJVs 0 - 1.2 V
JV_tio2_base = doJV(soleq_tio2_base.ion, 50e-3, 200, 1, 1, 0, 1.2, 3);
JV_tio2_low_rec = doJV(soleq_tio2_low_rec.ion, 50e-3, 200, 1, 1, 0, 1.2, 3);

%% Plot JVs
dfplot.JV(JV_tio2_base, 3)
figure(4)
hold on
dfplot.JV(JV_tio2_low_rec, 3)
ylim([-25e-3, 10e-3])
legend('Base, dk F', 'Base, dk R', 'Base, 1sun F', 'Base, 1sun R',...
    'Low rec, dk F', 'Low rec, dk R', 'Low rec, 1sun F', 'Low rec, 1sun R')
hold off
