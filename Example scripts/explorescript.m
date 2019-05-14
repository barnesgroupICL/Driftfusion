% An example script to demonstrate how to run a parameter exploration using
% the parallel computing toolbox
% Users will likely need to 

% Obtain the base parameters
par.ptpd = pc('input_files/ptpd.csv');

% For the first example we will run JV and steady-state Voc for 3 different
% active layer thicknesses and light intensities 

%ptpd_parex_dactive_light_33 = explore.explore2par(par.ptpd_Na18, {'dcell(1,4)','Int'},...
%    {[40e-7, 140e-7, 340e-7], logspace(-1,1,3)}, 200);
ptpd_parex_dETL_NDETL_light_33 = explore.explore2par(par.ptpd_Na18, {'dcell(1,7)','E0(7)'},...
    {[40e-7, 60e-7, 80e-7], [-3.7,-3.8,-3.9]}, 200);

save('ws_explore.mat')  % Save workspace