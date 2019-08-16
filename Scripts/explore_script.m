% An example script to demonstrate how to run a parameter exploration using
% the parallel computing toolbox
% Users will likely need to modify explore.explore2par
% Obtain the base parameters
par_ptpd = pc('input_files/ptpd_mapi_pcbm.csv');

% For the first example we will run JV and steady-state Voc for 3 different
% active layer thicknesses and light intensities
% ptpd_parex_dactive_light_33 = explore.explore2par(par_ptpd, {'dcell(1,4)','Int'},...
%     {[40e-7, 140e-7, 340e-7, 740e-7], logspace(-2,0,3)}, 200);
ptpd_parex_dactive_light_33 = explore.explore2par(par_ptpd, {'dcell(1,4)','Int'},...
    {[40e-7, 740e-7], [10]}, 200);

% Example plots
% 1 sun JV plot for different active layer thickness
explore.plotJV(ptpd_parex_dactive_light_33, [1,1,1], [0,1,0])
legend('d_{active} = 100 nm', 'd_{active} = 200 nm', 'd_{active} = 400 nm')

% Plot the reverse scan Voc as a function of thickness and light intensity as a surface
explore.plotsurf(ptpd_parex_dactive_light_33, 'Voc_r', 1, 0, 0)

% Plot the steady-state ion charge density profiles as a function of
% thickness at 1 sun
explore.plotprof_2D(ptpd_parex_dactive_light_33, 'a_f', [1,1,1], [0,1,0], 0,0)

% Save the workspace- this is commented out as the filepath should lead to
% a folder on your computer. It is not recommended to store large files in
% the Github repository as they will be cached even when deleted leading to
% a large folder size regardless of whether it contains data or not.
% save('/Users/Username/Data/temp.mat')
