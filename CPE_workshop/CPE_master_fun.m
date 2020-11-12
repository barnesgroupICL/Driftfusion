function sol = CPE_master_fun(params_filepath, output_filename, light_intensities, Vmax, Vmin)

% Created for the Centre for Processable Electronics (CPE) MRes course
% on semiconductor physics
% Authors: P.R.F. Barnes, P. Calado 2020
% Department of Physics, Imperial College London

%% Obtain cyclic voltammogram (CV) solutions and conductivity vs voltage (sigma)
sol = sigma_light(params_filepath, light_intensities, Vmax, Vmin, output_filename);


end