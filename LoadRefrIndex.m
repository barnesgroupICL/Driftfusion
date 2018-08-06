% Function LoadRefrIndex
% This function returns the complex index of refraction spectra, ntotal, for the
% material called 'name' for each wavelength value in the wavelength vector
% 'wavelengths'.  The material must be present in the index of refraction
% library 'Index_of_Refraction_library.xls'.  The program uses linear
% interpolation/extrapolation to determine the index of refraction for
% wavelengths not listed in the library.
function [n_interp, k_interp] = LoadRefrIndex(name,wavelengths)

%Data in IndRefr, Column names in IndRefr_names
[IndRefr,IndRefr_names]=xlsread('Index_of_Refraction_library.xls');

% Load index of refraction data in spread sheet, will crash if misspelled
file_wavelengths=IndRefr(:,strmatch('Wavelength',IndRefr_names));
n=IndRefr(:,strmatch(strcat(name,'_n'),IndRefr_names));
k=IndRefr(:,strmatch(strcat(name,'_k'),IndRefr_names));

% Interpolate/Extrapolate data linearly to desired wavelengths
n_interp=interp1(file_wavelengths, n, wavelengths, 'linear', 'extrap');
k_interp=interp1(file_wavelengths, k, wavelengths, 'linear', 'extrap');

%Return interpolated complex index of refraction data
%ntotal = n_interp+1i*k_interp; 

end