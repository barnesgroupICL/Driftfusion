function nid = ideality_from_Voc(structs_oc)
%IDEALITY_FROM_VOC - obtains the ideality factor from a fitting of VOC versus light intensity of stabilized solutions contained in a cell
% For generating the needed input, genIntStructs can be used.
%
% Syntax:  nid = ideality_from_Voc(structs_oc)
%
% Inputs:
%   STRUCTS_OC - a cell containing structures at open circuit and various
%     light intensities. The open circuit can be either obtained from
%     symmetric solutions (from genIntStructs) or from asymmetric solutions
%     with applied voltage (from genIntStructRealVoc)
%
% Outputs:
%   NID - the numeric ideality factor
%
% Example:
%   nid = ideality_from_Voc(genIntStructs(ssol_i_eq_SR, 1, 0.001, 4, false))
%     
%
% Other m-files required: pinana
% Subfunctions: none
% MAT-files required: none
%
% See also pindrift.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: May 2018

%------------- BEGIN CODE --------------

Int_array = NaN(1, size(structs_oc, 2));
Voc_array = NaN(1, size(structs_oc, 2));

for i = 1:size(structs_oc, 2)
    Int_array(i) = structs_oc{1, i}.p.Int;
    structs_oc{1, i}.p.figson = 0;
    [V_temp, ~, ~] = pinana(structs_oc{1, i});
    Voc_array(i) = V_temp(end);
end

disp(Int_array)
disp(Voc_array)

% removes the dark point
nonzero_index = find(Int_array);

p = structs_oc{1, 1}.p;
lm = fitlm(log(Int_array(nonzero_index)), Voc_array(nonzero_index));
nid = lm.Coefficients.Estimate(2) * p.q/(p.kB*p.T);

figure(1)
    hold off
    plot(lm)
    ylabel('VOC')
    xlabel('log(light intensity)')

%------------- END OF CODE --------------

