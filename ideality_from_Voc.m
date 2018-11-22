function nid = ideality_from_Voc(structs_oc)
%IDEALITY_FROM_VOC - obtains the ideality factor from a fitting of VOC versus light intensity of stabilized solutions
% example:
%   nid = ideality_from_Voc(genIntStructs(ssol_i_eq_SR, 1, 0.001, 4, false))

Int_array = NaN(1, size(structs_oc, 2));
Voc_array = NaN(1, size(structs_oc, 2));

for i = 1:size(structs_oc, 2)
    Int_array(i) = structs_oc{1, i}.p.Int;
    structs_oc{1, i}.p.figson = 0;
    [V_temp, ~, ~, ~] = pinana(structs_oc{1, i});
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


