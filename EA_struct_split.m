function [one, two, three] = EA_struct_split(in)

% sol_name: 'soleq_pedot_ion_stab_Int_0'
%                            Vdc: [3×1 double]
%                        periods: 20
%                           Freq: [3×61 double]
%                        tpoints: 1601
%                           tmax: [3×61 double]
%                            Int: [3×1 double]
%                         deltaV: 0.8000
%                      sun_index: 3
%                 AC_ExDC_E_bias: [3×61 double]
%                  AC_ExDC_E_amp: [3×61 double]
%                AC_ExDC_E_phase: [3×61 double]
%               AC_ExDC_E_i_bias: [3×61 double]
%                AC_ExDC_E_i_amp: [3×61 double]
%              AC_ExDC_E_i_phase: [3×61 double]
%                AC_Efield2_bias: [3×61 double]
%                 AC_Efield2_amp: [3×61 double]
%               AC_Efield2_phase: [3×61 double]
%              AC_Efield2_i_bias: [3×61 double]
%               AC_Efield2_i_amp: [3×61 double]
%             AC_Efield2_i_phase: [3×61 double]
%     AC_Efield_amp_squared_mean: [3×61 double]


list = ["one", "two", "three"];
names = fieldnames(in);

for i = 1:length(names)
    a=strcat('content = in.', names(i), ';');
    eval(a{1});
    if size(content,1) == 3
        disp(names(i))
        for j = 1:length(list)
            b=strcat(list(j), '.', names(i), ' = content(', num2str(j), ',:);');
            eval(b{1});
        end
    end
end

one.periods = in.periods;
two.periods = in.periods;
three.periods = in.periods;
one.tpoints = in.tpoints;
two.tpoints = in.tpoints;
three.tpoints = in.tpoints;
one.deltaV = in.deltaV;
two.deltaV = in.deltaV;
three.deltaV = in.deltaV;
one.sol_name = in.sol_name;
two.sol_name = in.sol_name;
three.sol_name = in.sol_name;

