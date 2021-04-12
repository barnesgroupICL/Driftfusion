function EA_IS_struct_split(in)

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



% >> IS_pedot_sc_20mV
% 
% IS_pedot_sc_20mV = 
% 
%   struct with fields:
% 
%                   sol_name: 'soleq_pedot_ion_stab_Int_0'
%                        Vdc: [3×1 double]
%                    periods: 20
%                       Freq: [3×61 double]
%                    tpoints: 801
%                       tmax: [3×61 double]
%                        Int: [3×1 double]
%                     deltaV: 0.0200
%                  sun_index: 3
%                  Jtot_bias: [3×61 double]
%                   Jtot_amp: [3×61 double]
%                 Jtot_phase: [3×61 double]
%              ion_disp_bias: [3×61 double]
%               ion_disp_amp: [3×61 double]
%             ion_disp_phase: [3×61 double]
%              cat_disp_bias: [3×61 double]
%               cat_disp_amp: [3×61 double]
%             cat_disp_phase: [3×61 double]
%              ani_disp_bias: [3×61 double]
%               ani_disp_amp: [3×61 double]
%             ani_disp_phase: [3×61 double]
%                     r_bias: [3×61 double]
%                      r_amp: [3×61 double]
%                    r_phase: [3×61 double]
%                 np_dt_bias: [3×61 double]
%                  np_dt_amp: [3×61 double]
%                np_dt_phase: [3×61 double]
%                        cap: [3×61 double]
%              impedance_abs: [3×61 double]
%               impedance_im: [3×61 double]
%               impedance_re: [3×61 double]
%               cap_ion_disp: [3×61 double]
%     impedance_ion_disp_abs: [3×61 double]
%      impedance_ion_disp_im: [3×61 double]
%      impedance_ion_disp_re: [3×61 double]
%               cap_cat_disp: [3×61 double]
%     impedance_cat_disp_abs: [3×61 double]
%      impedance_cat_disp_im: [3×61 double]
%      impedance_cat_disp_re: [3×61 double]
%               cap_ani_disp: [3×61 double]
%     impedance_ani_disp_abs: [3×61 double]
%      impedance_ani_disp_im: [3×61 double]
%      impedance_ani_disp_re: [3×61 double]
%                      cap_r: [3×61 double]
%            impedance_r_abs: [3×61 double]
%             impedance_r_im: [3×61 double]
%             impedance_r_re: [3×61 double]
%                  cap_np_dt: [3×61 double]
%        impedance_np_dt_abs: [3×61 double]
%         impedance_np_dt_im: [3×61 double]
%         impedance_np_dt_re: [3×61 double]



list = ["dark", "sun01", "sun1"];
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

dark.periods = in.periods;
sun01.periods = in.periods;
sun1.periods = in.periods;
dark.tpoints = in.tpoints;
sun01.tpoints = in.tpoints;
sun1.tpoints = in.tpoints;
dark.deltaV = in.deltaV;
sun01.deltaV = in.deltaV;
sun1.deltaV = in.deltaV;
dark.sol_name = in.sol_name;
sun01.sol_name = in.sol_name;
sun1.sol_name = in.sol_name;

assignin('base', [inputname(1) '_dark'], dark);
assignin('base', [inputname(1) '_01sun'], sun01);
assignin('base', [inputname(1) '_1sun'], sun1);



