function [G0_Arr, k_rad_Arr, R0_Arr, Eg, VocArr, DeltaVoc] = Eg_vs_Voc(EgArr, Jsc_vs_Eg)

Eg = 0.6:0.1:3;

for i=1:length(Eg)
    
    [R0, k_rad, Voc, G0] = calcR0(Eg(i), EgArr, Jsc_vs_Eg);
    VocArr(i) = Voc;
    DeltaVoc(i) = Eg(i) - VocArr(i);  
    R0_Arr(i) = R0;
    k_rad_Arr(i) = k_rad;
    G0_Arr = G0;
        
end

figure(1111)
plot(Eg, VocArr, Eg, DeltaVoc)
xlabel('Eg [eV]')
ylabel('Voc [V]')

end