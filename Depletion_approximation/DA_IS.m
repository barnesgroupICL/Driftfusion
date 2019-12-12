function [Cap, Farr] = DA_IS(par, Farr, Vac, periods, ppperiod)

coeff = zeros(3, length(Farr));

for i = 1:length(Farr)
    V_fun_arg = [par.Vapp, Vac, Farr(i), 0];        % [Vdc, Vac, Freq, Phase];
    coeff(:, i) = DA_sin(par, V_fun_arg, periods, ppperiod, 0, 0);
    
end

Zmag = Vac./coeff(2,:);
Cap = -sin(coeff(3,:))./(2*pi*Farr.*Zmag);

figure(600)
semilogx(Farr, Cap);
xlabel('Frequency [Hz]')
ylabel('Apparent Capacitance [Fcm-2]')

end