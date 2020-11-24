function export_ELx(output_filename, sol, t)

%% Export Energy level diagram
% sol is a Driftfusion solution structure
% t is the time desired for output

% Calculate energy levels as a function of x and t
[Ecb, Evb, Efn, Efp] = dfana.QFLs(sol);

% find the time point
p1 = find(sol.t <= t);
p1 = p1(end);

% Extract the row for time t
Ecb = Ecb(p1, :);
Evb = Evb(p1, :);
Efn = Efn(p1, :);
Efp = Efp(p1, :);

% Create the table and write to file
headers = {'Position  (cm)', 'Ecb (eV)', 'Evb (eV)', 'Efn (eV)', 'Efp (eV)'};
data = num2cell([sol.x', Ecb', Evb', Efn', Efp']);

output_cell = [headers; data];
T = cell2table(output_cell);
writetable(T, ['./Output_files/', output_filename, '.txt'], 'Delimiter', 'tab', 'WriteVariableNames', 0);

end