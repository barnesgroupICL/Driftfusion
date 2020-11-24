function export_solution(output_filename, sol, tplot)

%% Export Energy level diagram
% sol is a Driftfusion solution structure
% Vapp is the voltage for the desired output

% find the time point
p1 = find(sol.t <= tplot);
p1 = p1(end);

% Get the solution at desired time point
V = sol.u(p1,:,1);
n = sol.u(p1,:,2);
p = sol.u(p1,:,3);

% Create the table and write to file
headers = {'Position  (cm)', 'V (V)', 'n (cm-3)', 'p (cm-3)'};
data = num2cell([sol.x', V', n', p']);

output_cell = [headers; data];
T = cell2table(output_cell);
writetable(T, ['./Output_files/', output_filename, '.txt'], 'Delimiter', 'tab', 'WriteVariableNames', 0);

end