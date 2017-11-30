function all_stable = verifyStabilization(sol_matrix, t, time_fraction) % verify if the tmax provided was enough
% time_fraction has to be strictly greater than zero and strictly less than 1

names = ["electrons", "holes", "ions", "potential"];
log = [true, true, false, false]; % which values to check when representing in log10 form or not
all_stable = true;

% no need to calculate end_time for each of the 4 solutions: if they
% break they break at the same time
end_time = t(length(sol_matrix(:, 1, 1))); % using t(end) could get a time out of the solution when the computation broke before reaching the final time
[~, time_index] = min(abs(t - time_fraction * end_time)); % get time mesh index at a specified percentage of maximum time reached

for i = 1:length(sol_matrix(1, 1, :)) % for each component of the solution verify that didn't change significantly over the requested time
    profile_at_time = sol_matrix(time_index, :, i); % take profile of values at a certain time of evolution
    profile_end = sol_matrix(end, :, i); % take profile of values at the end of time
    if log(i) % for variables ranging on huge scales comparing the log values makes more sense
        profile_at_time = log10(profile_at_time);
        profile_end = log10(profile_end);
    end

    difference = sum(abs(profile_end - profile_at_time)); % sum up all the differences between the profiles
    threshold = 1e-4 * sum(abs(profile_end - mean(profile_end))); % sum up absolute values, ignore constant bias
    stable = difference <= threshold;

    if ~stable
        warning('pindrift:verifyStabilization', 'Comparing final solutions at %s s and %s s showed that the %s distribution did not reach stability. Consider trying with a greater tmax.', num2str(t(time_index)), num2str(t(end)), names(i));
    end
    all_stable = all_stable && stable;
end