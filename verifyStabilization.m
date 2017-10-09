function all_stable = verifyStabilization(sol_matrix, t, time_percentage) % verify if the tmax provided was enough

    names = ["electrons", "holes", "ions", "potential"];
    log = [true, true, false, false]; % which values to check when representing in log10 form or not
    all_stable = true;
    
    for i = 1:length(sol_matrix(1,1,:))

        [~,time_index]=min(abs(t-time_percentage*t(length(sol_matrix(:,1,i))))); % get time mesh index at 3/4 of maximum time reached

        profile_at_time=sol_matrix(time_index,:,i); % take profile of values at a certain time of evolution
        profile_end=sol_matrix(end,:,i); % take profile of values at the end of time

        if log(i) % for variables ranging on huge scales comparing the log values makes more sense
            profile_at_time=log10(profile_at_time);
            profile_end=log10(profile_end);
        end

        difference=sum(abs(profile_end-profile_at_time)); % sum up all the differences between the profiles
        threshold=1e-4*sum(abs(profile_end-mean(profile_end))); % sum up absolute values, ignore constant bias
        stable = difference <= threshold;

        if ~stable
            warning('pindrift:verifyStabilization','A tmax of %s s has not been enough for the %s distribution to reach stability. Consider trying with a greater tmax.', num2str(t(end)), names(i));
        end
        all_stable = all_stable && stable;
    end
end