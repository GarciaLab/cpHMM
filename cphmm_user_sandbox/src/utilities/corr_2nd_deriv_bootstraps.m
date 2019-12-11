function stds = corr_2nd_deriv_bootstraps(traces, max_delay, num_times, type)

% takes random selections of traces and calculates the standard deviation
% of the autocorrelations

    vals = cell([1 max_delay - 2]);

    for i = 1:num_times
        sample_idx = randi([1 length(traces)], 1, length(traces));
        samples = cell([1 length(traces)]);
        for j = 1:length(samples)
            samples{j} = traces{sample_idx(j)};
        end

        if type == 'r'
            corr = auto_corr_r_calc(samples, max_delay);
        else
            corr = auto_corr_m_calc(samples, max_delay);
        end
        c1st_deriv = corr(2:max_delay) - corr(1:max_delay - 1);
        c2nd_deriv = c1st_deriv(2:max_delay-1) - c1st_deriv(1:max_delay-2);
        for j = 1:length(c2nd_deriv)
            vals{j}(i) = c2nd_deriv(j);
        end
    end
    stds = zeros([1 length(vals)]);
    for i = 1:length(vals)
        stds(i) = std(vals{i});
    end
    %stds = vals;
end

