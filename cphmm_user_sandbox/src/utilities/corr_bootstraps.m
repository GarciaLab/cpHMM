% takes random selections of traces and calculates the standard deviation
% of the autocorrelations
function stds = corr_bootstraps(traces, max_delay, num_times, type)

    vals = cell([1 max_delay]);

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
        for j = 1:length(corr)
            vals{j}(i) = corr(j);
        end
    end
    stds = zeros([1 length(vals)]);
    for i = 1:length(vals)
        stds(i) = std(vals{i});
    end
end
