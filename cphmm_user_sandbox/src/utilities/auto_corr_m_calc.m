
function auto_corr_m = auto_corr_m_calc(traces, max_delay)

% returns the average autocovariance of the traces

% ----------------------calculates the global mean-------------------------

    global_mean = 0;
    count = 0;
    for i = 1:length(traces)
        global_mean = global_mean + sum(traces{i});
        count = count + length(traces{i});
    end
    global_mean = global_mean / count;
    
%--------------calculates autocovariance using auto_corr_r_calc------------
    new_traces = cell([1 length(traces)]);
    for i = 1:length(traces)
        new_traces{i} = traces{i} - global_mean;
    end
    auto_corr_m = auto_corr_r_calc(new_traces, max_delay);

end
