%Create Weighted Average Autocorrelation
function [wt_autocorr, a_boot_errors, wt_dd, dd_boot_errors] = weighted_autocorrelation(traces, lags, bootstrap,n_boots,trace_weights)
    %traces: array of traces with zeros preceeding and succeeding period of
    %        activity. Oriented column-wise
    %lags: num lags to use
    %bootstrap: binary var. If 1, returns array of bootstrap errors
    auto_array = zeros(lags+1,size(traces,2));
    time_steps = zeros(lags+1,size(traces,2));
    %Convert NaNs to zeros
    traces(isnan(traces)) = 0;
    t_filter = sum(traces>0)>1;
    traces = traces(:,t_filter);
    trace_weights = trace_weights(t_filter);
    if ~bootstrap
        n_boots = 1;
    end
    samples = zeros(lags+1,n_boots);    
    dd_samples = zeros(lags-1,n_boots);
    for b = 1:n_boots
        %If bootstrap errors are desired, take n_boots samples with
        %replacement, each the size of original array
        if bootstrap
            s_vec = 1:size(traces,2);
            s = randsample(s_vec,length(s_vec),true,trace_weights);
            trace_sample = traces(:,s);
        else
            trace_sample = traces;
        end        
        for col = 1:size(trace_sample,2)            
            trace = trace_sample(:,col);
            if sum(isnan(trace)) > 0 
                error('Problem with NaN filtering');
            end            
            %Isolate active portion
            trace_active = trace(find(trace,1):find(trace,1,'last'));
            if length(trace_active) < lags + 1
                warning('Length of input trace insufficient for specified number of lags')
                t_lags = length(trace_active) - 1;
            else
                t_lags = lags;
            end            
            auto_array(1:t_lags+1,col) = autocorr(trace_active,t_lags);
            time_steps(1:t_lags+1,col) = fliplr((length(trace_active)-t_lags):length(trace_active));
        end
        %Take weighted mean. Traces with length < lags should just not be
        %factored in for lags beyond their length
        numerator = sum(auto_array.*time_steps,2);
        denominator = sum(time_steps,2);
        samples(:,b) = numerator ./ denominator;
        dd_samples(:,b) = diff(diff(samples(:,b)));
    end
wt_autocorr = mean(samples,2);
wt_dd = mean(dd_samples,2);
a_boot_errors = std(samples')';
dd_boot_errors = std(dd_samples')';