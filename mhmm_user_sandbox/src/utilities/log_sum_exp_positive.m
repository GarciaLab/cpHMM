function logsum = log_sum_exp_positive (arr)
    % Calculates the sum of logs.
    % 
    % INPUTS
    % arr: list of log values
    % 
    % OUTPUTS
    % logsum: the sum of log values. Note: underflow issues are avoided.
    
    arr_max = max(arr(:)); 
    logsum = arr_max + log(sum(exp(arr(:)-arr_max)));
    if (isnan(logsum))
        logsum = -Inf;        
    end