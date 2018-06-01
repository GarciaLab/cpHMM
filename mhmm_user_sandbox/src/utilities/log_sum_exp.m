function logsum = log_sum_exp (arr, signs)
    % Calculates the sum of logs with specified signs.
    % 
    % INPUTS
    % arr: list of log values
    % sign: list of signs corresponding to each log values
    % 
    % OUTPUTS
    % logsum: a 2-element array of the log and the sign of the sum.
    %         Note: the calculation avoids underflow issues.
    
    arr_max = max(arr(:));
    term2_array = signs.*exp(arr-arr_max);
    term2 = sum(term2_array(:));
    logsum = [arr_max + log(abs(term2)), sign(term2)];
    if (isnan(logsum(1)))
        logsum = [-Inf, 0];        
    end