function logL_tot = likelihood_reduced_memory (fluo_values, v, noise, ...
                         pi0_log, A_log, K, w, kappa)

    % Returns the log likelihood of observing fluorescence data (single or
    % multiple traces), given model parameters.
    % 
    % INPUTS
    % fluo_values: cell data type with multiple time series fluorescence
    %              values that come from the same AP position
    % v: initial emission values
    % noise: initial gaussian noise
    % pi0_log: log of the initial naive state cmf at t=0
    % A_log: log of the initial transition probability matrix
    % K: number of naive states
    % w: memory
    % kappa: length of the MS2 loop in time steps
    %
    % OUTPUTS
    % logL_tot: log likelihood of observing the fluorescence data
    
    %%%%%%%%%%%%%%%%%%%%%%% Variable assignments %%%%%%%%%%%%%%%%%%%%%%%
    
    % number of traces
    n_traces = length(fluo_values);
    
    % logs of v and noise
    v_logs = log(abs(v));
    v_signs = sign(v);
    lambda_log = -2*log(noise);
    
    % logs of the fluorescence values &
    % the lengths of the fluorescence trace
    fluo_logs = cell([n_traces, 1]);
    fluo_signs = cell([n_traces, 1]);
    fluo_lengths = cell([n_traces, 1]);
    for i_tr = 1:n_traces
        fluo_logs{i_tr} = log(abs(fluo_values{i_tr}));
        fluo_signs{i_tr} = sign(fluo_values{i_tr});
        fluo_lengths{i_tr} = length(fluo_values{i_tr});
    end
    
    % list of states that can transition to each of 1:K^w compounds states
    allowed_to_list = zeros(K^w, K, 'int16');
    
    % list of the first digits (naive states) of each compound state
    digit_first_list = zeros(1, K^w, 'int16');
    
    % list of the number of each naive state encountered in compound states
    naive_count_list_MS2 = zeros(K^w, K);
    
    % mean fractions of MS2 loops transcribed at each elongation step
    ms2_coeff = ms2_loading_coeff(kappa,w);
    
    for i = 1:K^w
        allowed_to_list(i,:) = allowed_to(i, K, w);
        digit_first_list(i) = digit(i, 1, K, w);
        
        naive = compound_to_naive(i, K, w);
        for k = 1:K
            naive_count_list_MS2(i,k) = sum(ms2_coeff(naive == k));
        end
    end
    
    % variable used to account for adjustements at time points 1:(w-1)
    count_reduction = zeros([1,w-1]);
    for t = 1:(w-1)
        count_reduction(t) = sum(ms2_coeff((t+1):end));
    end
    
    % log of the naive count list
    naive_count_list_MS2_log = log(naive_count_list_MS2);
    
    % list of compound states that have a unique sequence of the number of
    % naive states contained in the base-K representation
    unique_naive_MS2_list = unique(naive_count_list_MS2,'rows');
    [n_unique, ~] = size(unique_naive_MS2_list);    
    
    % a list that maps the unique naive count combinations with compound
    % states
    naive_count_map = cell(n_unique, 1);   
    for i = 1:n_unique
        naive = unique_naive_MS2_list(i,:);
        naive_count_map{i} = find(ismember(naive_count_list_MS2,naive,'rows'));        
    end
    
    % list of possible compound states at all times
    time_max = max(cell2mat(fluo_lengths));
    possible_states_list = cell(time_max, 1);    
    for t = 1:time_max
        possible_states_list{t} = possible_states (t, K, w);        
    end
    
    % log likelihood of observing the fluorescence data
    logL_tot = 0;
    
    % row and column subscripts used for indexing A_log elements for
    % alpha matrix calculation
    A_log_alpha_rowSubs = digit_first_list(repmat(1:K^w, [1, K]));
    A_log_alpha_colSubs = digit_first_list(allowed_to_list(:));
    
    % A_log element indexing in a 2d slice used for alpha matrix
    % calculation
    A_log_alpha_subs_1d = A_log(sub2ind([K, K], A_log_alpha_rowSubs, ...
                                A_log_alpha_colSubs));
    alpha_A_log_list = reshape(A_log_alpha_subs_1d, [K^w K]);
    
    %%%%%%%%%%%%%%%%%%%%%% Likelihood calculation %%%%%%%%%%%%%%%%%%%%%%
    for i_tr = 1:n_traces
        
        % pre-allocation of the alpha coefficient matrix
        alpha_matrix_log = -Inf(K^w,fluo_lengths{i_tr});
        
        % list of log (X_t - V_t)^2 terms that appear in the emission pdf's
        difference_list_temp = zeros(K^w, fluo_lengths{i_tr});
        for i = 1:n_unique
            states = naive_count_map{i};
            difference_list_temp(states, :) = ...
                   repmat(difference_reduced_memory(fluo_logs{i_tr}, ...
                   fluo_signs{i_tr}, fluo_lengths{i_tr}, K, w, ...
                   count_reduction, states(1), v_logs, v_signs, ...
                   naive_count_list_MS2_log)', [length(states), 1]);
        end
        difference_list = difference_list_temp;
        eta_log_list = 0.5*(lambda_log - log(2*pi)) ...
            -0.5*exp(lambda_log + difference_list);
        
        % calculates the alpha matrix elements for t = 1
        for i = possible_states_list{1}
            alpha_matrix_log(i, 1) = eta_log_list(i, 1) + ...
                pi0_log(digit_first_list(i));
        end
        
        % calculates the alpha matrix elements for t > 1
        for t = 2:fluo_lengths{i_tr}
            % possible states at time t
            i_possible = possible_states_list{t};
            
            % list of terms that are added to find the alpha matrix
            % elements
            alpha_terms_list = alpha_A_log_list(i_possible, :) + ...
                reshape(alpha_matrix_log(allowed_to_list(i_possible,:),t-1), [], K);
            
            % local execution of the vectorized log_sum_exp_positive
            alpha_terms_max = max(alpha_terms_list, [], 2);
            alpha_terms_diff = alpha_terms_list - repmat(alpha_terms_max, [1, K]);
            alpha_terms_diff_sum_exp_log = log(sum(exp(alpha_terms_diff), 2));
            
            % assignment of alpha matrix element values
            alpha_matrix_log(i_possible,t) = alpha_terms_diff_sum_exp_log + ...
                alpha_terms_max + eta_log_list(i_possible, t);
        end
        
        % log-likelihood calculation
        logL_tot = logL_tot + ...
            log_sum_exp_positive(alpha_matrix_log(:,fluo_lengths{i_tr}));
    end