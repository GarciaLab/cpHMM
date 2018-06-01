function local_em_outputs = local_em_iid (fluo_values, v, noise, K, w, ...
                                          alpha, n_steps_max, eps)

    % Returns the estimated of the emission and noise parameters, assuming
    % that the data is i.i.d. and has no temporal dependence.
    % 
    % INPUTS
    % fluo_values: cell data type with multiple time series fluorescence
    %              values that come from the same AP position
    % v: initial emission values
    % noise: initial gaussian noise
    % K: number of naive states
    % w: memory
    % alpha: length of the MS2 loop in time steps
    % n_steps_max: maximum number of backward-forward iterations
    % eps: relative change in model parameters for which the backward-forward
    %      iterations stop and the algorithm is assumed to converge
    %
    % OUTPUTS
    % Inferred model parameters when:
    % - the relative change in model parameters is smaller than eps, or
    % - the number of backward-forward iterations exceeds n_steps_max
    % 
    % v_logs: logs of the absolute values of the emission parameters
    % v_sings: sings of the inferred emission parameters
    % lambda_log: log of the inferred lambda parameter defined as 1/noise^2

    % ---------------------- Variable assignments ----------------------
    
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
    
    % ---------------------- Accounting variables ----------------------

    % F_{zl} is the number of times the naive state z is encountered in the
    % naive state representation of the compound state l, adjusted by 
    % coefficients that account for the finite MS2 loop size
    F = zeros(K, K^w);

    ms2_coeff = ms2_loading_coeff(alpha,w);

    for i = 1:K^w
        naive = compound_to_naive(i, K, w);
        for k = 1:K
            F(k,i) = sum(ms2_coeff(naive == k));
        end
    end
    
    % ----------------------------- Lists ------------------------------
    
    % list of the first digits (naive states) of each compound state
    digit_first_list = zeros(1, K^w, 'int16');
    
    % list of the number of each naive state encountered in compound states
    naive_count_list_MS2 = zeros(K^w, K);
    
    for i = 1:K^w
        digit_first_list(i) = digit(i, 1, K, w);
        
        naive = compound_to_naive(i, K, w);
        for k = 1:K
            naive_count_list_MS2(i,k) = sum(ms2_coeff(naive == k));
        end
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
    
    % variable used to account for adjustements at time points 1:(w-1)
    count_reduction = zeros([1,w-1]);
    for t = 1:(w-1)
        count_reduction(t) = sum(ms2_coeff((t+1):end));
    end
    
    % calculation of the log of the F variable used in v maximization
    log_F_terms = cell([n_traces, 1]);
    for i_tr = 1:n_traces
        log_F_terms{i_tr} = cell([K, 1]);
        for n = 1:K
            log_F_terms{i_tr}{n} = repmat(log(F(n,:))', 1, fluo_lengths{i_tr});
        end
        
        temp_var = log_F_terms{i_tr}{1};
        for t = 1 : (w-1)
            temp_var(:,t) = log(abs(exp(temp_var(:,t)) - count_reduction(t)));
        end
        
        log_F_terms{i_tr}{1} = temp_var;
    end
    
    % logs and signs of the fluorescence data used in v maximization
    x_term_logs = cell([n_traces, 1]);
    x_term_signs = cell([n_traces, 1]);
    for i_tr = 1:n_traces
        x_term_logs{i_tr} = repmat(log(abs(fluo_values{i_tr})), K^w, 1);
        x_term_signs{i_tr} = repmat(sign(fluo_values{i_tr}), K^w, 1);
    end

    % ------------- Expectation-Maximizaiton iterations ----------------
    for baum_welch = 1:n_steps_max
        
        %%%%%%%%%%%%%%%%%%%%%%%%% EXPECTATION %%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % ------------- lists used in the expectation step -------------
        
        % list of log (X_t - V_t)^2 terms that appear in the emission pdf's
        difference_list = cell([n_traces, 1]);
        eta_log_list = cell([n_traces, 1]);
        for i_tr = 1:n_traces
            difference_list_temp = zeros(K^w, fluo_lengths{i_tr});
            for i = 1:n_unique
                states = naive_count_map{i};
                difference_list_temp(states, :) = ...
                   repmat(difference(fluo_logs{i_tr}, fluo_signs{i_tr}, fluo_lengths{i_tr}, K, w, alpha, ...
                   states(1), v_logs, v_signs, naive_count_list_MS2_log)', [length(states), 1]);
            end
            difference_list{i_tr} = difference_list_temp;
            eta_log_list{i_tr} = 0.5*(lambda_log - log(2*pi)) ...
                -0.5*exp(lambda_log + difference_list{i_tr});
        end
        
        smudge = 0;
        % --------------------- <S_t> calculation ----------------------
        p_s_log = cell([n_traces, 1]);
        for i_tr = 1:n_traces
            p_s_log{i_tr} = eta_log_list{i_tr};
            for t = 2: min([fluo_lengths{i_tr}, (w-1)])
                % possible states at time t
                i_possible = possible_states_list{t};
                i_impossible = setdiff(1:(K^w), i_possible);
                p_s_log{i_tr}(i_impossible,t) = -Inf;
            end
            for t = 1 : fluo_lengths{i_tr}
                smudge = smudge + log_sum_exp_positive(p_s_log{i_tr}(:,t));
                p_s_log{i_tr}(:,t) = p_s_log{i_tr}(:,t) - ...
                    log_sum_exp_positive(p_s_log{i_tr}(:,t));
                
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%% MAXIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%
        
        % ------------------------- lambda_log -------------------------
        lambda_log_old = lambda_log;

        arr = [];
        for i_tr = 1:n_traces
            term_ith = p_s_log{i_tr} + difference_list{i_tr};
            arr = [arr, term_ith(:)'];
        end
        arr_max = max(arr(:));
        lambda_log = log(sum(cell2mat(fluo_lengths))) - (arr_max + log(sum(exp(arr(:)-arr_max))));

        noise_log_old = -0.5 * lambda_log_old;
        noise_log = -0.5 * lambda_log;
        
        noise_change = exp(noise_log).*abs(exp(noise_log_old - noise_log) - 1);
        noise_rel_change = noise_change / exp(noise_log_old);

        % ----------------------------- v ------------------------------
        v_logs_old = v_logs;

        m_sign = ones(K, K);
        m_log = zeros(K, K);

        for m = 1:K
            for n = 1:K
                terms = [];
                for i_tr = 1:n_traces
                    terms_ith = p_s_log{i_tr} + log_F_terms{i_tr}{n} + log_F_terms{i_tr}{m};
                    terms = [terms, terms_ith(:)'];
                end
                terms_max = max(terms(:));
                if (terms_max == -Inf)
                    m_log(m, n) = -Inf;
                else
                    m_log(m, n) = terms_max + log(sum(exp(terms(:)-terms_max)));
                end
            end
        end

        b_sign = ones(1, K);
        b_log = zeros(1, K);
        
        for m = 1:K
            terms_b_log = [];
            terms_b_sign = [];
            for i_tr = 1:n_traces
                terms_b_log_ith = x_term_logs{i_tr} + p_s_log{i_tr} + ...
                    + log_F_terms{i_tr}{m};
                terms_b_log = [terms_b_log, terms_b_log_ith(:)'];
                terms_b_sign = [terms_b_sign, x_term_signs{i_tr}(:)'];
            end
            terms_b_sum = log_sum_exp(terms_b_log, terms_b_sign);
            b_log(m) = terms_b_sum(1);
            b_sign(m) = terms_b_sum(2);
        end

        v_updated = v_log_solve(m_log, m_sign, b_log, b_sign);
        v_logs = v_updated(1,:);
        v_signs = v_updated(2,:);

        v_norm_change = abs(norm(exp(v_logs_old)) - norm(exp(v_logs)));
        v_norm_rel_change = v_norm_change / norm(exp(v_logs_old));
        
        % ------------------- convergence criterion --------------------
        if (max([noise_rel_change, v_norm_rel_change]) < eps)
            break
        end
    end
    
    % -------------- collection of outputs into a cell -----------------
    local_em_outputs = struct('v_logs', v_logs, 'v_signs', v_signs, ...
                              'lambda_log', lambda_log);