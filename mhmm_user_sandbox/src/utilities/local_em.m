function local_em_outputs = local_em (fluo_values, ...
    v, noise, pi0_log, A_log, K, w, n_steps_max, eps)

    % The function returns the maximum likelihood estimates of the model
    % parameters, given multiple data sequences and initial conditions
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
    % n_steps_max: maximum number of backward-forward iterations
    % eps: relative change in model parameters for which the backward-forward
    %      iterations stop and the algorithm is assumed to converge
    %
    % OUTPUTS
    % Inferred model parameters when:
    % - the relative change in model parameters is smaller than eps,
    % - or the number of backward-forward iterations exceeds n_steps_max
    % 
    % pi0_log: log of the inferred initial state cmf
    % A_log: log of the inferred transition matrix
    % v_logs: logs of the absolute values of the emission parameters
    % v_sings: sings of the inferred emission parameters
    % lambda_log: log of the inferred lambda parameter defined as 1/noise^2
    % baum_welch: the number of forward-backward iterations done before
    %             convergence
    % log_likelihoods_total(1:baum_welch): log-likelihood of observing the
    %                                fluoresence data evaluated using the
    %                                inferred model parameters
    % p_s_log: expected log-probabilities of each compound state at all
    %          time points for each trace
    % p_ss_log: expected log-probabilities of joint compound states at all
    %           time points for each trace
    
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
    
    % total number of time points in the data set
    fluo_length_total = sum(cell2mat(fluo_lengths));
    
    % ---------------------- Accounting variables ----------------------
    
    % C_{zm} = 1 iff the first digit of compound state m is equal to z
    C = zeros(K, K^w);
    for i = 1:K
        columns = ((i-1)*K^(w-1) : ((i*K^(w-1))-1)) + 1;
        C(i, columns(:)) = 1;
    end

    % log of the accounting variable C
    logC = log(C);

    % F_{zm} is the number of times the naive state z is encountered in the
    % naive state representation of the compound state m
    F = zeros(K, K^w);
    for i = 1:K^w
        F(:,i) = naive_count(i, K, w);
    end
    
    % ----------------------------- Lists ------------------------------
    
    % list of states that can transition to each of 1:K^w compounds states
    allowed_to_list = zeros(K^w, K, 'int16');
    
    % list of states allowed to go from each of 1:K^w compound states
    allowed_from_list = zeros(K^w, K, 'int16');
    
    % list of the first digits (naive states) of each compound state
    digit_first_list = zeros(1, K^w, 'int16');
    
    % list of the number of each naive state encountered in compound states
    naive_count_list = zeros(K^w, K, 'int16');
    
    for i = 1:K^w
        allowed_to_list(i,:) = allowed_to(i, K, w);
        allowed_from_list(i,:) = allowed_from(i, K, w);
        digit_first_list(i) = digit(i, 1, K, w);
        naive_count_list(i, :) = naive_count (i, K, w);
    end          
    
    % log of the naive count list
    naive_count_log_list = log(double(naive_count_list));
    
    % list of compound states that have a unique sequence of the number of
    % naive states contained in the base-K representation
    unique_naive_list = unique(naive_count_list,'rows');
    [n_unique, ~] = size(unique_naive_list);    
    
    % a list that maps the unique naive count combinations with compound
    % states
    naive_count_map = cell(n_unique, 1);   
    for i = 1:n_unique
        naive = unique_naive_list(i,:);
        naive_count_map{i} = find(ismember(naive_count_list,naive,'rows'));        
    end
    
    % list of possible compound states at all times
    time_max = max(cell2mat(fluo_lengths));
    possible_states_list = cell(time_max, 1);    
    for t = 1:time_max
        possible_states_list{t} = possible_states (t, K, w);        
    end
    
    
    % calculation of p_ss_log indices in K^w x K^w x T representation
    ind_3d = cell([n_traces, 1]);
    for i_tr = 1:n_traces
        % list of additional 1d indices to account for multiple time
        % layers in the 3d representation of the p_ss_log matrix
        ind_addition_mat = repmat((1:(fluo_lengths{i_tr}-1))*K^(2*w), [K^(w+1), 1]);
        ind_addition_list = ind_addition_mat(:);
        
        % list of 1d indices of the positions of pairs in the 2d slice of
        % p_ss_log which correspond to valid col -> row transitions
        ind_positions_2d = sub2ind([K^w, K^w], allowed_from_list(:), repmat(1:K^w, [1, K])');
        
        % replication of the 1d pair indices along the time axis
        ind_positions_3d = repmat(ind_positions_2d, [1, (fluo_lengths{i_tr}-1)]);
        
        % combination of the pair indices with the additional 1d indices
        ind_3d{i_tr} = ind_addition_list + ind_positions_3d(:);
    end
    
    % calculation of p_ss_log indices for all traces used in A_log 
    % maximization
    A_log_maximization_ind = cell([n_traces, 1]);
    for i_tr = 1:n_traces
        A_log_maximization_ind{i_tr} = cell(K);
        
        ind_addition_mat = repmat((1:(fluo_lengths{i_tr}-1))*K^(2*w), [K^(w-1), 1]);
        ind_addition_list = ind_addition_mat(:);
        
        for m = 1:K
            for n = 1:K
                d_n_list = find(digit_first_list == n);
                d_m_list = find(digit_first_list == m);
                allowed_from_n = allowed_from_list(d_n_list, :);
                column_ind = ismember(allowed_from_n(1,:), d_m_list) == 1;
                ind_2d_A_log = sub2ind([K^w, K^w], allowed_from_n(:, column_ind), d_n_list');
                ind_3d_A_log = repmat(ind_2d_A_log, [1, fluo_lengths{i_tr}-1]);

                A_log_maximization_ind{i_tr}{m,n} = find(ismember(ind_3d{i_tr}, ...
                    ind_addition_list + ind_3d_A_log(:)));
            end
        end
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
            temp_var(:,t) = log(abs(exp(temp_var(:,t)) - (w-t)));
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

    
    % ------------------------ Pre-allocations -------------------------
    % log likelihoods of observing the fluorescence sequence at each EM
    % forward-backward iteration    
    log_likelihoods = cell([n_traces, 1]);
    for i_tr = 1:n_traces
        log_likelihoods{i_tr} = -Inf(1, n_steps_max);
    end
    log_likelihoods_total = -Inf(1, n_steps_max);
    
    % pre-allocation of alpha and beta coefficient matrices
    alpha_matrix_log = cell([n_traces, 1]);
    beta_matrix_log = cell([n_traces, 1]);
    for i_tr = 1:n_traces
        alpha_matrix_log{i_tr} = zeros(K^w,fluo_lengths{i_tr});                
        beta_matrix_log{i_tr} = zeros(K^w,fluo_lengths{i_tr});
    end

    % pre-allocation of the p_ss_log array used in the expectation step
    p_ss_log = cell([n_traces, 1]);
    for i_tr = 1:n_traces
        p_ss_log{i_tr} = -Inf([1, K^(w+1) * fluo_lengths{i_tr}]);
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
                   repmat(difference(fluo_logs{i_tr}, fluo_signs{i_tr}, fluo_lengths{i_tr}, K, w, ...
                   states(1), v_logs, v_signs, naive_count_log_list)', [length(states), 1]);
            end
            difference_list{i_tr} = difference_list_temp;
            eta_log_list{i_tr} = 0.5*(lambda_log - log(2*pi)) ...
                -0.5*exp(lambda_log + difference_list{i_tr});
        end
        
        % row and column subscripts used for indexing A_log elements for
        % alpha matrix calculation
        A_log_alpha_rowSubs = digit_first_list(repmat(1:K^w, [1, K]));
        A_log_alpha_colSubs = digit_first_list(allowed_to_list(:));
        
        % A_log element indexing in a 2d slice used for alpha matrix
        % calculation
        A_log_alpha_subs_1d = A_log(sub2ind([K, K], A_log_alpha_rowSubs, A_log_alpha_colSubs));
        alpha_A_log_list = reshape(A_log_alpha_subs_1d, [K^w K]);        
        
        % row and column subscripts used for indexing A_log elements for
        % alpha matrix calculation
        A_log_beta_rowSubs = digit_first_list(allowed_from_list(:));
        A_log_beta_colSubs = digit_first_list(repmat(1:K^w, [1, K]));
        
        % A_log element indexing in a 2d slice used for beta matrix
        % calculation
        A_log_beta_subs_1d = A_log(sub2ind([K, K], A_log_beta_rowSubs, A_log_beta_colSubs));
        beta_A_log_list = reshape(A_log_beta_subs_1d, [K^w K]);                       
        
        
        % ------------ alpha coefficient matrix calculation ------------ 
                
        for i_tr = 1:n_traces
            
            % initializes the alpha matrix as all zeros (logs all infinities)            
            alpha_matrix_log{i_tr}(:,:) = -Inf;
            
            % calculates the alpha matrix elements for t = 1
            for i = possible_states_list{1}
                alpha_matrix_log{i_tr}(i, 1) = eta_log_list{i_tr}(i, 1) + ...
                    pi0_log(digit_first_list(i));
            end
        
            % calculates the alpha matrix elements for t > 1
            for t = 2:fluo_lengths{i_tr}
                % possible states at time t
                i_possible = possible_states_list{t};

                % list of terms that are added to find the alpha matrix
                % elements
                alpha_terms_list = alpha_A_log_list(i_possible, :) + ...
                    reshape(alpha_matrix_log{i_tr}(allowed_to_list(i_possible,:),t-1), [], K);             

                % local execution of the vectorized log_sum_exp_positive
                alpha_terms_max = max(alpha_terms_list, [], 2); 
                alpha_terms_diff = alpha_terms_list - repmat(alpha_terms_max, [1, K]);
                alpha_terms_diff_sum_exp_log = log(sum(exp(alpha_terms_diff), 2));
                
                % fixed the 'nan' values for bases when some elements of
                % alpha_terms_max are '-inf'
                ind_inf = isnan(alpha_terms_diff_sum_exp_log);
                alpha_terms_diff_sum_exp_log(ind_inf) = 0;

                % assignment of alpha matrix element values
                alpha_matrix_log{i_tr}(i_possible,t) = alpha_terms_diff_sum_exp_log + ...
                    alpha_terms_max + eta_log_list{i_tr}(i_possible, t);
            end
        end

        % ------------ beta coefficient matrix calculation -------------
        
        for i_tr = 1:n_traces
            
            % initializes the alpha matrix as all zeros (logs all infinities)
            beta_matrix_log{i_tr}(:,:) = -Inf;       

            % assigns 1 to beta matrix elements at t=time (or 0's to log's)
            beta_matrix_log{i_tr}(:, fluo_lengths{i_tr}) = 0;

            % calculates the alpha matrix elements for t < time
            for t = (fluo_lengths{i_tr}-1: -1: 1)
                % possible states at time t
                i_possible = possible_states_list{t};

                % list of terms that are added to find the beta matrix elements
                beta_terms_list = beta_A_log_list(i_possible, :) + ...
                    reshape(eta_log_list{i_tr}(allowed_from_list(i_possible,:), t+1), [], K) + ...
                    reshape(beta_matrix_log{i_tr}(allowed_from_list(i_possible,:), t+1), [], K);

                % local execution of the vectorized log_sum_exp_positive
                beta_terms_max = max(beta_terms_list, [], 2);
                beta_terms_diff = beta_terms_list - repmat(beta_terms_max, [1, K]);
                beta_terms_diff_sum_exp_log = log(sum(exp(beta_terms_diff), 2));
                
                % fixed the 'nan' values for bases when some elements of
                % beta_terms_max are '-inf'
                ind_inf = isnan(beta_terms_diff_sum_exp_log);
                beta_terms_diff_sum_exp_log(ind_inf) = 0;

                % assignment of beta matrix element values
                beta_matrix_log{i_tr}(i_possible,t) = beta_terms_max + ...
                    beta_terms_diff_sum_exp_log;
            end
        end
        
        
        % --------------- log-likelihood calcuation --------------------
        log_likelihoods_total(baum_welch) = 0;
        for i_tr = 1:n_traces
            log_likelihoods{i_tr}(baum_welch) = ...
                log_sum_exp_positive(alpha_matrix_log{i_tr}(:,fluo_lengths{i_tr}));
            log_likelihoods_total(baum_welch) = log_likelihoods_total(baum_welch) + ...
                log_likelihoods{i_tr}(baum_welch);
        end
        
        % --------------------- <S_t> calculation ----------------------
        p_s_log = cell([n_traces, 1]);
        for i_tr = 1:n_traces
            p_s_log{i_tr} = alpha_matrix_log{i_tr} + beta_matrix_log{i_tr} ...
                - log_likelihoods{i_tr}(baum_welch);
        end
        
        % ------------------- <S_t, S_{t-1}> calculation ---------------

        for i_tr = 1:n_traces
            
            % replication of the alpha matrix K times along one of the 2d axes
            alpha_matrix_log_minus = alpha_matrix_log{i_tr}(:,1:(fluo_lengths{i_tr}-1));
            alpha_matrix_log_minus_rep = repmat(alpha_matrix_log_minus, [K, 1]);                        

            % 2d indexing of eta_log and beta_log matrices
            % note: the 2d slice components are accounted for through the
            %       1d representation of the matrix of possible transitions.
            %       Also, note that the 1d indexing of the p_ss_log matrix
            %       is done by rows, which is used in calculating the terms
            eta_log_list_minus = eta_log_list{i_tr}(allowed_from_list(:), 2:fluo_lengths{i_tr});
            beta_matrix_log_minus = beta_matrix_log{i_tr}(allowed_from_list(:), 2:fluo_lengths{i_tr});        

            % row and column subscripts used for indexing A_log elements
            A_log_rowSubs = digit_first_list(allowed_from_list(:));
            A_log_colSubs = digit_first_list(repmat(1:K^w, [1, K]));

            % A_log element indexing in a 2d slice
            A_log_sub_single = A_log(sub2ind([K, K], A_log_rowSubs, A_log_colSubs));

            % A_log element indexing in a 3d slice along the time axis
            A_log_colSubs_multi = repmat(A_log_sub_single, [1, fluo_lengths{i_tr}-1]);                

            % calcuation of p_ss_log elements
            p_ss_log{i_tr} = alpha_matrix_log_minus_rep(:) + ...
                eta_log_list_minus(:) + beta_matrix_log_minus(:) + ...
                A_log_colSubs_multi(:) - log_likelihoods{i_tr}(baum_welch);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%% MAXIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%

        % -------------------------- pi0_log ---------------------------
        pi0_old = exp(pi0_log);
        for m = 1:K
            pi0_terms = [];
            for i_tr = 1:n_traces
                pi0_terms = [pi0_terms, logC(m, :) + p_s_log{i_tr}(:,1)'];
            end
            pi0_log(m) = log_sum_exp_positive(pi0_terms) - log(n_traces);
        end
        
        pi0_norm_rel_change = ...
            abs(norm(pi0_old) - norm(exp(pi0_log)))/norm(pi0_old);

        % --------------------------- A_log ----------------------------
        A_old = exp(A_log);
        for m = 1:K
            for n = 1:K
                arr = [];
                for i_tr = 1:n_traces
                    p_ss_log_ith = p_ss_log{i_tr}(A_log_maximization_ind{i_tr}{m,n});
                    arr = [arr, p_ss_log_ith(:)'];
                end
                arr_max = max(arr(:));
                if (arr_max == -Inf)
                    A_log(m,n) = -Inf;
                else
                    A_log(m,n) = arr_max + log(sum(exp(arr(:)-arr_max)));
                end
            end
        end

        for n = 1:K
            arr = A_log(:,n);
            arr_max = max(arr(:));
            A_log(:,n) = A_log(:,n) - (arr_max + log(sum(exp(arr(:)-arr_max))));
        end

        A_norm_rel_change = ...
            abs(norm(A_old) - norm(exp(A_log)))/norm(A_old);
        
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

        v_updated = v_log_solve_test(m_log, m_sign, b_log, b_sign);
        v_logs = v_updated(1,:);
        v_signs = v_updated(2,:);

        v_norm_change = abs(norm(exp(v_logs_old)) - norm(exp(v_logs)));
        v_norm_rel_change = v_norm_change / norm(exp(v_logs_old));
        
        % --------------- change in logL per time step -----------------
        logL_norm_change = 0;
        if (baum_welch > 1)
            logL_norm_change = log_likelihoods_total(baum_welch) - ...
                log_likelihoods_total(baum_welch-1);
            logL_norm_change = abs(logL_norm_change) / fluo_length_total;
        end

        % ------------------- convergence criterion --------------------
        if (max([pi0_norm_rel_change, A_norm_rel_change, noise_rel_change, ...
                v_norm_rel_change, logL_norm_change]) < eps)
            log_likelihoods_total = log_likelihoods_total(1:baum_welch);
            break
        end
    end
    
    % -------------- collection of outputs into a cell -----------------
    local_em_outputs = struct('pi0_log', pi0_log, 'A_log', A_log, ...
        'v_logs', v_logs, 'v_signs', v_signs, 'lambda_log', lambda_log, ...
        'logL_bw', log_likelihoods_total, 'logL', max(log_likelihoods_total), ...
        'max_bw', baum_welch);
    local_em_outputs.p_s_log = p_s_log;
    local_em_outputs.p_ss_log = p_ss_log;