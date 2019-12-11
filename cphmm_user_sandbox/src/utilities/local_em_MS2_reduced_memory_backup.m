function local_em_outputs = local_em_MS2_reduced_memory_backup (fluo_values, ...
              v, noise, pi0_log, A_log, K, w, kappa, n_steps_max, eps,...
              use_backup, backup_file_path)

    % Returns the estimates of the model parameters, given multiple data 
    % sequences and initial conditions.
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
    % n_steps_max: maximum number of backward-forward iterations
    % eps: minimum relative change in model parameters need to stop the
    %      backward-forward algorithm iterations
    %
    % OUTPUTS
    % Inferred model parameters when:
    % - the relative change in model parameters is smaller than eps,
    % - or the number of backward-forward iterations exceeds n_steps_max
    % 
    % pi0_log: log of the inferred initial state cmf
    % A_log: log of the inferred transition matrix
    % v_logs: logs of the absolute values of the emission parameters
    % v_signs: signs of the inferred emission parameters
    % lambda_log: log of the inferred lambda parameter defined as 1/noise^2
    % n_iter: number of forward-backward iterations before convergence
    % logL: log-likelihood at the end of iterations
    % runtime: the total runtime of the EM algorithm
    % use_backup: binary flag. If 1 save individual em results to file
    % backup_file_path: path to save local EM results (useful in cases when
    % full inference fails)
    
    %%%%%%%%%%%%%%%%%%%%%%% Variable assignments %%%%%%%%%%%%%%%%%%%%%%%
    
    % begin the timer
    tic
    
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

    % F_{zl} is the number of times the naive state z is encountered in the
    % naive state representation of the compound state l, adjusted by 
    % coefficients that account for the finite MS2 loop size
    F = zeros(K, K^w);

    % mean fractions of MS2 loops transcribed at each elongation step
    ms2_coeff = ms2_loading_coeff(kappa,w);

    for i = 1:K^w
        naive = compound_to_naive(i, K, w);
        for k = 1:K
            F(k,i) = sum(ms2_coeff(naive == k));
        end
    end
    
    % variable used to account for adjustements at time points 1:(w-1)
    count_reduction = zeros([1,w-1]);
    for t = 1:(w-1)
        count_reduction(t) = sum(ms2_coeff((t+1):end));
    end

    
    % ----------------------------- Lists ------------------------------
    
    % list of states that can transition to each of 1:K^w compounds states
    allowed_to_list = zeros(K^w, K, 'int16');
    
    % list of states allowed to go from each of 1:K^w compound states
    allowed_from_list = zeros(K^w, K, 'int16');
    
    % list of the first digits (naive states) of each compound state
    digit_first_list = zeros(1, K^w, 'int16');
    
    % list of the number of each naive state encountered in compound states
    naive_count_list_MS2 = zeros(K^w, K);
    
    for i = 1:K^w
        allowed_to_list(i,:) = allowed_to(i, K, w);
        allowed_from_list(i,:) = allowed_from(i, K, w);
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
        naive_count_map{i} = find(ismember(naive_count_list_MS2, ...
                                  naive, 'rows'));
    end
    
    % list of possible compound states at all times
    time_max = max(cell2mat(fluo_lengths));
    possible_states_list = cell(time_max, 1);
    for t = 1:time_max
        possible_states_list{t} = possible_states (t, K, w);        
    end

    % calculation of p_ss_log indices in K^w x K^w x T representation
%     ind_positions_2d_mat = sub2ind([K^w, K^w], allowed_from_list(:), ...
%                                    repmat(1:K^w, [1, K])');
%     ind_positions_2d = ind_positions_2d_mat(:);
    
    ind_positions_2d = double(allowed_from_list(:)) + ...
                              K^w * (repmat(1:K^w, [1, K])'-1);
    
    ind_2d = cell([K, K]);
    for m = 1:K
        for n = 1:K
            d_n_list = find(digit_first_list == n);
            d_m_list = find(digit_first_list == m);
            allowed_from_n = allowed_from_list(d_n_list, :);
            column_ind = ismember(allowed_from_n(1,:), d_m_list) == 1;
%             ind_2d_A_log_mat = sub2ind([K^w, K^w], ...
%                                allowed_from_n(:, column_ind), d_n_list');
            ind_2d_A_log_mat = double(allowed_from_n(:, column_ind)) + ...
                               K^w * (d_n_list'-1);
            ind_2d{m,n} = find(ismember(ind_positions_2d, ...
                                           ind_2d_A_log_mat(:)));
        end
    end
    
    
    % ------------------------ Pre-allocations -------------------------
    
    % log likelihoods of observing the fluorescence sequence at each EM
    % forward-backward iteration    
    log_likelihoods = cell([n_traces, 1]);
    for i_tr = 1:n_traces
        log_likelihoods{i_tr} = -Inf(1, n_steps_max);
    end
    logL_tot = -Inf(1, n_steps_max);
    
    
    %%%%%%%%%%%%%%%% Expectation-Maximizaiton algorithm %%%%%%%%%%%%%%%%
    
    for baum_welch = 1:n_steps_max
        
        % row and column subscripts used for indexing A_log elements for
        % alpha matrix calculation
        A_log_alpha_rowSubs = digit_first_list(repmat(1:K^w, [1, K]));
        A_log_alpha_colSubs = digit_first_list(allowed_to_list(:));
        
        % A_log element indexing in a 2d slice used for alpha matrix
        % calculation
        A_log_alpha_subs_1d = A_log(sub2ind([K,K], A_log_alpha_rowSubs, ...
                                      A_log_alpha_colSubs));
        alpha_A_log_list = reshape(A_log_alpha_subs_1d, [K^w K]);
        
        % row and column subscripts used for indexing A_log elements for
        % alpha matrix calculation
        A_log_beta_rowSubs = digit_first_list(allowed_from_list(:));
        A_log_beta_colSubs = digit_first_list(repmat(1:K^w, [1, K]));
        
        % A_log element indexing in a 2d slice used for beta matrix
        % calculation
        A_log_beta_subs_1d = A_log(sub2ind([K, K], A_log_beta_rowSubs, ...
                                        A_log_beta_colSubs));
        beta_A_log_list = reshape(A_log_beta_subs_1d, [K^w K]);
        
        % total log-likelihood
        logL_tot(baum_welch) = 0;
        
        % Model parameter expectation terms that are updated when new
        % traces are taken into account
        pi0_terms = -Inf([1,K]);
        A_terms = -Inf([K,K]);
        lambda_terms = -Inf;
        v_M_terms = -Inf([K,K]);
        v_b_terms_log = -Inf([1,K]);
        v_b_terms_sign = ones([1,K]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%% Expectation %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for i_tr = 1:n_traces
            
            % calculation of the log of the F variable used in the
            % maximization of v
            log_F_terms = cell([K, 1]);
            for n = 1:K
                log_F_terms{n} = repmat(log(F(n,:))', 1, ...
                                         fluo_lengths{i_tr});
            end
            temp_var = log_F_terms{1};
            for t = 1 : (w-1)
                temp_var(:,t) = log(abs(exp(temp_var(:,t)) - ...
                                     count_reduction(t)));
            end
            log_F_terms{1} = temp_var;
            
            % calculation of p_ss_log indices for all traces used in A_log
            % maximization
            A_log_maximization_ind = cell(K);
            for m = 1:K
                for n = 1:K
                    ind_2d_rep = repmat(ind_2d{m,n}, [1, fluo_lengths{i_tr}-1]);
                    ind_addition_ls = repmat((1:(fluo_lengths{i_tr}-1))*K^(w+1), ...
                        [length(ind_2d{m,n}), 1]) - K^(w+1);
                    A_log_maximization_ind{m,n} = ind_2d_rep(:) + ind_addition_ls(:);
                end
            end
            
            % logs and signs of the fluorescence data used in v maximization
            x_term_logs = repmat(log(abs(fluo_values{i_tr})), K^w, 1);
            x_term_signs = repmat(sign(fluo_values{i_tr}), K^w, 1);
            
            
            % ------------- lists used in the expectation step ---------
            
            % list of log (X_t - V_t)^2 terms that appear in the emission pdf's
            difference_list_temp = zeros(K^w, fluo_lengths{i_tr});
            if K ~= length(v_signs)
                error('inconsitent state id and vec dims')
            end
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
            
            
            % ------------ alpha coefficient matrix calculation --------
            
            % pre-allocation of alpha and beta coefficient matrices
            alpha_matrix_log = zeros(K^w,fluo_lengths{i_tr});

            % initializes the alpha matrix as all zeros (logs as Infs)            
            alpha_matrix_log(:,:) = -Inf;
            
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
                
                % fixed the 'nan' values for bases when some elements of
                % alpha_terms_max are '-inf'
                ind_inf = isnan(alpha_terms_diff_sum_exp_log);
                alpha_terms_diff_sum_exp_log(ind_inf) = 0;

                % assignment of alpha matrix element values
                alpha_matrix_log(i_possible,t) = alpha_terms_diff_sum_exp_log + ...
                    alpha_terms_max + eta_log_list(i_possible, t);
            end
            
            % ------------ beta coefficient matrix calculation ---------
            
            % pre-allocation of alpha and beta coefficient matrices
            beta_matrix_log = zeros(K^w,fluo_lengths{i_tr});
    
            % initializes the alpha matrix as all zeros (logs all infinities)
            beta_matrix_log(:,:) = -Inf;

            % assigns 1 to beta matrix elements at t=time (or 0's to log's)
            beta_matrix_log(:, fluo_lengths{i_tr}) = 0;

            % calculates the alpha matrix elements for t < time
            for t = (fluo_lengths{i_tr}-1: -1: 1)
                % possible states at time t
                i_possible = possible_states_list{t};

                % list of terms that are added to find the beta matrix elements
                beta_terms_list = beta_A_log_list(i_possible, :) + ...
                    reshape(eta_log_list(allowed_from_list(i_possible,:), t+1), [], K) + ...
                    reshape(beta_matrix_log(allowed_from_list(i_possible,:), t+1), [], K);

                % local execution of the vectorized log_sum_exp_positive
                beta_terms_max = max(beta_terms_list, [], 2);
                beta_terms_diff = beta_terms_list - repmat(beta_terms_max, [1, K]);
                beta_terms_diff_sum_exp_log = log(sum(exp(beta_terms_diff), 2));
                
                % fixed the 'nan' values for bases when some elements of
                % beta_terms_max are '-inf'
                ind_inf = isnan(beta_terms_diff_sum_exp_log);
                beta_terms_diff_sum_exp_log(ind_inf) = 0;

                % assignment of beta matrix element values
                beta_matrix_log(i_possible,t) = beta_terms_max + ...
                    beta_terms_diff_sum_exp_log;
            end
            
            
            % --------------- log-likelihood calcuation --------------------
            log_likelihoods{i_tr}(baum_welch) = ...
                log_sum_exp_positive(alpha_matrix_log(:,fluo_lengths{i_tr}));
            logL_tot(baum_welch) = logL_tot(baum_welch) + ...
                log_likelihoods{i_tr}(baum_welch);
            
            
            % --------------------- <S_t> calculation ----------------------
            p_s_log = alpha_matrix_log + beta_matrix_log ...
                    - log_likelihoods{i_tr}(baum_welch);
            
            % ------------------- <S_t, S_{t-1}> calculation ---------------
            
            % replication of the alpha matrix K times along one of the 2d axes
            alpha_matrix_log_minus = alpha_matrix_log(:,1:(fluo_lengths{i_tr}-1));
            alpha_matrix_log_minus_rep = repmat(alpha_matrix_log_minus, [K, 1]);

            % 2d indexing of eta_log and beta_log matrices
            % note: the 2d slice components are accounted for through the
            %       1d representation of the matrix of possible transitions.
            %       Also, note that the 1d indexing of the p_ss_log matrix
            %       is done by rows, which is used in calculating the terms
            eta_log_list_minus = eta_log_list(allowed_from_list(:), 2:fluo_lengths{i_tr});
            beta_matrix_log_minus = beta_matrix_log(allowed_from_list(:), 2:fluo_lengths{i_tr});

            % row and column subscripts used for indexing A_log elements
            A_log_rowSubs = digit_first_list(allowed_from_list(:));
            A_log_colSubs = digit_first_list(repmat(1:K^w, [1, K]));

            % A_log element indexing in a 2d slice
            A_log_sub_single = A_log(sub2ind([K, K], A_log_rowSubs, A_log_colSubs));

            % A_log element indexing in a 3d slice along the time axis
            A_log_colSubs_multi = repmat(A_log_sub_single, [1, fluo_lengths{i_tr}-1]);                

            % calcuation of p_ss_log elements
            p_ss_log = alpha_matrix_log_minus_rep(:) + ...
                eta_log_list_minus(:) + beta_matrix_log_minus(:) + ...
                A_log_colSubs_multi(:) - log_likelihoods{i_tr}(baum_welch);
            
            % ---- Updates of the expectation terms for each trace -----

            % pi0_log terms
            for m = 1:K
                pi0_terms(m) = log_sum_exp_positive([pi0_terms(m), ...
                    logC(m, :) + p_s_log(:,1)']);
            end

            % A_log terms
            for m = 1:K
                for n = 1:K
                    p_ss_log_ith = p_ss_log(A_log_maximization_ind{m,n});
                    A_terms(m,n) = log_sum_exp_positive([A_terms(m,n), ...
                        p_ss_log_ith(:)']);
                end
            end

            % lambda_log terms
            term_ith = p_s_log + difference_list;
            lambda_terms = log_sum_exp_positive([lambda_terms, term_ith(:)']);
            
            % v terms
            for m = 1:K
                for n = 1:K
                    terms_ith = p_s_log + log_F_terms{n} + log_F_terms{m};
                    v_M_terms(m,n) = log_sum_exp_positive([v_M_terms(m,n), terms_ith(:)']);
                end
            end

            for m = 1:K
                terms_b_log_ith = x_term_logs + p_s_log + ...
                    + log_F_terms{m};
                tmp = log_sum_exp([v_b_terms_log(m), terms_b_log_ith(:)'], ...
                    [v_b_terms_sign(m), x_term_signs(:)']);
                v_b_terms_log(m) = tmp(1);
                v_b_terms_sign(m) = tmp(2);
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%% MAXIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%

        % -------------------------- pi0_log ---------------------------
        pi0_old = exp(pi0_log);
        pi0_log = pi0_terms - log(n_traces);
        
        pi0_norm_rel_change = ...
            abs(norm(pi0_old) - norm(exp(pi0_log)))/norm(pi0_old);
        
        % --------------------------- A_log ----------------------------
        A_old = exp(A_log);
        A_log = A_terms;

        for n = 1:K
            arr = A_log(:,n);
            arr_max = max(arr(:));
            A_log(:,n) = A_log(:,n) - (arr_max + log(sum(exp(arr(:)-arr_max))));
        end

        A_norm_rel_change = ...
            abs(norm(A_old) - norm(exp(A_log)))/norm(A_old);
        
        % ------------------------- lambda_log -------------------------
        lambda_log_old = lambda_log;
        
        lambda_log = log(sum(cell2mat(fluo_lengths))) - lambda_terms;

        noise_log_old = -0.5 * lambda_log_old;
        noise_log = -0.5 * lambda_log;
        
        noise_change = exp(noise_log).*abs(exp(noise_log_old - noise_log) - 1);
        noise_rel_change = noise_change / exp(noise_log_old);
        
        % ----------------------------- v ------------------------------
        v_logs_old = v_logs;
        
        m_sign = ones(K, K);
        m_log = v_M_terms;

        b_sign = v_b_terms_sign;
        b_log = v_b_terms_log;

        v_updated = v_log_solve(m_log, m_sign, b_log, b_sign);
        v_logs = v_updated(1,:);
        v_signs = v_updated(2,:);

        v_norm_change = abs(norm(exp(v_logs_old)) - norm(exp(v_logs)));
        v_norm_rel_change = v_norm_change / norm(exp(v_logs_old));
        
        % --------------- change in logL per time step -----------------
        logL_norm_change = 0;
        if (baum_welch > 1)
            logL_norm_change = logL_tot(baum_welch) - ...
                logL_tot(baum_welch-1);
            logL_norm_change = abs(logL_norm_change) / fluo_length_total;
        end

        
        % ------------------- convergence criterion --------------------
        if (max([pi0_norm_rel_change, A_norm_rel_change, noise_rel_change, ...
                v_norm_rel_change, logL_norm_change]) < eps)
            logL_tot = logL_tot(1:baum_welch);
            break
        end
    end
    
    % end the timer
    runtime = toc;
    
    % -------------- collection of outputs into a cell -----------------
    local_em_outputs = struct('pi0_log', pi0_log, 'A_log', A_log, ...
        'v_logs', v_logs, 'v_signs', v_signs, 'lambda_log', lambda_log, ...
        'logL_bw', logL_tot, 'logL', max(logL_tot), ...
        'n_iter', baum_welch, 'runtime', runtime);
    if use_backup
        save(backup_file_path, 'local_em_outputs')
    end