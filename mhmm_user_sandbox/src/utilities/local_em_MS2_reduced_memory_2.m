function local_em_outputs = local_em_MS2_reduced_memory_2 (fluo_values, ...
              v, noise, pi0_log, A_log, K, w, kappa, n_steps_max, eps)

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
    
    %%%%%%%%%%%%%%%%%%%%%%% Variable assignments %%%%%%%%%%%%%%%%%%%%%%%
    
    % begin the timer
    tic
    
    % number of traces
    n_traces = length(fluo_values);
    
    % model parameters
    lambda = 1/noise^2;
    A = exp(A_log);
    pi0 = exp(pi0_log);
    
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
    ind_positions_2d = double(allowed_from_list(:)) + ...
                              K^w * (repmat(1:K^w, [1, K])'-1);
    
    ind_2d = cell([K, K]);
    for m = 1:K
        for n = 1:K
            d_n_list = find(digit_first_list == n);
            d_m_list = find(digit_first_list == m);
            allowed_from_n = allowed_from_list(d_n_list, :);
            column_ind = ismember(allowed_from_n(1,:), d_m_list) == 1;
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
        
        % row and column subscripts used for indexing A elements for
        % alpha matrix calculation
        A_alpha_rowSubs = digit_first_list(repmat(1:K^w, [1, K]));
        A_alpha_colSubs = digit_first_list(allowed_to_list(:));
        
        % A element indexing in a 2d slice used for alpha matrix
        % calculation
        A_alpha_subs_1d = A(sub2ind([K,K], A_alpha_rowSubs, ...
                                      A_alpha_colSubs));
        alpha_A_list = reshape(A_alpha_subs_1d, [K^w K]);
        
        % row and column subscripts used for indexing A elements for
        % alpha matrix calculation
        A_beta_rowSubs = digit_first_list(allowed_from_list(:));
        A_beta_colSubs = digit_first_list(repmat(1:K^w, [1, K]));
        
        % A_log element indexing in a 2d slice used for beta matrix
        % calculation
        A_beta_subs_1d = A(sub2ind([K, K], A_beta_rowSubs, ...
                                        A_beta_colSubs));
        beta_A_list = reshape(A_beta_subs_1d, [K^w K]);
        
        % total log-likelihood
        logL_tot(baum_welch) = 0;
        
        % Model parameter expectation terms that are updated when new
        % traces are taken into account
        pi0_terms = zeros([K,1]);
        A_terms = zeros([K,K]);
        lambda_terms = 0;
        v_M_terms = zeros([K,K]);
        v_b_terms = zeros([1,K]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%% Expectation %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for i_tr = 1:n_traces
            
            % calculation of the log of the F variable used in the
            % maximization of v
            F_terms = cell([K, 1]);
            for n = 1:K
                F_terms{n} = repmat(F(n,:)', 1, fluo_lengths{i_tr});
            end
            temp_var = F_terms{1};
            for t = 1 : (w-1)
                temp_var(:,t) = abs(temp_var(:,t) - count_reduction(t));
            end
            F_terms{1} = temp_var;
            
            % calculation of p_ss indices for all traces used in A
            % maximization
            A_maximization_ind = cell(K);
            ind_addition_ls = repmat((1:(fluo_lengths{i_tr}-1))*K^(w+1), ...
                        [length(ind_2d{1,1}), 1]) - K^(w+1);
            for m = 1:K
                for n = 1:K
                    ind_2d_rep = repmat(ind_2d{m,n}, [1, fluo_lengths{i_tr}-1]);
                    A_maximization_ind{m,n} = ind_2d_rep(:) + ind_addition_ls(:);
                end
            end
            
            % ------------- lists used in the expectation step ---------
            
            % list of (X_t - V_t)^2 terms that appear in the emission pdf's
            difference_sq_list = zeros(K^w, fluo_lengths{i_tr});
            for i = 1:n_unique
                states = naive_count_map{i};
                
%                 v_multi = v.*naive_count_list_MS2(states(1),:);
%                 v_multi_time = repmat(v_multi, [fluo_lengths{i_tr}, 1]);
%                 
%                 naive_counts = naive_count_list_MS2(states(1),:);
%                 state1_count = naive_counts(1);
%                 
%                 if state1_count > 0
%                     for t = 1:(w-1)
%                         v_multi_time(t, 1) = v_multi_time(t, 1)/state1_count*...
%                             abs(state1_count - count_reduction(t));
%                     end
%                 end
%                 
%                 difference_tmp = [fluo_values{i_tr}', -v_multi_time];
%                 difference = sum(difference_tmp,2).^2;
%                 
%                 difference_list(states, :) = ...
%                     repmat(difference', [length(states), 1]);
                
                v_multi = v.*naive_count_list_MS2(states(1),:);
                
                naive_counts = naive_count_list_MS2(states(1),:);
                state1_count = naive_counts(1);
                
                difference_tmp = fluo_values{i_tr}'-sum(v_multi);
                if state1_count > 0
                    for t = 1:(w-1)
                        difference_tmp(t) = difference_tmp(t) + ...
                            v_multi(1) - v_multi(1)/state1_count * ...
                              abs(state1_count - count_reduction(t));
                    end
                end
                
                difference_sq = difference_tmp.^2;
                
                difference_sq_list(states, :) = ...
                    repmat(difference_sq', [length(states), 1]);
            end
            
            eta_list = sqrt(lambda/(2*pi))*...
                exp(-0.5*lambda*difference_sq_list);
            
%             eta_list = sqrt(lambda/(2*pi))*...
%                 exp(-0.5*lambda*(difference_sq_list-max(difference_sq_list)));
            
            % ------------ alpha coefficient matrix calculation --------
            
            % pre-allocation of the alpha matrix
            alpha_matrix = zeros(K^w,fluo_lengths{i_tr});
            
            % scaling coefficients for the alpha matrix
            coeff_alpha = zeros(1,fluo_lengths{i_tr});
            
            % calculates the alpha matrix elements for t = 1
            for i = possible_states_list{1}
                alpha_matrix(i, 1) = eta_list(i, 1)* ...
                    pi0(digit_first_list(i));
            end
            coeff_alpha(1) = sum(alpha_matrix(possible_states_list{1}, 1));
            alpha_matrix(possible_states_list{1}, 1) = ...
                alpha_matrix(possible_states_list{1}, 1) / coeff_alpha(1);
        
            % calculates the alpha matrix elements for t > 1
            for t = 2:fluo_lengths{i_tr}
                % possible states at time t
                i_possible = possible_states_list{t};

                % list of terms that are added to find the alpha matrix
                % elements
                alpha_terms_list = alpha_A_list(i_possible, :).* ...
                    reshape(alpha_matrix(allowed_to_list(i_possible,:),t-1), ...
                        length(i_possible), K);
                
                % assignment of alpha matrix element values
                alpha_ls_tmp = sum(alpha_terms_list,2).* ...
                                 eta_list(i_possible, t);
                coeff_alpha(t) = sum(alpha_ls_tmp);
                alpha_matrix(i_possible,t) = alpha_ls_tmp / coeff_alpha(t);
            end
            
            % ------------ beta coefficient matrix calculation ---------
            
            % pre-allocation of the beta coefficient matrix
            beta_matrix = zeros(K^w,fluo_lengths{i_tr});
            
            % scaling coefficients for the beta matrix
            coeff_beta = zeros(1,fluo_lengths{i_tr});

            % assigns 1 to beta matrix elements at t=time
            beta_matrix(:, fluo_lengths{i_tr}) = 1;
            coeff_beta(fluo_lengths{i_tr}) = sum(beta_matrix(:, fluo_lengths{i_tr}));
            beta_matrix(:, fluo_lengths{i_tr}) = ...
                beta_matrix(:, fluo_lengths{i_tr})/coeff_beta(fluo_lengths{i_tr});

            % calculates the alpha matrix elements for t < time
            for t = (fluo_lengths{i_tr}-1: -1: 1)
                % possible states at time t
                i_possible = possible_states_list{t};

                % list of terms that are added to find the beta matrix elements
%                 beta_terms_list = beta_A_list(i_possible, :).* ...
%                     reshape(eta_list(allowed_from_list(i_possible,:), t+1), [], K).* ...
%                     reshape(beta_matrix(allowed_from_list(i_possible,:), t+1), [], K);
                beta_terms_list = beta_A_list(i_possible, :).* ...
                    reshape(eta_list(allowed_from_list(i_possible,:), t+1).* ...
                    beta_matrix(allowed_from_list(i_possible,:), t+1), length(i_possible), K);
                
                % assignment of beta matrix element values
                beta_ls_tmp = sum(beta_terms_list,2);
                coeff_beta(t) = sum(beta_ls_tmp);
                beta_matrix(i_possible,t) = beta_ls_tmp / coeff_beta(t);
            end
            
            % --------------- log-likelihood calcuation --------------------
            log_likelihoods{i_tr}(baum_welch) = ...
                sum(log(coeff_alpha)) + log(sum(alpha_matrix(:,fluo_lengths{i_tr})));
            
%             log_likelihoods{i_tr}(baum_welch) = ...
%                 sum(log(coeff_alpha)) - ...
%                 0.5*lambda*sum(max(difference_sq_list)) + ...
%                 log(sum(alpha_matrix(:,fluo_lengths{i_tr})));
	
            logL_tot(baum_welch) = logL_tot(baum_welch) + ...
                log_likelihoods{i_tr}(baum_welch);
            
            % --------------------- <S_t> calculation ----------------------
            p_s = alpha_matrix.*beta_matrix;
            if (K>1)
                p_s = p_s./sum(p_s);
            end
                
            % ------------------- <S_t, S_{t-1}> calculation ---------------
            
            % replication of the alpha matrix K times along one of the 2d axes
            alpha_matrix_minus = alpha_matrix(:,1:(fluo_lengths{i_tr}-1));
            alpha_matrix_minus_rep = repmat(alpha_matrix_minus, [K, 1]);

            % 2d indexing of eta_log and beta_log matrices
            % note: the 2d slice components are accounted for through the
            %       1d representation of the matrix of possible transitions.
            %       Also, note that the 1d indexing of the p_ss_log matrix
            %       is done by rows, which is used in calculating the terms
            eta_list_minus = eta_list(allowed_from_list(:), 2:fluo_lengths{i_tr});
            beta_matrix_minus = beta_matrix(allowed_from_list(:), 2:fluo_lengths{i_tr});

            % row and column subscripts used for indexing A_log elements
            A_rowSubs = digit_first_list(allowed_from_list(:));
            A_colSubs = digit_first_list(repmat(1:K^w, [1, K]));

            % A element indexing in a 2d slice
            A_sub_single = A(sub2ind([K, K], A_rowSubs, A_colSubs));

            % A element indexing in a 3d slice along the time axis
            A_colSubs_multi = repmat(A_sub_single, [1, fluo_lengths{i_tr}-1]);                

            % calcuation of p_ss elements
            p_ss = alpha_matrix_minus_rep(:).* ...
                eta_list_minus(:).* beta_matrix_minus(:).* ...
                A_colSubs_multi(:);
            p_ss_reshape = reshape(p_ss, [K^(w+1),fluo_lengths{i_tr}-1]);
            p_ss_normalize = repmat(sum(p_ss_reshape, 1),[K^(w+1),1]);
            p_ss = p_ss./p_ss_normalize(:);
            
            % ---- Updates of the expectation terms for each trace -----

            % pi0 terms
            for m = 1:K
                pi0_terms(m) = pi0_terms(m) + sum(C(m,:).*p_s(:,1)');
            end

            % A terms
            for m = 1:K
                for n = 1:K
                    A_terms(m,n) = A_terms(m,n) + sum(p_ss(A_maximization_ind{m,n}));
                end
            end

            % lambda terms
            lambda_terms = lambda_terms + sum(sum(p_s.*difference_sq_list));
            
            % v terms
            for m = 1:K
                for n = 1:K
                    terms_ith = p_s.*F_terms{n}.*F_terms{m};
                    v_M_terms(m,n) = v_M_terms(m,n) + sum(sum(terms_ith));
                    
                    v_M_terms(m,n) = v_M_terms(m,n) - sum(sum(terms_ith));
                    v_M_terms(m,n) = v_M_terms(m,n) + sum(terms_ith(:));
                end
            end

            x_term = repmat(fluo_values{i_tr}, K^w, 1);
            
            for m = 1:K
                terms_b_ith = x_term.*p_s.*F_terms{m};
                v_b_terms(m) = v_b_terms(m) + sum(sum(terms_b_ith));
                
                v_b_terms(m) = v_b_terms(m) - sum(sum(terms_b_ith));
                v_b_terms(m) = v_b_terms(m) + sum(terms_b_ith(:));
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%% MAXIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%

        % ---------------------------- pi0 -----------------------------
        pi0_old = pi0;
        pi0 = pi0_terms/n_traces;
        
        pi0_norm_rel_change = ...
            abs(norm(pi0_old) - norm(pi0))/norm(pi0_old);
        
        % --------------------------- A_log ----------------------------
        A_old = A;
        A = A_terms;

        for n = 1:K
            A(:,n) = A(:,n)/sum(A(:,n));
        end

        A_norm_rel_change = ...
            abs(norm(A_old) - norm(A))/norm(A_old);
        
        % ------------------------- lambda_log -------------------------
        noise_old = 1/sqrt(lambda);
        
        lambda = sum(cell2mat(fluo_lengths))/lambda_terms;
        noise = 1/sqrt(lambda);
        
        noise_change =abs(noise-noise_old);
        noise_rel_change = noise_change / noise_old;
        
        % ----------------------------- v ------------------------------
        v_old = v;
        
        m_sign = ones(K, K);
        m_log = log(v_M_terms);

        b_sign = sign(v_b_terms);
        b_log = log(abs(v_b_terms));

        v_updated = v_log_solve(m_log, m_sign, b_log, b_sign);
        v_logs = v_updated(1,:);
        v_signs = v_updated(2,:);
        v = v_signs.*exp(v_logs);

        v_norm_change = abs(norm(v_old) - norm(v));
        v_norm_rel_change = v_norm_change / norm(v_old);
        
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
    local_em_outputs = struct('pi0_log', log(pi0)', 'A_log', log(A), ...
        'v_logs', log(abs(v)), 'v_signs', sign(v), 'lambda_log', log(lambda), ...
        'logL_bw', logL_tot, 'logL', max(logL_tot), ...
        'n_iter', baum_welch, 'runtime', runtime);