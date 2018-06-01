function viterbi_out = viterbi (fluo_values, v, noise, pi0_log, ...
                                     A_log, K, w, alpha_var)

    % Given a time series of fluorescence values and the inferred model
    % parameters, returns the most likely path.
    % 
    % INPUTS
    % fluo_values: an array of time series fluorescence values
    % v: emission values
    % noise: Gaussian noise
    % pi0_log: log of the initial naive state pmf
    % A_log: log of the transition probability matrix
    % K: number of promoter states
    % w: memory
    % alpha: length of the MS2 loop in time steps
    %
    % OUTPUTS
    % viterbi_out: structure with the Viterbi algorithm results
    % viterbi_out.z_viterbi: Viterbi path of naive states
    % viterbi_out.v_viterbi: Viterbi path of the emission states
    % viterbi_out.fluo_viterbi: Fluorescence trace corresponding to the
    %                           Viterbi path
    
    % length of the fluorescence series
    seq_length = length(fluo_values);

    % fraction of the full MS2 sequence transcribed at each elongation step
    ms2_coeff = ms2_loading_coeff(alpha_var,w);
    
    % ----------------------------- Lists ------------------------------
    
    % list of states allowed to go from each of 1:K^w compound states
    allowed_from_list = zeros(K^w, K, 'int16');
    
    % list of the first digits (naive states) of each compound state
    digit_first_list = zeros(1, K^w, 'int16');
    
    % list of the number of each naive state encountered in compound states
    naive_count_list_MS2 = zeros(K^w, K);
    
    for i = 1:K^w
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
    
    % list that maps the unique naive count combinations to compound states
    naive_count_map = cell(n_unique, 1);   
    for i = 1:n_unique
        naive = unique_naive_MS2_list(i,:);
        naive_count_map{i} = find(ismember(naive_count_list_MS2,naive,'rows'));        
    end
    
    % list of possible compound states at all times
    possible_states_list = cell(seq_length, 1);    
    for t = 1:seq_length
        possible_states_list{t} = possible_states (t, K, w);        
    end
    
    % ---------------------- Variable assignments ----------------------
    % logs of v and noise
    v_logs = log(abs(v));
    v_signs = sign(v);
    lambda_log = -2*log(noise);
    
    % log and sign of the fluorescence values
    fluo_log = log(abs(fluo_values));
    fluo_sign = sign(fluo_values);
    
    % accounting variable
    B = zeros([K^w, K^w]);
    for i = 1:K^w
        for j = 1:K^w
            if find(allowed_from_list(j,:) == i)
                B(i,j) = 1;
            end
        end
    end
    
    % --------------------- Log of emission terms ---------------------
    difference_list = zeros(K^w, seq_length);
    
    for i = 1:n_unique
        states = naive_count_map{i};
        difference_list(states, :) = ...
            repmat(difference(fluo_log, fluo_sign, seq_length, K, w, alpha_var, ...
            states(1), v_logs', v_signs', naive_count_list_MS2_log)', [length(states), 1]);        
    end
    eta_log = 0.5*(lambda_log - log(2*pi)) ...
        -0.5*exp(lambda_log + difference_list);
    
    % --------------------- Viterbi path calculation -------------------
    
    % viterbi variable in log scale
    V_log = zeros([seq_length, K^w]);
    V_log(:,:) = -Inf;
    
    for j = possible_states (1, K, w)
        V_log(1, j) = eta_log(j, 1) + pi0_log(digit_first_list(j));
    end
    
    for t = 2:seq_length
        prev_states = possible_states_list{t-1};
        for i = possible_states_list{t}
            V_log(t, i) = eta_log(i, t) + max(V_log(t-1, prev_states) + ...
                A_log(digit_first_list(i), digit_first_list(prev_states)) + ...
                log(B(i, prev_states)));
        end
    end
    
    % Viterbi path in compound states
    S_viterbi = zeros([1, seq_length]);
    
    ind_max = find(V_log(seq_length, :) == max(V_log(seq_length, :)));
    S_viterbi(seq_length) = ind_max(1);
    
    for t = (seq_length-1):-1:1
        prev_states = possible_states_list{t};
        terms = V_log(t, prev_states) + log(B(S_viterbi(t+1), prev_states)) + ...
            A_log(digit_first_list(S_viterbi(t+1)), digit_first_list(prev_states));
        ind_max = find(terms == max(terms));
        S_viterbi(t) = prev_states(ind_max(1));
        if (length(ind_max) > 1)
            warning('Viterbi path is not uniquely defined');
        end
    end
    
    % viterbi path in naive states
    z_viterbi = digit_first_list(S_viterbi);
    
    % corresponding emission values
    v_viterbi = v(z_viterbi);

    % create a shifted emission matrix for fluorescence calculation
    emissions_mat = zeros(seq_length, w);
    for j = 1:w
        i_start = min([j, seq_length]);
        i_end_emission = seq_length-j+1;
        emissions_mat(i_start:end,j) = v_viterbi(1:i_end_emission);
    end

    % fluorescence with MS2 loops taken into account
    fluo_viterbi = ms2_coeff * transpose(emissions_mat);
    
    % output structure of the Viterbi results
    viterbi_out = struct('z_viterbi', z_viterbi, ...
                         'v_viterbi', v_viterbi, ...
                         'fluo_viterbi', fluo_viterbi);