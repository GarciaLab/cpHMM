function synthetic_data = synthetic_prob_regime2(seq_length, kappa, K, w, ...
                             v, kon_min, kon_max, koff, deltaT, noise, pi0)
    
    % Generates a fluorescence sequence using a discrete transition model
    % and the given model parameters.
    % 
    % INPUTS
    % seq_length: length of the sequence
    % kappa: length of the MS2 loop in time steps
    % K: number of naive states
    % w: memory        
    % A: transition matrix
    % v: emission values
    % noise: Gaussian noise
    % pi0: initial cmf of naive states
    % 
    % OUTPUTS
    % synthetic_data: structure that contains the synthetic data info
    %   .fluo: fluorescence sequence with Gaussian noise and the MS2 loop
    %          effect ignored
    %   .fluo_MS2: fluorescence sequence with Gaussian noise and taking the
    %              MS2 loop effect into account
    %   .fluo_no_noise: fluorescence seqeunce with no Gaussian noise and
    %                   the MS2 effect ignored
    %   .fluo_MS2_no_noise: fluorescence sequence with no Gaussian noise
    %                       and the MS2 effect taken into account
    %   .naive_states: sequence of generated naive states
    
    naive_states = zeros(1, seq_length);
    
    % sample from the initial state cmf for the first state
    naive_states(1) = randsample(1:K,1,true,pi0);
    
    R_const = [-2*kon_max,        koff,           0; ...
                2*kon_max,    -(kon_max+koff),   2*koff; ...
                0.000,        kon_max,      -2*koff];
    A_const = rate_to_prob(R_const, deltaT);
    t_const = round(seq_length/3);
    
    % generate a sequence of naive states using transition probabilities
    for t = 2:t_const
        naive_states(t) = ...
            randsample(1:K,1,true,A_const(:,naive_states(t-1)));
    end
    for t = (t_const + 1) : seq_length
        kon = kon_max + (kon_min-kon_max)*(t-t_const)/(seq_length - t_const);
        R_linear = [-2*kon,        koff,           0; ...
                     2*kon,    -(kon+koff),   2*koff; ...
                     0.000,        kon,      -2*koff];
        A_linear = rate_to_prob(R_linear, deltaT);
        
        naive_states(t) = ...
            randsample(1:K,1,true,A_linear(:,naive_states(t-1)));
    end
  
    emissions = v(naive_states);
    
    % create a shifted emission matrix for fluorescence calculation
    % X1  0  0  0  0
    % X2 X1  0  0  0
    % X3 X2 X1  0  0
    % X4 X3 X2 X1  0
    % X5 X4 X3 X2 X1
    % ..............
    emissions_mat = zeros(seq_length, w);
    for j = 1:w
        i_start = min([j, seq_length]);
        i_end_emission = seq_length-j+1;
        emissions_mat(i_start:end,j) = emissions(1:i_end_emission);
    end

    % scaling coefficients that account for the presence of MS2 loops
    coeff_MS2 = ms2_loading_coeff(kappa, w);
    
    % fluorescence with MS2 loops taken into account
    fluo_MS2 = coeff_MS2 * transpose(emissions_mat);
    
    % use w to calculate aggregate emission states
    e_sum = cumsum(emissions);
    fluo = cat(2,e_sum(1:w),(e_sum((w+1):end) - e_sum(1:(end-w))));
    
    % generate Gaussian noise
    gauss_noise = normrnd(0,noise,1,seq_length);
    
    % add the Gaussian noise
    fluo_noise = fluo + gauss_noise;
    fluo_MS2_noise = fluo_MS2 + gauss_noise;
    
    synthetic_data = struct('fluo', fluo_noise, 'fluo_MS2', fluo_MS2_noise, ...
        'fluo_no_noise', fluo, 'fluo_MS2_no_noise', fluo_MS2, ...
        'naive_states', naive_states, 'emissions', emissions);