function param_init = initialize_random_with_priors (K, noise, v)
    % Returns random initialization values for the HMM model, given
    % approximate estimates of the noise and emission parameters from the
    % i.i.d. model.
    % 
    % INPUTS
    % K: number of naive states        
    % noise: Gaussian noise
    % v: emission values
    % 
    % OUTPUTS
    % param_init: a structure with randomly initialized values for the 
    %             model parameters
    % param_init.pi0: initial probability distribution of the naive states
    % param_init.A: transition matrix
    % param_init.v: emission values
    % param_init.noise: Gaussian noise

    % structure array to store the initialization parameters
    param_init = struct;
    
    % random pi0 generation
    pi0 = gamrnd(ones(1,K),1,1,K)';
    pi0 = pi0/sum(pi0);
    param_init.pi0 = pi0;
    
    % random A generation
    A = zeros([K,K]);
    for j = 1:K
        A(:,j) = gamrnd(ones(1,K),1,1,K)';
        A(:,j) = A(:,j)/sum(A(:,j));
    end
    param_init.A = A;
    
    % random noise generation
    noise_min = 0.5 * noise;
    noise_max = 2.0 * noise;
    noise_range = noise_max - noise_min;
    param_init.noise = noise_min + rand() * noise_range;
    
    % random v generation
    v_max = max(v) * (0.7 + 0.6*rand());
    if K==1
        param_init.v = rand()*v_max;
    else
        param_init.v = [0, sort(rand(1, K-2)), 1] * v_max;
    end