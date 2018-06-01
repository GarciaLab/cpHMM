function param_init = initialize_random (K, w, fluo_data)
    % Returns random initialization values for the model parameters using
    % the fluorescence data from multiple traces.
    % 
    % INPUTS
    % K: number of naive states    
    % w: memory
    % fluo_data: call variable of the fluorescence sequences for all traces
    % 
    %
    % OUTPUTS
    % param_init: a structure with randomly initialized values for the 
    %             model parameters
    % param_init.pi0: initial probability distribution of naive states
    % param_init.A: transition matrix
    % param_init.v: emission values
    % param_init.noise: Gaussian noise

    % random pi0 generation
    pi0 = gamrnd(ones(1,K),1,1,K)';
    pi0 = pi0/sum(pi0);
    
    % random A generation
    A = zeros([K,K]);
    for j = 1:K
        A(:,j) = gamrnd(ones(1,K),1,1,K)';
        A(:,j) = A(:,j)/sum(A(:,j));
    end
    
    fluo_values_list = [];
    n_traces = length(fluo_data);
    for i_tr = 1:n_traces
        fluo_values_list = [fluo_values_list, ...
            fluo_data{i_tr}(:)'];
    end
    
    % random v generation
    v = sort(rand(1,K)) * (2/w) * max(fluo_values_list(:));
    
    % random noise genertion
    noise = rand() * mean(fluo_values_list(:));
    
    % combine the initialized values into a list
    param_init = struct;
    param_init.pi0 = pi0;
    param_init.A = A;
    param_init.v = v;
    param_init.noise = noise;