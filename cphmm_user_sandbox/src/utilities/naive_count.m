function counts = naive_count (state, K, w)
    % Returns the number of each naive state at the given compound state.
    % 
    % INPUTS
    % state: compound state
    % K: number of naive states
    % w: memory
    % 
    % OUTPUTS
    % counts : the number of times each naive state is accessed in 
    %          the given compound state.
    
    naive = compound_to_naive(state, K, w);
    counts = zeros(1,K);
    for i = 1:K
        counts(i) = sum(naive==i);
    end