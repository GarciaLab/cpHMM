function naive = compound_to_naive(state, K, w)
    % Given a compound state, returns the corresponding naive state 
    % sequence.
    % 
    % INPUTS
    % state: compound state
    % K: number of naive states
    % w: memory
    % 
    % OUTPUTS
    % naive: naive state sequence of length w (memory), where each naive 
    %        state can have values in the range [1:K].
    
    remainders = mod(state-1, K.^(w - (0:w-1)));
    naive = floor(remainders./ K.^(w - (0:w-1) - 1)) + 1;