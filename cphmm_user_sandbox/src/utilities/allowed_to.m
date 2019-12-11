function l = allowed_to(state, K, w)   
    % List of all compound states that can transition into the given 
    % compound state.
    % 
    % INPUTS
    % state: compound state
    % K: number of naive states
    % w: memory
    % 
    % OUTPUTS
    % l : the list of all compound states that can transition into the 
    %     given compound state. Note: the fact that should be an overlap of 
    %     (w-1) naive states in the past and present is used in the 
    %     calculations.
    
    first_naive = floor( (state-1) / K^(w-1));
    tail_naive = state-1 - first_naive * K^(w-1);
    l = 1 + (K*tail_naive : (K*tail_naive + K - 1));