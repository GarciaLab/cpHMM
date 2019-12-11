function states = possible_states(t, K, w)
    % Returns a list of all possible compound states at time t.
    % 
    % INPUTS
    % t: time
    % K: number of naive states
    % w: memory
    % 
    % OUTPUTS
    % states: all possible compounds states at time t
    
    if (t >=w)
        states = 1:K^w;
    else
        states = 1 + (0:(K^t-1)) * K^(w-t);
    end