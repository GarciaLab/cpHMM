function d = digit(state, position, K, w)
    % Returns the naive state at the given position.
    % 
    % Given a compound state, finds the corresponding naive state sequence,
    % and returns the naive state at the given position.
    % 
    % INPUTS
    % state: compound state
    % position: position in naive state sequence
    % K: number of naive states
    % w: memory
    % 
    % OUTPUTS
    % d: the naive state at the specified position
    
    remainder = mod(state-1, K^(w - position + 1));
    d = floor(remainder / K^(w - position)) + 1;