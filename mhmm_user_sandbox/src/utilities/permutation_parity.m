function parity = permutation_parity (lst)        
    % INPUTS
    % lst: a permutation of numbers from 1 to max(lst)
    %
    % OUTPUTS
    % parity: the parity of the list of numbers from 1 to max(lst)
    
    I = speye(length(lst));
    parity = det(I(:,lst));