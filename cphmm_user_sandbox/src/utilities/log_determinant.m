function [log_det] = log_determinant (m_lg, m_sgn)
    % The determinant of the matrix given in log form.
    % 
    % INPUTS
    % m_lg: log of the absolute values of the matrix elements
    % m_sgn: signs of the matrix elements
    % 
    % OUTPUTS
    % log_det: a 2-element array of the sign and the log of the determinant
    
    [n, ~] = size(m_lg);     
    perm_list = perms(1:n);
    [perm_num, ~] = size(perm_list);
       
    logs = zeros(perm_num, 1);
    signs = zeros(perm_num, 1);

    for i = 1:perm_num        
        rows = 1:n;
        cols = perm_list(i, :);
        
        ind = sub2ind(size(m_sgn), rows, cols);

        signs(i) = permutation_parity(cols) * prod(m_sgn(ind));
        logs(i) = sum(m_lg(ind));        
    end

    log_det = log_sum_exp(logs, signs);