function v_sln = v_log_solve (m_lg, m_sgn, b_lg, b_sgn)
    % Solves the system Mv=b for v.
    % 
    % INPUTS
    % m_lg: log of the absolute value of the matrix M
    % m_sgn: signs of the elements of the matrix M
    % b_lg: log of the absolute value of the vector b
    % b_sgn: signs of the elements of vector b
    % 
    % OUTPUTS
    % v: the solution of the system Mv=b for the unknown v. Returns v in 
    %    log form and avoids underflow issues.
    
    [n, ~] = size(m_lg);
    v_lgs = zeros(n, 1);
    v_sgns = zeros(n, 1);

    det_m = log_determinant(m_lg, m_sgn);    
    
    for j = 1:n
        m_log_j = m_lg;
        m_log_j(:,j) = b_lg;

        m_sgn_j = m_sgn;
        m_sgn_j(:,j) = b_sgn;

        det_j = log_determinant(m_log_j, m_sgn_j);
        v_lgs(j) = det_j(1) - det_m(1);
        v_sgns(j) = det_j(2) * det_m(2);        
    end
    
    v_sln = vertcat(transpose(v_lgs), transpose(v_sgns));