function steady = occupancies(A)
    % Given the transition probability matrix, returns the steady state
    % occupancies of the states
    
    eps = 1e-4;
    [V, D] = eig(A);
    [n, ~] = size(A);
    [~, j] = ind2sub([n,n], find( abs(D-1) < eps ));
    steady = V(:,j)./sum(V(:,j));