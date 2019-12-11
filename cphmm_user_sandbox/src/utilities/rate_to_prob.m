function A = rate_to_prob (R, deltaT)
    % Calculates the transition probability matrix using the transition 
    % rate matrix
    % 
    % INPUTS
    % R: transition rate matrix
    % deltaT: time interval between consecutive transitions
    % 
    % OUTPUTS
    % A: transition probability matrix
    
    A = expm(R*deltaT);