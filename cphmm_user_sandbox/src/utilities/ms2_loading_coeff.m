function [coeff] = ms2_loading_coeff (alpha, w)
    % Returns the mean fractions of the full MS2 loop transcribed at each
    % step of the elongation window
    % 
    % INPUTS
    % alpha: length of the MS2 loop in time steps
    % w: system memory
    % 
    % OUTPUTS
    % coeff: an array of size w of the mean fractions of the full MS2 loop
    %        transcribed at each step of the elongation window
    
    
    % throw an error when the size of the MS2 loop is larger than memory
    % or when it's negative
    if (alpha > w || alpha < 0)
        error('The condition 0 <= alpha <= w is not met');
    end
    
    coeff = ones(1,w);
    
    if alpha > 0
        alpha_ceil = ceil(alpha);
        alpha_floor = floor(alpha);

        coeff(1:alpha_floor) = ((1:alpha_floor)-0.5)/alpha;
        coeff(alpha_ceil) = (alpha_ceil-alpha) + ...
            (alpha^2 - (alpha_ceil-1)^2)/(2*alpha);
    end