function coeff_int = ms2_loading_coeff_integral (alpha, w, deltaT, t1, t2)
    
    % Finds the scaling coefficient kappa which accounts for the reduced
    % fluorescence caused by the finite size of the MS2 loops
    %
    % INPUTS
    % alpha: length of the MS2 loop in time steps
    % deltaT: time resolution [sec]
    % w: memory
    % t1: start of the promoter activity, counting from moment of
    %     polymerase initiation (0 <= t1 <= w*deltaT)
    % t2: end of the promoter activity, counting from the moment of
    %     polymerase initiation (0 <= t1 <= t2 <= w*deltaT)
    % 
    % OUTPUTS
    % coeff_int: scaling coefficient
    
    % small number used in checking conditions
    eps = 10^-5;
    
    % throw an error if the conditions for t1 and t2 are not met
    if (t1 < 0 || t1 -t2 > eps || t2 - w*deltaT > eps)
        error('The condition 0 <= t1 <= t2 <= w*deltaT is not met');
    end
    
    % throw an error when the size of the MS2 loop is larger than memory
    % or when it's negative
    if (alpha > w || alpha < 0)
        error('The condition 0 <= alpha <= w is not met');
    end
    
    % MS2 loop transcription time
    t_MS2 = alpha*deltaT;
    
    % the interval is in the MS2 transcription region
    if (t2 <= t_MS2)
        coeff_int = (t2^2 - t1^2)/(2*t_MS2);
        
    % the interval is after the MS2 transcription region
    elseif (t1 >= t_MS2)
        coeff_int = t2-t1;
        
    % the interval starts at the transcription region but ends after it
    else
        coeff_int = 0;
        coeff_int = coeff_int + (t_MS2^2 - t1^2)/(2*t_MS2);
        coeff_int = coeff_int + (t2-t_MS2);
    end