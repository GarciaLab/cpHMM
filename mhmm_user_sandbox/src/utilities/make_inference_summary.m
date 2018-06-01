function [header, summary] = make_inference_summary (global_out)

%-----------------------------Make Summary File-----------------------%
    
    % filter for min # data points
    % global_out = global_out(data >= dtMin);

    grID = global_out(1).grID;
    K = global_out(1).K;
    deltaT = global_out(1).deltaT;
    ms = [global_out.ms]';
    
    if min(ms) == 1
        gr = [global_out.gr]';
    else
        gr = floor([global_out.gr]');
    end

    % generate arrays with relevant inference variables
    pi0 = [global_out.pi0]';
    v = [global_out.v]';
    noise = [global_out.noise]';
    A = [global_out.A]';
    R = [global_out.R]';
    
    % keep track of which groups required rate fitting 
    R_fit_flag = zeros(size(R,1),1);
    
    for i = 1:size(R,1)
        Rvec = R(i,:);
        global_out(i).R_orig = Rvec;
        
        % if rates contain complex values or more than K negative rates,
        % find best-fitting "non-degenerate" rate matrix 
        if ~isreal(Rvec)||(sum(Rvec<0)>K)
            out = prob_to_rate_fit_sym(reshape(A(i,:),K,K), deltaT, 'gen', .005, 1);
            R(i,:) = reshape(out.R_out,1,[]);            
            R_fit_flag(i) = 1;
        end
    end

    gr_count = size(A,1);
    for i = 1:gr_count
       % rank initaition rates by size
       [v(i,:), ranked_v] = sort(v(i,:));
       
       % reorder outputs accordingly
       pi0(i,:) = pi0(i,ranked_v);
       Amat = reshape(A(i,:),K,K);
       A(i,:) = reshape(Amat(ranked_v, ranked_v),1,[]);
       Rmat = reshape(R(i,:),K,K);
       R(i,:) = reshape(Rmat(ranked_v, ranked_v),1,[]);       
    end

    gr = double(gr);
    summary = [gr, ms, pi0, v, noise, A, R, R_fit_flag];

    header = [{grID} {'ms'}];
    for i = 1:K
        header = [header {['pi' int2str(i)]}];
    end
    for i = 1:K
        header = [header {['v' int2str(i)]}];
    end
    header = [header {'noise'}];

    for j = 1:K
        for i = 1:K
            header = [header {['A' int2str(i) int2str(j)]}];
        end
    end

    for j = 1:K
        for i = 1:K
            header = [header {['R' int2str(i) int2str(j)]}];
        end
    end
    
    header = [header 'RfitFlag'];