function out = prob_to_rate_fit_sym (A, dT, type, prior, u_bound)
    % This function uses nonlinear least squares fitting to find rate
    % matrix of proper form that corresponds to the probbility matrix 
    % obtained from HMM as closely as possible.
    
    % Three rate matrix formats: 
    %       'gen': no assumptions about rate values. 6 rates obtained as
    %              independent parameters
    %       'ladder': Assumes system cannot jump between low and high
    %                 states
    %       '2state': Assumes unitary Off and On rate for matrix
    %       construction
    % 
    % INPUTS:
    % A: probability matrix to be fit
    % dT: Time step between observations
    % type: string indicating which type of rate matrix to fit
    % prior: initialization value for rate fitting MUST BE POSITIVE
    % u_bound: upper bound for rates obtained through fitting
    
    K = size(A,1);
    k = sym('k%d_%d', [K K]);
    
    if strcmp(type,'gen')
        typX = ones(1,K^2-K)*prior;
    elseif strcmp(type,'ladder')
        typX = ones(1,K^2-K-2)*prior;
    elseif strcmp(type,'2state')
        typX = ones(1,2)*prior;
    end
    
    % Specify options to be incorporated into least squares fitting
    options = optimset('MaxFunEvals', 1000 ,'TolFun',1e-15, 'TolX',1e-15,...
                 'TypicalX',typX);
             
    if strcmp(type,'gen')
        k_init = k - diag(diag(k));
        k_fun = matlabFunction(k_init - diag(sum(k_init)));
        k_gen_cost = @(x) expm(call_on(k_fun,x)*dT) - A;
        
        k_bound_l = zeros(1,K^2-K);
        k_bound_u = ones(1,K^2-K)*u_bound;
        k_prior = (ones(1,K^2-K)*prior);
    end

    if strcmp(type,'ladder') && (K > 2)
        k_init = k - diag(diag(k));
        k_init(1,K) = 0;
        k_init(K,1) = 0;
        
        k_fun = matlabFunction(k_init - diag(sum(k_init)));
        k_gen_cost = @(x) expm(call_on(k_fun,x)*dT) - A;
        
        k_prior = (ones(1,K^2-K-2)*prior);
        k_bound_l = zeros(1,K^2-K-2);
        k_bound_u = ones(1,K^2-K-2)*u_bound;    
    end

    if strcmp(type,'2rate') && (K ==3)
        k_init = k - diag(diag(k));
        k_init(1,K) = 0;
        k_init(K,1) = 0;
        k_init(2,1) = 2*k_init(3,2);
        k_init(2,3) = 2*k_init(1,2);

        k_fun = matlabFunction(k_init  - diag(sum(k_init)));
        k_gen_cost = @(x) expm(call_on(k_fun,x)*dT) - A;
        
        k_bound_l = zeros(1,2);
        k_bound_u = ones(1,2)*u_bound;
        k_prior = (ones(1,2)*prior);
    end
    
    [x1,~,out.residuals] = lsqnonlin(k_gen_cost,k_prior,k_bound_l,k_bound_u,options); 
    input = num2cell(x1);
    out.R_out = k_fun(input{:});