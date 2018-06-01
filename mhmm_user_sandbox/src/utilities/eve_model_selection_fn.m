function output_cv = eve_model_selection_fn(fluo_all_ap, bin, n_local, ...
                     n_steps_max, eps, n_cv, K_cv_ls, cv_sample_size,test_fraction, w)

% INPUTS:
% fluo_all_ap: structure variable of all traces
% bin: the AP position for which model selection is to be done
%
% n_local: number of localEM runs per each CV set
% n_steps_max: maximum number of EM iterations
% eps: tolerance parameter for inference convergence
%
% n_cv: number of MCCV subsets
% K_cv_ls: list of state counts to be used in model selection
% cv_sample_size: size of total cv sample (if data set is smaller will use
% full set)
% test_fraction: fraction of data points to be used in the test sets
% w: (int) number of time steps for which transcript remains on gene

%%

% indices of traces with the specified AP
ind_region = [fluo_all_ap.binID] == bin;

% subset of the structure corresponding to the specified AP
fluo_region = fluo_all_ap(ind_region);

% cell variable of traces in the specified AP
fluo_data = cell([1, length(fluo_region)]);
for i = 1:length(fluo_region)
    fluo_data{i} = fluo_region(i).fluo;
end

% length of the MS2 loop in time steps
kappa = fluo_region(1).alpha_frac*w;

% total number of localEM runs
n_local_total = n_local * length(K_cv_ls) * n_cv;

% total number of points to test (total points in region or specified samp
% size...whichever is less)
n_points_total = min(cv_sample_size, length([fluo_region.fluo]));

% number of points in the testing dataset (ideal count)
n_points_test = round(test_fraction * n_points_total);

% trace indices for the testing dataset for each CV partition
test_ind = cell([n_cv, 1]);

% number of data points in each testing partition
n_points_test_ind = zeros([n_cv, 1]);

% trace indices for the training dataset for each CV partition
train_ind = cell([n_cv, 1]);

%Create index for sampling
sample_index = 1:length(fluo_region);

%% ---------------- Index Training and Testing Datasets ----------------
for ind_cv = 1:n_cv
    n_sr = 0;
    sample_id_list = [];
    s_index_iter = sample_index;
    while n_sr < n_points_total
        s_id = randsample(s_index_iter,1,false);
        sample_id_list = [sample_id_list s_id];
        n_sr = n_sr + length(fluo_region(s_id).fluo);
        s_index_iter = s_index_iter(s_index_iter~=s_id);
    end
    %Sample ids comprise full pool for each iteration
    ind_full_remaining = sample_id_list;
    ind_test_current = [];
    while n_points_test_ind(ind_cv) < n_points_test
        ind = datasample(ind_full_remaining, 1);
        ind_full_remaining = setdiff(ind_full_remaining, ind);
        n_points_test_ind(ind_cv) = n_points_test_ind(ind_cv) + length(fluo_region(ind).fluo);
        ind_test_current = sort([ind_test_current, ind]);
    end
    test_ind{ind_cv} = ind_test_current;
    train_ind{ind_cv} = setdiff(sample_id_list, ind_test_current);
end

%% ------------------------- Index Local Runs --------------------------
output_cv_local = struct;

local_count = 0;
for ind_cv = 1:n_cv
    for K_cv = K_cv_ls
        for i_local = 1:n_local
            local_count = local_count + 1;
            output_cv_local(local_count).ind_cv = ind_cv;
            output_cv_local(local_count).K_cv = K_cv;
            output_cv_local(local_count).i_local = i_local;
            output_cv_local(i_local).local_out = [];
        end
    end
end

%% -------------- Perform i.i.d. runs for initialization ---------------
noise_iid_cell = cell(max(K_cv_ls));
v_iid_cell = cell(1,max(K_cv_ls));

for K_cv = K_cv_ls
    % random initialization of model parameters
    param_init = initialize_random (K_cv, w, fluo_data);
    local_iid_out = local_em_iid_reduced_memory (fluo_data, param_init.v, ...
        param_init.noise, K_cv, w, kappa, n_steps_max, eps);
    
    noise_iid_cell{K_cv} = 1/sqrt(exp(local_iid_out.lambda_log));
    v_iid_cell{K_cv} = sort(exp(local_iid_out.v_logs));
end
%% ------------------ Perform Training LocalEM Runs --------------------
pool = parpool(24);
parfor i = 1:n_local_total
    ind_cv = output_cv_local(i).ind_cv;
    K_cv = output_cv_local(i).K_cv;
    
    test_ind_local = test_ind{ind_cv};
    train_ind_local = train_ind{ind_cv};
    
    % training set    
    fluo_data_train = cell([length(train_ind_local), 1]);
    for j = 1:length(train_ind_local)
        fluo_data_train{j} = fluo_region(train_ind_local(j)).fluo;
    end
    % testing set
    fluo_data_test = cell([length(test_ind_local), 1]);
    for j = 1:length(test_ind_local)
        fluo_data_test{j} = fluo_region(test_ind_local(j)).fluo;
    end
    
    noise_iid = noise_iid_cell{K_cv};
    v_iid = v_iid_cell{K_cv};   
    % random initialization of model parameters
    param_init = initialize_random_with_priors(K_cv, noise_iid, v_iid);
    pi0_log_init = log(param_init.pi0);
    A_log_init = log(param_init.A);
    v_init = param_init.v;
    noise_init = param_init.noise;
    % localEM call
    output_cv_local(i).local_out = local_em_MS2_reduced_memory (fluo_data_train, ...
        v_init, noise_init, pi0_log_init', A_log_init, K_cv, w, ...
        kappa, n_steps_max, eps);
end

%% ----------------------- Find the Global Optima ------------------------

output_cv = struct;
output_size = 0;
logL_max_cv = zeros([length(K_cv_ls), n_cv]);

for ind_cv = 1:n_cv
    
    test_ind_local = test_ind{ind_cv};
    train_ind_local = train_ind{ind_cv};
    
    % training set
    fluo_data_train = cell([length(train_ind_local), 1]);
    for j = 1:length(train_ind_local)
        fluo_data_train{j} = fluo_region(train_ind_local(j)).fluo;
    end
    % testing set
    fluo_data_test = cell([length(test_ind_local), 1]);
    for j = 1:length(test_ind_local)
        fluo_data_test{j} = fluo_region(test_ind_local(j)).fluo;
    end
    
    for K_cv = K_cv_ls
        
        logL_max = -Inf;
        ind_local = find((ind_cv==[output_cv_local.ind_cv]).*...
            (K_cv==[output_cv_local.K_cv]));
        for i = ind_local
            local_out = output_cv_local(i).local_out;
            if local_out.logL > logL_max
                logL_max = local_out.logL;
                A_log_inf = local_out.A_log;
                A_inf = exp(A_log_inf);
                v_inf = exp(local_out.v_logs).*local_out.v_signs;
                lambda_inf = exp(local_out.lambda_log);
                noise_inf = 1/sqrt(lambda_inf);
                pi0_log_inf = local_out.pi0_log;
                pi0_inf = exp(pi0_log_inf);
            end
        end
        
        logL_max_cv(K_cv, ind_cv) = likelihood_reduced_memory (fluo_data_test, v_inf, ...
            noise_inf, pi0_log_inf, A_log_inf, K_cv, w, kappa);
        
        output_size = output_size + 1;
        
        output_cv(output_size).A = A_inf;
        output_cv(output_size).v = v_inf;
        output_cv(output_size).noise = noise_inf;
        output_cv(output_size).pi0 = pi0_inf;
        output_cv(output_size).K = K_cv;
        output_cv(output_size).ind_cv = ind_cv;
        output_cv(output_size).n_points_test_ind = n_points_test_ind(ind_cv);
        output_cv(output_size).n_points_test = n_points_test;
        output_cv(output_size).binID = bin;
        
        output_cv(output_size).logL_test = logL_max_cv(K_cv, ind_cv);
        output_cv(output_size).logL_test_norm = ...
            logL_max_cv(K_cv, ind_cv) * n_points_test / n_points_test_ind(ind_cv);
        output_cv(output_size).logL_train = logL_max;
    end
end