% Script to Conduct Windowed HMM Inference 
% This script reads in interpolated time series data and conducts mhmmm
% inference for the memory (w), state count (K), and data groupings
% (bin_groups) specified
%%%---------------------------------------------------------------------%%%
% Key variables that must be present in data structure:
    % 1) time_interp: interpolated time vector
    % 2) fluo_interp: interpolated fluorescence vector (experimentally observed
    %    quanity)
    % 3) group_vec_interp: vector indicating trace group membership at each
    %    time point. For instance, in Drosophila, this is some quantity derived
    %    from the position a trace along the AP axis.
    % 4) alpha_frac: fractional length of MS2 loops relative to full
    %    transcript. For instance, if MS2 length is 1302 and transcript
    %    length is 6444, alpha_frac = 1302/6444
    % 5) ParticleID: unique trace (particle) identifier
    
%%%---------------------------------------------------------------------%%%
clear 
close all
addpath('../utilities'); % Route to utilities folder

%-------------------------------System Vars-------------------------------%
bin_groups = [-1,1]; % cell array containing region grouping for inference
w = 6; % Memory
K = 2; % State to use for inference
project = 'example_trace_set';  %identifier pointing to desired data set

%------------------Define Inference Variables------------------------------%
n_localEM = 25; % set number of distinct local runs (do not drop below 20)
n_steps_max = 500; % set max steps per inference (500 typically sufficient)
eps = 1e-4; % set convergence criteria (keep at 1e-4)
MaxWorkers = 4; % maximum number of parpool workers that can be supported (will vary by machine)

%----------------------------Bootstrap Vars-------------------------------%
n_bootstrap = 10; % number of unique bootstrap samples to use for error estimation
sample_size = 5000; % N data points to use 
min_dp_per_inf = 1250; % inference will be aborted if fewer present (do not go below 1000)                                         

%-------------------Load Data and Set Write Paths-------------------------%
datapath = ['../../dat/' project '/']; %Path to inference data
dataname = 'trace_struct_final.mat'; %name of inference set

% Load data for inference into struct named: trace_struct_final
load([datapath dataname]);
alpha = trace_struct_final(1).alpha_frac*w; % Rise Time for MS2 Loops in time steps
Tres = trace_struct_final(1).Tres; % Time Resolution

% get vector of group IDs 
group_id_vec = [trace_struct_final.group_id];
bin_groups = unique(group_id_vec);
% Set write path (inference results are now written to external directory)
out_dir =  ['../../out/' project '/w' num2str(w) '_t' num2str(Tres)...
    '_alpha' num2str(round(alpha*10)) '/K' num2str(K) '/']; 
mkdir(out_dir);

%% Conduct Inference
rng('shuffle'); % ensure random number generator initialized at random point
for g = 1:length(bin_groups) % loop through different AP groups
    bin_id = bin_groups(g); % get groups for present iteration                   
    for b = 1:n_bootstrap % iterate through bootstrap replicates
        iter_start = now;
        % structures of store inference info
        local_struct = struct; % track results of each local run            
        output = struct; % final output structure saved to file

        % Use current time as unique inference identifier (kind of unwieldly)
        inference_id = num2str(round(10e5*now));

        % Generate filenames            
        fName_sub = ['inference_set_w' num2str(w) '_K' num2str(K) ...
            '_bin' num2str(bin_id) '_t' inference_id];                
        out_file = [out_dir '/' fName_sub];            

        % Extract fluo_data
        trace_filter = bin_id==group_id_vec;            
        trace_ind = find(trace_filter);
        inference_set = trace_struct_final(trace_filter); % find eligible traces
       
        skip_flag = 0;
        if ~isempty(inference_set)
            set_size = length([inference_set.fluo_interp]);                 
        end
        if isempty(inference_set)
            skip_flag = 1;
        elseif set_size < min_dp_per_inf                    
            skip_flag = 1;                    
        end
        if skip_flag
            warning('Too few data points. Skipping')                
        else % proceed with inference if there is sufficient data
            sample_index = 1:length(inference_set);
            
            % % reset bootstrap size to be on order of set size for small bins
            ndp = 0;    
            sample_ids = [];                                
            if set_size < sample_size
                sample_size = ceil(set_size/100)*100;
            end
            
            % randomly select traces to use for this bootstrap
            while ndp < sample_size
                tr_id = randsample(sample_index,1);
                sample_ids = [sample_ids tr_id];
                ndp = ndp + length(inference_set(tr_id).time_interp);
            end                
            fluo_data = cell([length(sample_ids), 1]); % traces used for inference
            sample_particles = [inference_set(sample_ids).ParticleID];
            for tr = 1:length(sample_ids)
                fluo_data{tr} = inference_set(sample_ids(tr)).fluo_interp;                                        
            end 

            % random initialization of parameters for iid inference
            param_init = initialize_random (K, w, fluo_data);
            % approximate inference assuming iid data for param initialization     
            % iid results used as basis for full inference
            local_iid_out = local_em_iid_reduced_memory(fluo_data, param_init.v, ...
                                param_init.noise, K, w, alpha, n_steps_max, eps);
            noise_iid = 1/sqrt(exp(local_iid_out.lambda_log));
            v_iid = exp(local_iid_out.v_logs);         

            % we can do local EM runs in parallel
            % first check to see if distributed computing is installed
            pool_check = ver('distcomp');
            if ~isempty(pool_check) 
                p = gcp('nocreate');
                if isempty(p)
                    parpool(MaxWorkers); 
                elseif p.NumWorkers > MaxWorkers
                    delete(gcp('nocreate')); % if pool with too many workers, delete and restart
                    parpool(MaxWorkers);
                end
            end
            parfor i_local = 1:n_localEM % Parallel Local EM                
                % random initialization of model parameters
                param_init = initialize_random_with_priors(K, noise_iid, v_iid);
                % get intial values
                pi0_log_init = log(param_init.pi0);
                A_log_init = log(param_init.A);
                v_init = param_init.v;                        
                noise_init = param_init.noise;                    
                %--------------------LocalEM Call-------------------------%
                local_out = local_em_MS2_reduced_memory_truncated(fluo_data, ...
                    v_init, noise_init, pi0_log_init', A_log_init, K, w, ...
                    alpha, n_steps_max, eps);                    
                %---------------------------------------------------------%                
                % Save Results 
                local_struct(i_local).inference_id = inference_id;
                local_struct(i_local).subset_id = i_local;
                local_struct(i_local).logL = local_out.logL;                
                local_struct(i_local).A = exp(local_out.A_log);
                local_struct(i_local).v = exp(local_out.v_logs).*local_out.v_signs;
                local_struct(i_local).r = exp(local_out.v_logs).*local_out.v_signs / Tres;                                
                local_struct(i_local).noise = 1/exp(local_out.lambda_log);
                local_struct(i_local).pi0 = exp(local_out.pi0_log);
                local_struct(i_local).total_steps = local_out.n_iter;               
                local_struct(i_local).soft_struct = local_out.soft_struct;               
            end 

            [logL, max_index] = max([local_struct.logL]); % Get index of best result                    
            % Save parameters from most likely local run
            output.pi0 =local_struct(max_index).pi0;                        
            output.r = local_struct(max_index).r(:);                
            output.noise = local_struct(max_index).noise;
            output.A = local_struct(max_index).A(:);
            output.A_mat = local_struct(max_index).A;                            
            % Info about run time
            output.total_steps = local_struct(max_index).total_steps;                                  
            output.total_time = 100000*(now - iter_start);                                                         
            % other inference characteristics
            output.group_id = bin_id; % inf group info                                                      
            output.iter_id = b;                                                 
            output.particle_ids = sample_particles;
            output.N = ndp;                
            output.w = w;
            output.alpha = alpha;
            output.Tres = Tres;                 
        end
        output.skip_flag = skip_flag;
        save([out_file '.mat'], 'output');           
    end      
end    
