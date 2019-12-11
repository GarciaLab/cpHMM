% Script to Conduct Windowed HMM Inference on Experimental Data
close all
clear 
addpath('../utilities'); % Route to utilities folder
savio = 0; % Specify whether inference is being conducted on Savio Cluster
ap_ref_cell = {-7:-4,-3:-2,-1:1,2:3,4:7}; 
if savio
    %Get environment variable from job script
    savio_groups = {str2num(getenv('SLURM_ARRAY_TASK_ID'))};    
    bin_groups = cell(1,length(savio_groups));
    for i =1:length(bin_groups)
        bin_groups{i} = ap_ref_cell{savio_groups{i}};
    end
else
    bin_groups = {-7:-4,-3:-2,-1:1,2:3,4:7};
end
warning('off','all') %Shut off Warnings
%-------------------------------System Vars-------------------------------%
w = 7; % Memory
Tres = 20; % Time Resolution
K = 3; % State to use for inference
stop_time_inf = 50; % Specify cut-off time for inference
min_dp = 10; % min length of traces to include
clipped = 1; % if 0 use "full" trace with leading and trailing 0's
fluo_field = 1; % specify which fluo field to (1 or 3)
inference_times = 25*60;%(7.5:2.5:40)*60;%fliplr((25:2.5:40)*60);
t_window = 50*60; % determines width of sliding window
clipped_ends = 1; % if one, remove final w time steps from traces
off_traces_flag = 0; % if 1 filter for only traces that are observed to transition into quiescence
                     % if 2 filter and back-align
project = 'mHMMeve2_weka_inf_2018_05_07';  %identifier pointing to desired iteration of data set
%------------------Define Inference Variables------------------------------%
n_localEM = 25; % set num local runs (default=25)
n_steps_max = 500; % set max steps per inference (default=500)
eps = 1e-4; % set convergence criteria (keep at 1e-4)
%----------------------------Bootstrap Vars-------------------------------%
dp_bootstrap = 1; % if 1 use dp bootstrapping
set_bootstrap = 0; % if 1 use set bootstrapping
n_bootstrap = 10; % number of bootstraps (overridden for set bootstrapping)
sample_size = 8000; % N data points to use 
min_dp_per_inf = 1000; % inference will be aborted if fewer present                                         

% max num workers
if savio
    MaxWorkers = 24;
else
    MaxWorkers = 25;
end
%-------------------Load Data and Set Write Paths-------------------------%
datapath = ['../../dat/' project '/']; %Path to inference data
dataname = ['inference_traces_' project '_dT' num2str(Tres) '.mat']; %name of inference set
% Load data for inference into struct named: trace_struct_final
load([datapath dataname]);
load([datapath 'analysis_data_final.mat'])
alpha = trace_struct_final(1).alpha_frac*w; % Rise Time for MS2 Loops in time steps

if set_boostrap && dp_bootstrap
    error('both set and dp bootstrapping set to 1. Choose one.')
end

if set_bootstrap
    d_type = '_set';
elseif dp_bootstrap
    d_type = '_dp';
end

% Set write path (inference results are now written to external directory)
out_suffix =  ['/' project '/w' num2str(w) '_t' num2str(Tres)...
    '_alpha' num2str(round(alpha*10)) '_f' num2str(fluo_field) '_cl' num2str(clipped) ...
    '_no_ends' num2str(clipped_ends) '_off' num2str(off_traces_flag) '/K' num2str(K) ...
    '_tw' num2str(t_window/60) d_type '/']; 
% set write path
if savio
    out_prefix = '/global/scratch/nlammers/'; %hmmm_data/inference_out/';
else    
    out_prefix = '../../out/';
end
out_dir = [out_prefix out_suffix];
mkdir(out_dir);

if clipped_ends % enforce clipped ends 
    end_clip = w + 1;
else
    end_clip = 0;
end

% apply time filtering 
trace_struct_filtered = [];
for i = 1:length(trace_struct_final)
    temp = trace_struct_final(i);
    if clipped
        time = temp.time_interp(w+1:end-end_clip); % we must ignore first w + 1 time points for windowed inference
        if fluo_field == 1            
            fluo = temp.fluo_interp(w+1:end-end_clip); 
        elseif fluo_field == 3           
            fluo = temp.fluo_interp3(w+1:end-end_clip);
        else
            error('unknown fluo field')
        end
    else
        time = temp.time_full;
        if fluo_field == 1            
            fluo = temp.fluo_full;            
        elseif fluo_field == 3            
            fluo = temp.fluo_full3;
        else
            error('unknown fluo field')
        end
        f_ind = find(fluo,1) + w;
        time = time(f_ind:end);        
        fluo = fluo(f_ind:end);
        fluo = fluo(time<=stop_time_inf*60);
        time = time(time<=stop_time_inf*60);
    end
%     tf = time((time>=0)&(time<stop_time_inf*60));
    if length(time) >= min_dp
        temp.fluo = fluo;
        temp.time = time;
        trace_struct_filtered = [trace_struct_filtered temp];
    end
end

trace_struct_filtered = trace_struct_filtered([trace_struct_filtered.inference_flag]==1);
set_index = [trace_struct_filtered.setID];
set_vec = unique(set_index);
if set_bootstrap % override n_bootstrap input if set bootstrapping
    n_bootstrap = length(set_vec);
end

% get last obs time for each set
set_vec = [trace_struct_filtered.setID];
set_index = unique(set_vec);
last_time_vec_set = zeros(1,length(set_index));
for i = 1:length(set_vec)
    set_id = trace_struct_filtered(i).setID;
    lt_curr = last_time_vec_set(set_index==set_id);
    last_time_trace = max(trace_struct_filtered(i).time_interp);
    lt_curr = max(lt_curr,last_time_trace);
    last_time_vec_set(set_index==set_id) = lt_curr;
end

% generate list of all traces for which we observed silencing event
short_list = [];
for i = 1:length(trace_struct_filtered)
    set_id = set_vec(i);
    lt_set = last_time_vec_set(set_index==set_id);
    t_max = max(trace_struct_filtered(i).time_interp);    
    pID = trace_struct_filtered(i).ParticleID;    
    if t_max < (lt_set-300+w*Tres) && ismember(pID,on_off_particles)
        short_list = [short_list i];
    end
end
index_vec = 1:length(trace_struct_filtered);
if off_traces_flag == 1
    trace_struct_filtered = trace_struct_filtered(ismember(index_vec,short_list)); 
elseif off_traces_flag == 2 % back-align
    trace_struct_filtered = trace_struct_filtered(ismember(index_vec,short_list)); 
    for i = 1:length(trace_struct_filtered)
        tt = trace_struct_filtered(i).time;
        tt = tt + 40*60-max(tt); % back-align
        trace_struct_filtered(i).time = tt;
    end
end
% trace start and end times for each trace
first_time_vec = [];
last_time_vec = [];
for i = 1:length(trace_struct_filtered)
    first_time_vec = [first_time_vec min(trace_struct_filtered(i).time)];
    last_time_vec = [last_time_vec max(trace_struct_filtered(i).time)];
end
%% Conduct Inference
% structure array to store the analysis data
local_meta = struct; % Store local_EM results
init_meta = struct; % Store initiation info

for g = 1:length(bin_groups) % loop through different AP groups
    bin_list = bin_groups{g}; % get groups for present iteration         
    for t = 1:length(inference_times)
        t_inf = inference_times(t);
        t_start = t_inf - t_window/2;
        t_stop = t_inf + t_window/2;
        % generate time-resolved AP reference vector
        ap_ref_vec = NaN(1,length(trace_struct_filtered));
        for i = 1:length(trace_struct_filtered)
            tt = trace_struct_filtered(i).time;
            ap = trace_struct_filtered(i).rel_ap_vector_interp;            
            tr_ap = round(mean(ap(tt>=t_start&tt<t_stop))); % use average position         
            ap_ref_vec(i) = tr_ap;
        end        
        for b = 1:n_bootstrap % iterate through bootstrap replicates
            iter_start = now;
            s = (g-1)*n_bootstrap + b;
            local_struct = struct;
            init_struct = struct;
            output = struct;
            % if using set bootstrapping, select set to hold out
            if set_bootstrap
                boot_set = set_vec(b);
            end
            % Use current time as unique inference identifier (kind of unwieldly)
            inference_id = num2str(round(10e5*now));
            if set_bootstrap
                inference_id = [inference_id '_set' num2str(boot_set)];
            end
            % Generate filenames            
            fName_sub = ['eveSet_w' num2str(w) '_K' num2str(K) ...
                '_bin' num2str(round(10*bin_list(1))) '_' num2str(round(10*bin_list(end))) ...
                '_time' num2str(round(t_inf/60)) '_t' inference_id];                
            out_file = [out_dir '/' fName_sub];            
            
            % Extract fluo_data
            trace_filter = ismember(ap_ref_vec,bin_list) & (last_time_vec-w*Tres) >= t_start ...
                    & (first_time_vec+w*Tres) < t_stop;
            if set_bootstrap
                trace_ind = find(([trace_struct_filtered.setID]~=boot_set)&...
                    trace_filter);
            else
                trace_ind = find(trace_filter);
            end
            inference_set = [];
            for m = 1:length(trace_ind)
                temp = trace_struct_filtered(trace_ind(m));
                tt = temp.time;
                ft = temp.fluo;
                temp.time = tt(tt>=t_start & tt < t_stop);
                temp.fluo = ft(tt>=t_start & tt < t_stop);
                if sum(temp.fluo>0) > 1 % exclude strings of pure zeros
                    inference_set = [inference_set temp];
                end
            end
            skip_flag = 0;
            set_size = length([inference_set.fluo]);                 
            if isempty(inference_set)
                skip_flag = 1;
            elseif set_size < min_dp_per_inf                    
                skip_flag = 1;                    
            end
            if skip_flag
                warning('Too few data points. Skipping')                
            else                 
                sample_index = 1:length(inference_set);
                if dp_bootstrap                        
                    ndp = 0;    
                    sample_ids = [];                    
                    %Reset bootstrap size to be on order of set size for small bins
                    if set_size < sample_size
                        sample_size = ceil(set_size/1000)*1000;
                    end
                    while ndp < sample_size
                        tr_id = randsample(sample_index,1);
                        sample_ids = [sample_ids tr_id];
                        ndp = ndp + length(inference_set(tr_id).time);
                    end
                    fluo_data = cell([length(sample_ids), 1]);    
                    time_data = cell([length(sample_ids), 1]);    
                    sample_particles = [inference_set(sample_ids).ParticleID];
                    for tr = 1:length(sample_ids)
                        fluo_data{tr} = inference_set(sample_ids(tr)).fluo;                    
                        time_data{tr} = inference_set(sample_ids(tr)).time;                    
                    end            
                else % Take all relevant traces if not bootstrapping
                    fluo_data = cell([length(inference_set), 1]);            
                    for tr = 1:length(inference_set)
                        fluo_data{tr} = inference_set(tr).fluo;
                        time_data{tr} = inference_set(tr).time;                    
                    end
                end                
                % Random initialization of model parameters
                param_init = initialize_random (K, w, fluo_data);
                % Approximate inference assuming iid data for param initialization                
                local_iid_out = local_em_iid_reduced_memory_truncated (fluo_data, param_init.v, ...
                                    param_init.noise, K, w, alpha, n_steps_max, eps);
                noise_iid = 1/sqrt(exp(local_iid_out.lambda_log));
                v_iid = exp(local_iid_out.v_logs);         
                %%% parpool business
                p = gcp('nocreate');
                if isempty(p)
                    parpool(MaxWorkers); %6 is the number of cores the Garcia lab server can reasonably handle per user.
                elseif p.NumWorkers > MaxWorkers
                    delete(gcp('nocreate')); % if pool with too many workers, delete and restart
                    parpool(MaxWorkers);
                end
                parfor i_local = 1:n_localEM % Parallel Local EM                
                    % Random initialization of model parameters
                    param_init = initialize_random_with_priors(K, noise_iid, v_iid);
                    % Get Intial Values
                    pi0_log_init = log(param_init.pi0);
                    A_log_init = log(param_init.A);
                    v_init = param_init.v;                        
                    noise_init = param_init.noise;
                    % Record
                    init_struct(i_local).A_init = exp(A_log_init);                
                    init_struct(i_local).v_init = v_init;
                    init_struct(i_local).noise_init = noise_init;                
                    init_struct(i_local).subset_id = i_local;
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
%                         local_struct(i_local).total_time = local_out.runtime;
                    local_struct(i_local).total_steps = local_out.n_iter;               
                    local_struct(i_local).soft_struct = local_out.soft_struct;               
                end
                local_meta(s).init = init_struct;
                local_meta(s).local = local_struct;
                [logL, max_index] = max([local_struct.logL]); % Get index of best result                    
                % Save parameters from most likely local run
                output.pi0 =local_struct(max_index).pi0;                        
                output.r = local_struct(max_index).r(:);
                output.off_traces = off_traces_flag;
                output.noise = local_struct(max_index).noise;
                output.A = local_struct(max_index).A(:);
                output.A_mat = local_struct(max_index).A;                            
                % Info about run time
                output.total_steps = local_struct(max_index).total_steps;                                  
                output.total_time = 100000*(now - iter_start);            
                % Save inference ID variables
                output.APbin = min(bin_list):max(bin_list);
                output.boot_set = NaN;
                if set_bootstrap
                    output.boot_set = boot_set;
                else
                    output.boot_set = NaN;
                end
                % other inference characteristics
                output.t_window = t_window;
                output.t_inf = t_inf;
                output.fluo_type = fluo_field;
                output.dp_bootstrap_flag = dp_bootstrap;
                output.set_bootstrap_flag = set_bootstrap;
                output.iter_id = b;
                output.start_time_inf = 0;
                output.stop_time_inf = stop_time_inf;                            
                output.clipped = clipped;
                output.off = off_traces_flag;
                output.clipped_ends = clipped_ends;
                if dp_bootstrap || set_bootstrap                    
                    output.particle_ids = sample_particles;
                    output.N = ndp;
                end                
                output.w = w;
                output.alpha = alpha;
                output.deltaT = Tres; 
                % save inference data used
            end
            output.skip_flag = skip_flag;
            save([out_file '.mat'], 'output');           
        end  
    end
end    
