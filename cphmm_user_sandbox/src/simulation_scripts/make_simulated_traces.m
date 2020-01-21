% This script uses the Gillespie Algorithm to generate a set of simulated
% MS2 traces that can be to test the cpHMM inference pipeline
% Author: Nick Lammers

% clear workspace
clear
close all

% add path to utilities folder (contains useful functions)
addpath('../utilities/');

project = 'example_trace_set';
outpath = ['../../dat/' project '/'];
mkdir(outpath);
% Set simulation parameters
K = 2; % number of promoter states
w = 6; % system memory
Tres = 20; % time resolution
t_MS2 = 30; % time it takes to transcribe the MS2 loop [sec]
alpha = t_MS2./ Tres ; % alpha: length of the MS2 loop in time steps

% group-specific parameters
inference_groups = [1 2];
inference_labels = {'test1', 'test2'};

% define group-specific parameteres
R_vec = {[-0.0150,  0.060; 0.0150, -0.060],.5*[-0.0150,  0.060; 0.0150, -0.060]}; % transition rates
pi0_vec = {[0.3, 0.7],[0.3, 0.7]}; % initial state pmf
noise_vec = [100, 150]; % background noise [a.u.]
r_emission_vec = {[0, 60, 130],[0, 60, 130]*1.5}; % emission rate [a.u. / sec]

%%% Analysis parameters 
n_points_total = 10000; % total number of time points in a pooled data set
trace_duration = 2400; % time of the process [sec]
seq_length = round(trace_duration/Tres);
n_traces = round(n_points_total / seq_length);

% structure to store traces and synthetic parameter values
synthetic_parameters = struct;
synthetic_parameters.R = R_vec;
synthetic_parameters.inference_groups = inference_groups;
synthetic_parameters.inference_labels = inference_labels;
synthetic_parameters.noise = noise_vec;
synthetic_parameters.r_emission = r_emission_vec;
synthetic_parameters.process_time = trace_duration;
synthetic_parameters.K = K;
synthetic_parameters.w = w;
synthetic_parameters.alpha_frac = t_MS2 / Tres / w;
synthetic_parameters.t_MS2 = t_MS2;
synthetic_parameters.Tres = Tres;
synthetic_parameters.pi0 = pi0_vec;
synthetic_parameters.seq_length = seq_length;
synthetic_parameters.n_traces = n_traces;
synthetic_parameters.n_points_total = n_points_total;

%%% ---------------------- Generate synthetic data ----------------------
trace_struct_final = struct;
i_pass = 1;
for i = 1:length(inference_groups)
    i_group = inference_groups(i);
    i_label = inference_labels{i};
    for tr = 1:n_traces
        fluo_gill = synthetic_rate_gillespie(seq_length, alpha, ...
            K, w, R_vec{i}, Tres, r_emission_vec{i}, noise_vec(i), pi0_vec{i});

        trace_struct_final(i_pass).Tres = Tres;
        trace_struct_final(i_pass).fluo_interp = fluo_gill.fluo_MS2; % simulated fluorescence
        trace_struct_final(i_pass).time_interp = Tres:Tres:Tres*seq_length; % simulation time
        trace_struct_final(i_pass).promoter_trajectory = fluo_gill.naive_states;
        trace_struct_final(i_pass).promoter_transition_times = fluo_gill.transition_times;     
        trace_struct_final(i_pass).alpha_frac = synthetic_parameters.alpha_frac; % length of MS2 loops (as fraction of total gene length)
        trace_struct_final(i_pass).group_id = i_group; 
        trace_struct_final(i_pass).group_label = i_label; % group membership
        trace_struct_final(i_pass).ParticleID = i_pass; % unique particle identifier
        i_pass = i_pass + 1;
    end  
end
save([outpath '/trace_struct_final.mat'],'trace_struct_final')
save([outpath '/synthetic_parameters.mat'],'synthetic_parameters')