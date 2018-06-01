% Script to Generate Stripe Profile Arrays and Extract Time On/Off info

% Assumes that Instantaneous fluorescence is a reasonable proxy for
% instantaneous productivity (lagged by w/2). See Appendix 5

addpath('../utilities')
close all
clear 

% ID variables
project = 'mHMMeve2_weka_inf_2018_05_07';
TracePath = ['../../dat/' project '/inference_traces_' project '_dT20.mat'];
load(TracePath)
NucleusPath = ['../../dat/' project '/inference_nuclei_' project '_dT20.mat'];
load(NucleusPath)

DataPath = ['../../dat/' project '/'];
% global param values
n_boots = 1;

%%% ------------------------Stripe Profile Analysis----------------------- %%
TrAPVec = round([trace_struct_final.MeanAP]);
NucAPVec = round([nucleus_struct_final.MeanAP]);
ap_index = -7:7;

InterpGrid = trace_struct_final(1).InterpGrid; % Time grid used for interpolation
TrParticles = [trace_struct_final.ParticleID]; 
TrSets = [trace_struct_final.setID];
NucSets = [nucleus_struct_final.setID];
NucParticles = [nucleus_struct_final.ParticleID];
set_index = unique(TrSets);
nc14_vec = [nucleus_struct_final.on_off_time_flag]; % if 1 nucleus was present for all of nc14
FOV_array = reshape([nucleus_struct_final.in_frame_flags],length(InterpGrid),length(nucleus_struct_final));


%%% Make arrays that aggregate nuclear and spot tracking info
nc_act_array_mf = NaN(size(FOV_array)); % mean field (avg among active)
nc_act_array_full = NaN(size(FOV_array)); % full (avg among all extant nuclei)
trace_mapping_vec = NaN(1,length(nucleus_struct_final));
trace_analysis_vec = zeros(1,length(nucleus_struct_final)); % designates traces that can be used for on/off analyses
for i = 1:length(nucleus_struct_final)
    ParticleID = nucleus_struct_final(i).ParticleID;
    FOV_vec = FOV_array(:,i);
    if isnan(ParticleID)
        nc_act_array_mf(:,i) = NaN; % nuclei without activity are ignored
        nc_act_array_full(FOV_vec==1,i) = 0; % in full case count in mean where present
    else
        temp = trace_struct_final(ParticleID==TrParticles);
        time = temp.time_interp;
        fluo = temp.fluo_interp;
        filter = ismember(InterpGrid,time); % filter for active time points
        nc_act_array_mf(filter,i) = fluo;
        nc_act_array_full(filter,i) = fluo;
        trace_mapping_vec(i) = find(ParticleID==TrParticles,1);
        trace_analysis_vec(i) = temp.on_time_flag & temp.off_time_flag &...
                                nc14_vec(i);
        nc14_vec(i) = nc14_vec(i) & temp.on_time_flag & temp.off_time_flag; % remove nuclei/particles with problem start or stop
    end
end

%%% Calculate average on and off times for each AP region
on_boot_mat = NaN(n_boots,length(ap_index));
off_boot_mat = NaN(n_boots,length(ap_index));
on_time_vec = []; % record individual o times
off_time_vec = []; % record individual off times
on_off_ap_vec = [];
on_off_sets = [];
on_off_particles = [];
for a = 1:length(ap_index)
    ap = ap_index(a);
    % grab relevant traces
    tr_filter = trace_analysis_vec&ismember(trace_mapping_vec,find(TrAPVec==ap));
    ap_traces = nc_act_array_full(:,tr_filter);
    on_off_sets = [on_off_sets NucSets(tr_filter)];
    on_off_particles = [on_off_particles NucParticles(tr_filter)];
    ap_on_times = [];
    ap_off_times = [];
    for p = 1:size(ap_traces,2)
        ap_on_times = [ap_on_times InterpGrid(find(~isnan(ap_traces(:,p)),1))];
        ap_off_times = [ap_off_times InterpGrid(find(~isnan(ap_traces(:,p)),1,'last'))];
    end
    
    % take bootstrap samples
    sample_index = 1:length(ap_on_times);
    for n = 1:n_boots
        boot_samp = randsample(sample_index,length(sample_index),true);
        on_boot_mat(n,a) = mean(ap_on_times(boot_samp));
        off_boot_mat(n,a) = mean(ap_off_times(boot_samp));
    end
    on_time_vec = [on_time_vec ap_on_times];
    off_time_vec = [off_time_vec ap_off_times];
    on_off_ap_vec = [on_off_ap_vec repelem(ap,length(ap_on_times))];    
end

% results_struct = struct;
mean_ap_on = nanmean(on_boot_mat);
ste_ap_on = nanstd(on_boot_mat);
mean_ap_off = nanmean(off_boot_mat);
ste_ap_off = nanstd(off_boot_mat);

%%% Perform bootstrap sampling of approximate nuclear mRNA level for Each AP region
% n_boots = 10;
lambda = 1/(7*60); % decay rate (7 min)
decay_kernel = fliplr(exp(-lambda*InterpGrid)); % for decay simulation
fraction_on = NaN(length(InterpGrid),length(ap_index),n_boots);
full_mRNA_decay = NaN(length(InterpGrid),length(ap_index),n_boots); % proxy for observed scenario
full_mRNA_no_decay = NaN(length(InterpGrid),length(ap_index),n_boots); % effect of decay alone
on_mRNA_decay = NaN(length(InterpGrid),length(ap_index),n_boots); % everything but fraction on
on_mRNA_no_decay = NaN(length(InterpGrid),length(ap_index),n_boots); % less decay
ss_bin_mRNA_decay = NaN(length(InterpGrid),length(ap_index),n_boots); % combines ss and switching
ss_bin_mRNA_no_decay = NaN(length(InterpGrid),length(ap_index),n_boots); % combines ss and switching
binary_mRNA_decay = NaN(length(InterpGrid),length(ap_index),n_boots); % effect of switching + decay
binary_mRNA_no_decay = NaN(length(InterpGrid),length(ap_index),n_boots); % effecto of switching alone
rate_mRNA = NaN(length(InterpGrid),length(ap_index),n_boots); % effect of switching + decay
ss_mRNA_rate = NaN(n_boots,length(ap_index)); % effect of switching alone
f_on_mat = NaN(n_boots,length(ap_index));
% ref names
data_names = {'full_mRNA_decay', 'full_mRNA_no_decay',...
    'ss_bin_mRNA_decay','ss_bin_mRNA_no_decay','rate_mRNA',...
    'binary_mRNA_decay','binary_mRNA_no_decay','on_mRNA_decay',...
    'ss_bin_mRNA_decay', 'ss_bin_mRNA_no_decay', 'on_mRNA_no_decay'};
%%% iterate through ap positions, conduct bootstrap sampling
for a = 1:length(ap_index)
    ap = ap_index(a);
    % filtered arrays
    ap_nuclei_mf = nc_act_array_mf(:,nc14_vec&NucAPVec==ap); % used for mf analyses 
    ap_nuclei_full = nc_act_array_full(:,nc14_vec&NucAPVec==ap); % used for full analyses
    samp_vec = 1:size(ap_nuclei_full,2); % convenience vector for sampling
    for n = 1:n_boots
        % take bootstrap samples
        mf_sample = randsample(samp_vec,length(samp_vec),true);
        full_sample = randsample(samp_vec,length(samp_vec),true);        
        boot_mf = ap_nuclei_mf(:,mf_sample);
        boot_full = ap_nuclei_full(:,full_sample);        
        active_boot_full = boot_full(:,max(boot_full)>0);
        % model binary switch behavior
        on_off_boot_mat = zeros(length(InterpGrid),size(active_boot_full,2)); % binary matrix for switching analyses
        for i = 1:size(active_boot_full,2)
            on_off_boot_mat(find(active_boot_full(:,i)>0,1):find...
                (active_boot_full(:,i)>0,1,'last'),i) = 1;
        end
        n_steps = sum(on_off_boot_mat(:));
        ss_mRNA_rate(n,a) = nanmean(boot_mf(:)); % equal weight to each time point (not necessarily each nucleus)
        f_on_mat(n,a) = sum(max(boot_full)>0)/size(boot_full,2); % fraction turned on ever
        % now move through time
        for i = 1:length(InterpGrid)
            t = InterpGrid(i);
            fraction_on(i,a,n) = sum(boot_full(i,:)>0)/size(boot_full,2);            
            ss_wt_vec = nansum(active_boot_full(1:i,:),2)./nansum(on_off_boot_mat(1:i,:),2);
            rate_mRNA(i,a,n) = nanmean(ss_wt_vec);
            % no decay (just take sum)            
            binary_mRNA_no_decay(i,a,n) = nansum(mean(on_off_boot_mat(1:i,:),2));
            ss_bin_mRNA_no_decay(i,a,n) = nansum(nanmean(on_off_boot_mat(1:i,:).*active_boot_full(1:i,:),2));                
            full_mRNA_no_decay(i,a,n) = nansum(nanmean(boot_full(1:i,:),2));
            on_mRNA_no_decay(i,a,n) = nansum(nanmean(active_boot_full(1:i,:),2));
            % for decay case we must weight past time steps with decay
            % kernel
            binary_mRNA_decay(i,a,n) = nansum(decay_kernel(end-i+1:end)'.*nanmean(on_off_boot_mat(1:i,:),2));  
            ss_bin_mRNA_decay(i,a,n) = nansum(decay_kernel(end-i+1:end)'.*nanmean(on_off_boot_mat(1:i,:).*active_boot_full(1:i,:),2));  
            full_mRNA_decay(i,a,n) = nansum(decay_kernel(end-i+1:end)'.*nanmean(boot_full(1:i,:),2));
            on_mRNA_decay(i,a,n) = nansum(decay_kernel(end-i+1:end)'.*nanmean(active_boot_full(1:i,:),2));
        end        
    end
end

results_struct = struct;
results_struct.mean_ss_rate = nanmean(ss_mRNA_rate)/sum(nanmean(ss_mRNA_rate));
results_struct.ste_ss_rate = nanstd(ss_mRNA_rate)/sum(nanmean(ss_mRNA_rate));
results_struct.mean_fraction_on = nanmean(f_on_mat)/sum(nanmean(f_on_mat));
results_struct.ste_fraction_on = nanstd(f_on_mat)/sum(nanmean(f_on_mat));

for i = 1:length(data_names)
    data = eval(data_names{i});        
    results_struct.(['mean_' data_names{i}]) = nanmean(data,3)./repmat(sum(nanmean(data,3),2),1,length(ap_index));
    results_struct.(['ste_' data_names{i}]) = nanstd(data,[],3)./repmat(sum(nanmean(data,3),2),1,length(ap_index));                                                    
end
clear trace_struct_final
clear nucleus_struct_final
save([DataPath 'analysis_data_final.mat'])
