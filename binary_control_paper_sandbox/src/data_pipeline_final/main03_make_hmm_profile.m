% Script to generate "Analogue" Profile Using Inference Results
% Will use posterior decoded matrices to improve spatial resolution to 1 AP
addpath('../utilities/');
close all
clear 
%------------------------------Set System Params--------------------------%
w = 7; %memory assumed for inference
K = 3; %states used for final inference
Tres = 20; %Time Resolution
alpha = 1.4; % MS2 rise time in time steps
fluo_field = 1; % type of spot integration used
clipped_traces = 1; % if 0, traces are taken to be full length of nc14
t_window = 50*60;
no_ends_flag = 1;
off_traces_flag = 0;
%-----------------------------ID Variables--------------------------------%
% write_csv  = 1; % if 1 write to csv (only for ss)
% id variables
datatype = 'weka';
d_type = '_dp';
project = 'mHMMeve2_weka_inf_2018_05_07'; %project identifier
suffix = '';
%Generate filenames and writepath
id_string =  ['/w' num2str(w) '_t' num2str(Tres)...
    '_alpha' num2str(round(alpha*10)) '_f' num2str(fluo_field) '_cl' num2str(clipped_traces) ...
    '_no_ends' num2str(no_ends_flag) '_off' num2str(off_traces_flag) '/K' num2str(K) ...
    '_tw' num2str(t_window/60) d_type suffix '/']; 
OutPath = ['../../dat/' project '/' id_string];
FigPath = ['../../fig/experimental_system/' project '/' id_string];
mkdir(OutPath)
mkdir(FigPath)
%%% path to raw data
% DropboxFolder = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)\hmmm_data\inference_out\';
DropboxFolder = 'E:/Nick/Dropbox (Garcia Lab)/hmmm_data/inference_out/';
folder_path =  [DropboxFolder '/' project '/' id_string];

%---------------------------------Read in Files---------------------------%
files = dir(folder_path);
filenames = {};
for i = 1:length(files)
    if ~isempty(strfind(files(i).name,['w' num2str(w)])) && ...
       ~isempty(strfind(files(i).name,['K' num2str(K)]))
        filenames = [filenames {files(i).name}];
    end
end

if isempty(filenames)
    error('No file with specified inference parameters found')
end


%Iterate through result sets and concatenate into 1 combined struct
glb_all = struct;
f_pass = 1;
% p_mat = NaN(length(filenames),20);
% dup_list = [];
for f = 1:length(filenames)
    % load the eve validation results into a structure array 'output'    
    load([folder_path filenames{f}]);
%     if output.skip_flag == 1
    if length(fieldnames(output)) < 2
        continue
    end
%     for i = 1:length(filenames)
%         if sum(ismember(output.particle_ids(1:20),p_mat(i,:)))==20
%             dup_list = [dup_list f];
%         end
%     end
%     p_mat(f,:) = output.particle_ids(1:20);
    for fn = fieldnames(output)'
        glb_all(f_pass).(fn{1}) = output.(fn{1});
    end
    glb_all(f_pass).source = filenames{f};        
    f_pass = f_pass + 1
end
%%
%%% Load inference traces
datapath = ['../../dat/' project '/']; %Path to raw data
% generate read and write names
dataname = ['inference_traces_' project '_dT' num2str(Tres) '.mat'];
% Load data for inference into struct named: trace_struct_final
load([datapath dataname]);

%%%------------------Define Indexing Vectors----------------------------%%%
tr_particle_vec = [trace_struct_final.ParticleID];
tr_ap_vec = round([trace_struct_final.MeanAP]);
ap_index = -7:7;
fluo_cell = cell(1,length(ap_index));
state_cell = cell(1,length(ap_index));
f_ct_cell = cell(1,length(ap_index));
% Make alpha correction kernel
alpha_kernel = [];
times = 0:Tres:w*Tres;
for n = 1:w
    t1 = times(n);
    t2 = times(n+1);
    alpha_kernel = [alpha_kernel ms2_loading_coeff_integral(alpha, w, Tres, t1, t2)];
end
alpha_kernel = fliplr(alpha_kernel)/Tres;
alpha_mat = repmat(alpha_kernel,K,1)';
iter = 0;
for i = 1:length(glb_all)
    particle_list = glb_all(i).particle_ids;    
    for p = 1:length(particle_list)
        ap = tr_ap_vec(tr_particle_vec==particle_list(p));
        ap_vector = round(trace_struct_final(tr_particle_vec==particle_list(p)).rel_ap_vector_interp);
        ap_vector = ap_vector(w+1:end-w-1);       
        if sum(ismember(ap_vector,ap_index))==0
            continue
        end
        fluo = glb_all(i).fluo_data{p}(w+1:end);        
        promoter_path = exp(glb_all(i).soft_struct.p_z_log_soft{p})';
        f_ct = zeros(size(promoter_path,1)-w,K);
        for j = w+1:size(promoter_path,1)
            f_ct(j-w,:) = sum(promoter_path(max(1,j-w+1):j,:).*alpha_mat(max(1,w-j+1):end,:));
        end     
        ap_unique = unique(ap);
        for a = 1:length(ap_unique)
            ap = ap_unique(a);
            if sum(ap==ap_index)==0
                continue
            end
            fc = f_ct_cell{ap==ap_index};
            pp = state_cell{ap==ap_index};
            ff = fluo_cell{ap==ap_index};
            if isempty(pp)
                pp = promoter_path(ap_vector==ap,:);
                ff = fluo(ap_vector(w+1:end)==ap)';
                fc = f_ct(ap_vector(w+1:end)==ap,:);
            else
                pp = [pp ; promoter_path(ap_vector==ap,:)];
                ff = [ff ; fluo(ap_vector(w+1:end)==ap)'];
                fc = [fc ; f_ct(ap_vector(w+1:end)==ap,:)];
            end
            state_cell{ap==ap_index} = pp;
            fluo_cell{ap==ap_index} = ff;        
            f_ct_cell{ap==ap_index} = fc;        
        end        
    end
    iter = iter + 1
end
%% Calculate emission values and occupancies
r_mat = NaN(length(ap_index),K);
occ_mat = NaN(length(ap_index),K);
for a = 1:length(ap_index)
    fluo_vec = fluo_cell{a};
    state_mat = state_cell{a};
    state_count_mat = f_ct_cell{a};
    F_square = zeros(K,K);
    b_agg = zeros(1,K);    
    for k = 1:K
        F_square(k,:) = sum(state_count_mat.*repmat(state_count_mat(:,k),1,K));
        b_agg(k) = sum(fluo_vec.*state_count_mat(:,k));
    end
    % solve system                
    r_mat(a,:) = linsolve(F_square,b_agg')';     
    occ_mat(a,:) = sum(state_mat)/sum(state_mat(:));
end
production_vec = sum(r_mat.*occ_mat,2);
hmm_profile.production_vec = production_vec;
hmm_profile.r_mat = r_mat;
hmm_profile.occ_mat = occ_mat;
save([OutPath 'mean_profile_hmm.mat'],'hmm_profile')