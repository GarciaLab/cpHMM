% Script to clean and interpolate fluorescence trace
% Also clean up nucleus data-set
clear 
close all
%%%Set Cleaning and Summary Parameters
min_dp = 10; % minimum # dp acceptable for an inference trace
Tres_interp = 20; % time resolution
InterpGrid = 0:Tres_interp:60*50;
min_bin_size = 1;
max_bin_size = 3;
FOV_edge_padding = 10; % pixels
%------------------------Import Raw Trace Set------------------------%
%ID's of sets to include
include_vec = [9,10,11,12,13,19,20,22,23,24,26];
project = 'mHMMeve2_weka_inf_2018_05_07';
print_traces = 1; %Output PNG files for traces?

%---------------------------Set Paths-------------------------------------%
TracePath = ['../../dat/' project '/' 'raw_traces_01' project '.mat'];
NucleusPath = ['../../dat/' project '/' 'ellipse_info_01' project '.mat'];
OutPath = ['../../dat/' project '/'];
FigPath = ['../../fig/experimental_system/' project '/preprocessing'];
TraceSavePath = [FigPath '/traces/'];
mkdir(TraceSavePath);
mkdir(OutPath);
% Save Names
CleanTraceName = ['inference_traces_' project '_dT' num2str(Tres_interp) '.mat'];
CleanNucleusName = ['inference_nuclei_' project '_dT' num2str(Tres_interp) '.mat'];

%----------------Load Traces and Perform First Pass QC--------------------%
%Load raw traces (saved in struct titled "trace_struct")
load(TracePath);
trace_struct = new_trace_struct(ismember([new_trace_struct.setID],include_vec));
load(NucleusPath);
nucleus_struct = new_nucleus_struct(ismember([new_nucleus_struct.setID],include_vec));

%%% Cleaning Parameters
big_jump1 = prctile([trace_struct.fluo],99);
jump_threshold1 = big_jump1/1.5; % this should be a conservative threshold for single time step increase
big_jump3 = prctile([trace_struct.fluo3],99);
jump_threshold3 = big_jump3/1.5; % this should be a conservative threshold for single time step increase
index_vec = 1:length(trace_struct); % convenience ref vector
field_names = fieldnames(trace_struct);

trace_struct_final = []; % Structure to store traces after first round of cleaning
jump_ct = 0;
blip_ct1 = 0;
blip_ct3 = 0;
for i = 1:length(trace_struct) 
    temp = trace_struct(i);
    trace1 = temp.fluo; %Load full trace, including intervening NaN's
    trace3 = temp.fluo3; %3 slice trace should be identical to single in presence/absence
    time = temp.time;      
    quality_flag = 1;
    if sum(~isnan(trace1)) < min_dp
        quality_flag = 0;
    end
    %Null assumption is that all clusters of 6 or more NaNs are 0s. Single
    %,double, or triple NaNs are assumed to have been missed nonzero dps
    trace1_nans = isnan(trace1);      
    %Look for clusters of 6 or more NaNs
    kernel = [1,1,1,1,1];
    tn_conv = conv(kernel,trace1_nans);
    tn_conv = tn_conv(3:end-2);
    z_ids = find(tn_conv==5);
    z_ids = unique([z_ids-1 z_ids z_ids+1]); % get set of z_ids    
    trace1(z_ids) = 0; % set clusters to zeros
    trace3(z_ids) = 0;
    trace1(trace1<0) = 0; % deal with negative values
    trace3(trace3<0) = 0;
    % find single dp "blips"
    tr_dd1 = abs([0 diff(diff(trace1)) 0]);  % 1 slice
    trace1(tr_dd1>2*jump_threshold1) = NaN;
    blip_ct1 = blip_ct1 + sum(tr_dd1>2*jump_threshold1);
    tr_dd3 = abs([0 diff(diff(trace3)) 0]);  % 3 slice  
    trace3(tr_dd3>2*jump_threshold3) = NaN;
    blip_ct3 = blip_ct3 + sum(tr_dd3>2*jump_threshold3);
    
    % interpolate remaining NaNs    
    query_points1 = time(isnan(trace1));%InterpGrid((InterpGrid>=min(time))&(InterpGrid<=max(time)));
    interp_t1 = time(~isnan(trace1));
    interp_f1 = trace1(~isnan(trace1));
    new_f1 = interp1(interp_t1,interp_f1,query_points1);  
    trace1(ismember(time,query_points1)) = new_f1;
    
    query_points3 = time(isnan(trace3));
    interp_t3 = time(~isnan(trace3));
    interp_f3 = trace3(~isnan(trace3));
    new_f3 = interp1(interp_t3,interp_f3,query_points3);  
    trace3(ismember(time,query_points3)) = new_f3;
    
    %%% flag traces with unreasonably large rises or falls    
    tr_d1 = diff(trace1);
    tr_d3 = diff(trace3);
    if max(abs(tr_d1)) >= jump_threshold1 || max(abs(tr_d3)) >= jump_threshold3
        jump_ct = jump_ct + 1;
        quality_flag = 0;
    end
    
    % Interpolate to standardize spacing
    t_start = InterpGrid(find(InterpGrid>=time(1),1));
    t_stop = InterpGrid(find(InterpGrid<=time(end),1,'last'));
    time_interp = t_start:Tres_interp:t_stop;
    trace1_interp = interp1(time,trace1,time_interp);
    trace3_interp = interp1(time,trace3,time_interp);
    interp_fields = {'xPos','yPos','ap_vector','rel_ap_vector'};
    % interpolate other vector fields
    for j = 1:length(interp_fields)
        int_vec = temp.(interp_fields{j});
        int_time = temp.time;
        cp_frames = temp.cp_frames;
        all_frames = temp.all_frames;
        int_time = int_time(ismember(all_frames,cp_frames));
        temp.([interp_fields{j} '_interp']) = interp1(int_time,int_vec,time_interp);
    end     
    temp.fluo_interp = trace1_interp;
    temp.fluo_interp3 = trace3_interp;
    temp.time_interp = time_interp;
    temp.inference_flag = quality_flag;
    trace_struct_final = [trace_struct_final temp];    
end

% Find last obs times for each set
last_times = [];
for s = 1:length(include_vec)
    setID = include_vec(s);
    last_times = [last_times max([trace_struct_final([trace_struct_final.setID]...
        ==setID).time_interp])];
end
time_ceiling = InterpGrid(find(InterpGrid<=min(last_times),1,'last'));
%%% Add a few useful fields
for i = 1:length(trace_struct_final)
    fluo_interp = trace_struct_final(i).fluo_interp;
    fluo_interp3 = trace_struct_final(i).fluo_interp3;
    time_interp = trace_struct_final(i).time_interp;    
    setID = trace_struct_final(i).setID;        
    time_full = InterpGrid;
    fluo_full1 = zeros(1,length(time_full));
    fluo_full3 = zeros(1,length(time_full));
    fluo_full1(ismember(time_full,time_interp)) = fluo_interp;
    fluo_full3(ismember(time_full,time_interp)) = fluo_interp3;
    trace_struct_final(i).fluo_full = fluo_full1; %"unclipped" version
    trace_struct_final(i).fluo_full3 = fluo_full3; %"unclipped" version
    trace_struct_final(i).time_full = time_full;    
    trace_struct_final(i).N = length(fluo_interp);
    trace_struct_final(i).dT = Tres_interp;                
    trace_struct_final(i).alpha_frac = 1302/6544;
    trace_struct_final(i).InterpGrid = InterpGrid;
    trace_struct_final(i).FOV_edge_padding = FOV_edge_padding;
    trace_struct_final(i).MeanAPOrig = mean(trace_struct_final(i).ap_vector);   
end

%------------------------- Clean Ellipse Set -----------------------------%

% padding = 10; % # frames at end or beginning of nc14 during which 
trace_particle_vec = [trace_struct_final.ParticleID];

nucleus_struct_final = [];
rm_counts = 0;
pixel_per_ap_vec = unique([trace_struct_final.PixelperAP]);
set_vec = unique([trace_struct_final.setID]);
for i = 1:length(nucleus_struct)
    temp = nucleus_struct(i);
    setID = temp.setID;
    lt = last_times(include_vec==setID);
    nc_times = temp.time;         
    temp.InterpGrid = InterpGrid;
    % interpolate relevant fields
    time = temp.time;    
    dT = diff(time);
    temp.quality_flag = 1;
    if max(dT) > 120 || length(time) < min_dp 
        temp.quality_flag = 0;
    end
    %%% assign on and off time flags
    temp.on_time_flag = min(time) <= 6*60 && quality_flag;
    temp.off_time_flag = max(time) >= time_ceiling - 120 && quality_flag;
    
    frames = temp.frames;    
    t_start = InterpGrid(find(InterpGrid>=time(1),1));    
    t_stop = InterpGrid(find(InterpGrid<=time(end),1,'last'));
    time_interp = t_start:Tres_interp:t_stop;
    frames_interp = linspace(min(frames), max(frames), length(time_interp));
    tracking_flags = ismember(floor(frames_interp),frames);
    interp_fields = {'time','xPos','yPos', 'rel_ap_vector'};
    for j = 1:length(interp_fields)
        int_vec = temp.(interp_fields{j});
        int_time = temp.time;
        interp = interp1(int_time,int_vec,time_interp);
        interp(~tracking_flags) = NaN;
        temp.([interp_fields{j} '_interp']) = interp;
    end
    % check for traking issues--step-over-step disp > .25 AP
    nt = temp.time_interp;
    n_ap = temp.rel_ap_vector_interp;
    disp_nc = max(abs(diff(n_ap(nt>600))));
    if isempty(disp_nc)
        ap_shift_flag = 1;
    else
        ap_shift_flag = disp_nc < .25;
    end
    temp.ap_shift_flag = ap_shift_flag;
    temp.on_off_time_flag = temp.on_time_flag&temp.off_time_flag&ap_shift_flag;    
    FOV_flags = (temp.xPos_interp > FOV_edge_padding) & (temp.xPos_interp < 512- FOV_edge_padding)...
        &(temp.yPos_interp > FOV_edge_padding) & (temp.yPos_interp < 256 - FOV_edge_padding);    
    edge_flags_all = NaN(size(InterpGrid));
    edge_flags_all(ismember(InterpGrid,time_interp)) = FOV_flags;
    temp.FOV_flags = edge_flags_all;
    temp.in_frame_flags = ~isnan(edge_flags_all);
    temp.FOV_edge_padding = FOV_edge_padding;
    % assign trace indices
    ParticleID = temp.ParticleID;
    trace_ind = [];
    if ~isnan(ParticleID)
        trace_ind = find(trace_particle_vec==ParticleID);
        if isempty(trace_ind) % catch cases in which particle has been removed for QC reasons
            error('Inconsistent indexing between Trace and Nucelus Structures')
            trace_ind = NaN;
        else
            ncIDCheck = trace_struct_final(trace_ind).ncID;
            if ncIDCheck ~= temp.ncID
                error('Inconsistent Particle and NC Identifiers')
            end
        end
    else
        trace_ind = NaN;
    end
    if length(trace_ind) > 1
        error('Degenerate Identifiers')
    end
    temp.TraceIndex = trace_ind;    
    nucleus_struct_final = [nucleus_struct_final temp];
end

%%%%%%%%%%%%%%%%%%%%% Trace and Nucleus On/Off Screen
%%% Screen for traces for which start and/or stop times may not be reliable
%%% Most common causes of this are FOV edge effects
%%% Check for Particles that: 
                % 1) start or end ``high'' 
                % 2) or for which first or last particle obs == first or
                %    last nucleus obs
                                

% trace index vec
nucleus_trace_vec = [nucleus_struct_final.TraceIndex];
mismatches = 0;
for i = 1:length(trace_struct_final)
    lt = last_times(include_vec==trace_struct_final(i).setID);
    ParticleID = trace_struct_final(i).ParticleID;
    ncID = find(nucleus_trace_vec==i);
    pIDCross = nucleus_struct_final(ncID).ParticleID;
    if pIDCross ~= ParticleID
        error('Inconsistent Identifiers')
    end   
%     trace_struct_final(i).inference_flag = trace_struct_final(i).inference_flag && nc_quality;
    fluo = trace_struct_final(i).fluo_interp;
    time = trace_struct_final(i).time_interp;
    ap_vec = trace_struct_final(i).rel_ap_vector_interp;
    ap_shift_flag = abs(sum(diff(trace_struct_final(i).rel_ap_vector))) < 4;
    % start time screens
    start_jump_flag = fluo(1) < jump_threshold1; % traces that start bright are suspect
    start_nc_flag = time(1) ~= nucleus_struct_final(ncID).time_interp(1);
    start_early_flag = time(1) < 15*60;
    trace_struct_final(i).on_time_flag = start_jump_flag&&start_nc_flag&&start_early_flag;
    % stop time screens
    stop_jump_flag = ~(fluo(end) > jump_threshold1 && time(end) < time_ceiling); % only problematic when both are true
    stop_nc_flag = ~(time(end) == nucleus_struct_final(ncID).time_interp(end)&& ...
                        time(end) < time_ceiling);
    trace_struct_final(i).off_time_flag = stop_jump_flag&&stop_nc_flag;
    trace_struct_final(i).ap_shift_flag = ap_shift_flag;
    trace_struct_final(i).on_off_time_flag = start_jump_flag&&stop_jump_flag&&...
                            start_nc_flag&&stop_nc_flag&&ap_shift_flag;
end


% save
save([OutPath CleanTraceName],'trace_struct_final');
save([OutPath CleanNucleusName],'nucleus_struct_final');
%%
cm = jet(128);
%%% Make N DP Histogram
ap_vec = round([trace_struct_final.MeanAP]);
inf_vec = [trace_struct_final.inference_flag];
ap_range = -7:7;
n_dp = NaN(1,length(ap_range));
for i = 1:length(ap_range)
    n_dp(i) = length([trace_struct_final(inf_vec==1&ap_vec==ap_range(i)).fluo_interp]);
end
dp_fig = figure;
hold on
plot(ap_range,n_dp,'Color','black')
scatter(ap_range,n_dp,'MarkerFaceColor',cm(30,:),'MarkerEdgeColor','black')
xlabel('relative AP position')
ylabel('N dp')
saveas(dp_fig, [FigPath, 'n_data_points_ap.png'],'png');
%%
%%% If desired, save individual trace plots
if print_traces
    for i = 1:length(trace_struct_final)
        on_flag =  trace_struct_final(i).on_time_flag;            
        off_flag =  trace_struct_final(i).off_time_flag;            
        ParticleID = trace_struct_final(i).ParticleID;
        t_fig = figure('Visible','off');
        cm = jet(64);
        hold on
        plot(trace_struct_final(i).time / 60, trace_struct_final(i).fluo...
                , '-o','Linewidth',1.5,'Color',cm(15,:)/1.5)
%         plot(trace_struct_final(i).time / 60, trace_struct_final(i).fluo3...
%                 , '-','Linewidth',1.5,'Color',cm(30,:)/1.5)
        plot(trace_struct_final(i).time_interp / 60, trace_struct_final(i).fluo_interp,...
            '-s', 'Linewidth',1.5,'Color',cm(30,:))
%         plot(trace_struct_final(i).time_interp / 60, trace_struct_final(i).fluo_interp3,...
%             '-s', 'Linewidth',1.5,'Color',cm(30,:))
        inf_flag = trace_struct_final(i).inference_flag;
        title(['Original vs. Final: Trace ' num2str(ParticleID) ' (inf flag=' num2str(inf_flag) ...
            ', off flag=' num2str(off_flag) ', on flag='  num2str(on_flag) ')'])            
        legend('Raw', 'Interpolated');
        xlabel('Minutes');
        grid on        
        saveas(t_fig, [TraceSavePath, 'trace_' num2str(ParticleID) '.png'],'png');
        close all
    end
end