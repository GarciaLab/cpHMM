% Script to Compile Data Sets and Find Stripe Centers
close all
clear 
%--------------------------Set Path Specs, ID Vars------------------------%
FolderPath = 'E:/Nick/Dropbox (Garcia Lab)/mHMM/weka/';
% FolderPath = 'D:/Data/Nick/LivemRNA/LivemRNAFISH/Dropbox (Garcia Lab)/mHMM/weka/';
project = 'mHMMeve2_weka_inf_2018_05_07'; %Project Identifier
%%% folders
fig_path = ['../../fig/experimental_system/' project '/preprocessing/'];
data_path = ['../../dat/' project '/']; % data mat directory
%%% fig subfolders
ap_pos_path = [fig_path 'ap_positioning/'];
fluo_path = [fig_path 'fluo_stats/'];
%%% assign save names
trace_name = [data_path 'raw_traces_01' project]; % names for compiled trace struct
nucleus_name = [data_path 'ellipse_info_01' project]; % names for compiled elipse struct

%%% make filepaths
mkdir(data_path);
mkdir([ap_pos_path '/stripe_fits']);
mkdir(fluo_path);
%%% cleaning params
keyword = '20sec'; % Keyword to ensure only sets from current project are pulled
include_vec = [9:13 19 20 22:24 26]; %data set numbers to include
% show_ap_fit_figs = 0;
snippet_size = 15; % particles within snippet/2+1 are at risk for tracking issues
% pre_post_padding = 10; % max mun frames for which nucleus can be MIA at start or end
%--------------------------Obtain Relevant Filepaths----------------------%
% store set names
dirinfo = dir(FolderPath);
dirinfo(~[dirinfo.isdir]) = []; %remove non-directories
cp_filenames = {}; % particles
ap_filenames = {}; % ap info
nc_filenames = {}; % nuclei
fov_filenames = {}; % fov info
set_nums = [];
for d = 1 : length(dirinfo)
    thisdir = dirinfo(d).name;
    % Skip files lacking project keyword 
    if isempty(strfind(thisdir,keyword)) 
        continue
    end
    % Remove sets not in include_vec
    set_num_start_ind = strfind(thisdir,'_');
    set_num_start_ind = set_num_start_ind(end);
    set_num = str2num(thisdir(set_num_start_ind+1:end));    
    if sum(set_num==include_vec) ~= 1 
        continue
    end    
    set_nums = [set_nums set_num];
    % append file paths
    cp_filenames = [cp_filenames {[thisdir '/CompiledParticles.mat']}];    
    ap_filenames = [ap_filenames {[thisdir '/APDetection.mat']}];    
    nc_filenames = [nc_filenames {[thisdir '/' thisdir '_lin.mat']}];           
    fov_filenames = [fov_filenames {[thisdir '/FrameInfo.mat']}];           
end

trace_struct = struct; % Generate data structure to store extracted trace sets
nucleus_struct = []; % structure to store nucleis info

%%% compile traces and nuclei from across experimental sets
j_pass = 0; % indexing variable to keep track of iterations
total_matched = 0;
for i = 1:length(cp_filenames) % Loop through filenames    
    % read in raw files
    load([FolderPath ap_filenames{i}]) % AP Info   
    load([FolderPath nc_filenames{i}]) % Ellipse Info
    load([FolderPath fov_filenames{i}]) % FrameInfo Info
    yDim = FrameInfo(1).LinesPerFrame;
    xDim = FrameInfo(1).PixelsPerLine;
    raw_data = load([FolderPath cp_filenames{i}]); % Particles    
    setID = set_nums(i);    
    % get angle between the x-axis and the AP-axis 
    APAngle = round(atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)))*360 / (2*pi));    
    %Correction for if APAngle is in quadrants II or III
    if coordPZoom(1)-coordAZoom(1) < 0       
        APAngle = APAngle + 180;
    end    
    % pull trace and nuclei variables
    time_raw = raw_data.ElapsedTime*60; % time vector            
    traces_raw = raw_data.AllTracesVector; % Array with a column for each trace    
    frames_raw = 1:length(time_raw); % Frame list    
    first_frame = raw_data.nc14; % Get frame that marks start of nc14
    last_frame = frames_raw(end); % last frame in set
    % filter trace mat and time
    traces_clean = traces_raw(first_frame:end,:);
    time_clean = time_raw(first_frame:end);    
    time_clean = time_clean - min(time_clean); % Normalize to start of nc14
    frames_clean = frames_raw(first_frame:end);    
    % compile nucleus info
    s_cells = struct;
    e_pass = 1;    
    for e = 1:length(schnitzcells)
        e_frames = schnitzcells(e).frames;
        nc14_frames = e_frames(ismember(e_frames,frames_clean));
        if length(nc14_frames) > 2% skip nuclei not in nc14
            x = schnitzcells(e).cenx;            
            y = schnitzcells(e).ceny;                
            s_cells(e_pass).xPos = x(ismember(nc14_frames,frames_clean));
            s_cells(e_pass).yPos = y(ismember(nc14_frames,frames_clean));             
            %Will be set to particle position for nuclei with matching
            %particle
            s_cells(e_pass).xPosParticle = NaN;
            s_cells(e_pass).yPosParticle = NaN;
            s_cells(e_pass).frames = nc14_frames';
            s_cells(e_pass).N = length(nc14_frames);
            s_cells(e_pass).Nucleus = e;                        
            s_cells(e_pass).ncID = eval([num2str(setID) '.' sprintf('%04d',e)]);
            %Will be set to mean particle position for nuclei with matiching
            %particles
            s_cells(e_pass).xMean = mean(x);
            s_cells(e_pass).yMean = mean(y);
            % time and set info
            s_cells(e_pass).time = time_clean(ismember(frames_clean,nc14_frames));                        
            s_cells(e_pass).setID = set_nums(i);
            e_pass = e_pass + 1;
        end
    end
    %%% Now Particles
    e_index = [s_cells.Nucleus]; % Index vector to cross-ref w/ particles        
    fn = cp_filenames{i}; % Get filename to store in struct           
    % iterate through traces
    particles = raw_data.CompiledParticles; % extract particle set
    j_init = j_pass;
    for j = 1:size(traces_clean,2)        
        raw_trace = traces_clean(:,j);     
        trace_start = find(~isnan(raw_trace),1);
        trace_stop = find(~isnan(raw_trace),1,'last');        
        %Creat versions with all intervening frames present (missing frames
        %appear as NaNs
        trace_full = raw_trace(trace_start:trace_stop)';               
        time_full = time_clean(trace_start:trace_stop);                
        frames_full = frames_clean(trace_start:trace_stop);
        % skip small fragments
        short_flag = 0;
        if length(raw_trace(~isnan(raw_trace))) < 3            
            continue
        elseif length(raw_trace(~isnan(raw_trace))) < 5            
            short_flag = 1;
        end
        j_pass = j_pass + 1;
        
        % Pull variables from particle structure        
        ap_positions_raw = particles(j).APpos;
        fov_xPos_raw = particles(j).xPos;
        fov_yPos_raw = particles(j).yPos;
        fluo3_raw = particles(j).Fluo3;
        fluo5_raw = particles(j).Fluo5;
        pt_frames = particles(j).Frame;
        % Only take values from CP struct corresponding to frames in
        % filtered frames
        NaN_vec = NaN(1,length(frames_full));        
        ap_positions = ap_positions_raw(ismember(pt_frames,frames_full));        
        xPos = fov_xPos_raw(ismember(pt_frames,frames_full));        
        yPos = fov_yPos_raw(ismember(pt_frames,frames_full));
        fluo3 = NaN_vec;
        fluo3(ismember(frames_full,pt_frames)) = fluo3_raw(ismember(pt_frames,frames_full));
        fluo5 = NaN_vec;
        fluo5(ismember(frames_full,pt_frames)) = fluo5_raw(ismember(pt_frames,frames_full));
        pt_frames = pt_frames(ismember(pt_frames,frames_full));
        % look for edge issues
        edge_frames = (((xDim-xPos) <= 1.5*snippet_size/2)|(xPos <= 1.5*snippet_size/2)|...
                          ((yDim-yPos) <= 1.5*snippet_size/2)|(yPos <= 1.5*snippet_size/2))>0;
        trace_struct(j_pass).edge_flag = max(edge_frames)&&~short_flag;        
        % Record info in trace struct
        trace_struct(j_pass).APAngle = APAngle;
        trace_struct(j_pass).cp_frames = pt_frames; % Keep track of frame correspondence
        trace_struct(j_pass).all_frames = frames_full;
        trace_struct(j_pass).nc14 = first_frame;
        trace_struct(j_pass).last_frame = last_frame;        
        trace_struct(j_pass).xPos = xPos;
        trace_struct(j_pass).yPos = yPos;
        trace_struct(j_pass).ap_vector = ap_positions;
        trace_struct(j_pass).fluo = trace_full;
        trace_struct(j_pass).fluo3 = fluo3;
        trace_struct(j_pass).fluo5 = fluo5;
        trace_struct(j_pass).time = time_full;
        trace_struct(j_pass).FluoError = particles(j).FluoError; %Estimated error in bkg subtraction     
        % For nuclei corresponding to particles, add fluo info and revise
        % mean position to align with particle pos
        nucleus = particles(j).Nucleus; 
        % Identifier variables        
        trace_struct(j_pass).Nucleus = nucleus;
        ncID = eval([num2str(setID) '.' sprintf('%04d',nucleus)]);
        trace_struct(j_pass).ncID = ncID;
        particle = particles(j).OriginalParticle;        
        pID = eval([num2str(setID) '.' sprintf('%04d',particle)]);                    
        trace_struct(j_pass).ParticleID = pID;
        trace_struct(j_pass).setID = setID;
        trace_struct(j_pass).source_path = fn;                        
    end
    % Add Particle Info to Nuclei
    for j = j_init+1:j_pass
        Nucleus = trace_struct(j).Nucleus;
        nc_ind = find(e_index==Nucleus);
        if length(nc_ind) ~= 1
            error('Error: Problem with Particle-Nucleus Crossref')
        end
        total_matched = total_matched + 1;
        s_cells(nc_ind).xPosParticle = trace_struct(j).xPos;
        s_cells(nc_ind).yPosParticle = trace_struct(j).yPos;
        s_cells(nc_ind).xMean = mean(trace_struct(j).xPos);
        s_cells(nc_ind).yMean = mean(trace_struct(j).yPos);                
        s_cells(nc_ind).ParticleID = trace_struct(j).ParticleID;                
    end
    nucleus_struct = [nucleus_struct  s_cells];        
end
for i = 1:length(nucleus_struct) % assign NaNs for ParticleID to silent nuclei
    if isempty(nucleus_struct(i).ParticleID)
        nucleus_struct(i).ParticleID = NaN;
    end
end

%%% Look for trace fragments that belong together. Stitch them up
%%% This is premised on the trace start times being sorted in ascending
%%% order!

nc_id_vec = [nucleus_struct.ncID];
set_index = [trace_struct.setID];
set_vec = unique(set_index);
% stitch together overlaps
remove_indices = [];
match_indices = [];
dupe_indices = [];
for s = 1:length(set_vec)
    set_trace_struct = trace_struct(set_index==set_vec(s));     
    base_indices = find(set_index==set_vec(s));
    % vectors used to asses fragment proximity
    start_t = [];
    stop_t = [];
    start_x = [];
    stop_x = [];
    start_y = [];
    stop_y = [];
    for i = 1:length(set_trace_struct)
        time = set_trace_struct(i).time;
        xPos = set_trace_struct(i).xPos;
        yPos = set_trace_struct(i).yPos;
        % add start and stop info
        start_t = [start_t time(1)];
        stop_t = [stop_t time(end)];
        start_x = [start_x xPos(1)];
        stop_x = [stop_x xPos(end)];
        start_y = [start_y yPos(1)];
        stop_y = [stop_y yPos(end)];
    end
    %%% enforce ascending sort order
    [start_t, si] = sort(start_t);
    stop_t = stop_t(si);
    start_x = start_x(si);
    stop_x = stop_x(si);
    start_y = start_y(si);
    stop_y = stop_y(si);
    base_indices = base_indices(si);   
    t_mat = repmat(start_t,length(start_t),1) - repmat(stop_t',1,length(stop_t));
    x_mat = repmat(start_x,length(start_x),1) - repmat(stop_x',1,length(stop_x));
    y_mat = repmat(start_y,length(start_y),1) - repmat(stop_y',1,length(stop_y));
    logic_mat = (sqrt((x_mat.^2 + y_mat.^2))<=5)&(t_mat>0)&(t_mat<120); % look for spatially and temporally proximate fragments
    logic_mat(eye(length(stop_y))==1) = 0; % remove diagonals
    overlap_ids = find(logic_mat) ;    
    tr_col = floor((overlap_ids-1)/length(stop_x))+1; % convert from linear to 2D
    tr_row = overlap_ids - (tr_col-1)*length(stop_x);
    tr_col = base_indices(tr_col);
    tr_row = base_indices(tr_row);        
    if length(unique(tr_col))~=length(tr_col) || length(unique(tr_row))~=length(tr_row)           
        warning('Duplicate Trace Fragments Detected. Removing.')
        col_mat = repmat(tr_col,length(tr_col),1)==repmat(tr_col',1,length(tr_col));
        row_mat = repmat(tr_row,length(tr_row),1)==repmat(tr_row',1,length(tr_row));        
        row_mat(eye(size(row_mat,1))==1)=0;
        col_mat(eye(size(col_mat,1))==1)=0;
        [~, row_ind] = find(max(row_mat));
        [~, col_ind] = find(max(col_mat));
        rm_vec = 1:length(tr_col);        
        dupe_indices = [dupe_indices tr_col(row_ind) tr_row(col_ind)];        
        tr_row = tr_row(~ismember(rm_vec,[row_ind col_ind]));
        tr_col = tr_col(~ismember(rm_vec,[row_ind col_ind]));                
    end    
    cat_fields = {'fluo','fluo3','fluo5','time','ap_vector','xPos','yPos'...
                    'cp_frames','all_frames'};
    for j = 1:length(tr_col)
        ID1 = min(tr_col(j),tr_row(j)); % keep earlier as base
        ID2 = max(tr_col(j),tr_row(j));
        if ismember(ID1,remove_indices)            
            ID1 = match_indices(ID1==remove_indices);            
        end
        % take ID info from earlier trace (smaller ID)
        base = trace_struct(ID1);
        extra = trace_struct(ID2);        
        for f = 1:length(cat_fields)
            base.(cat_fields{f}) = [base.(cat_fields{f}) extra.(cat_fields{f})];
        end
        base.edge_flag = max(base.edge_flag,extra.edge_flag);        
        % assign nucleus that is better fit
%         xp1 = nucleus_struct(nc_id_vec == base.ncID).xPos;
%         yp1 = nucleus_struct(nc_id_vec == base.ncID).yPos;
%         t1 = nucleus_struct(nc_id_vec == base.ncID).time;
%         xp2 = nucleus_struct(nc_id_vec == extra.ncID).xPos;
%         yp2 = nucleus_struct(nc_id_vec == extra.ncID).yPos;
%         t2 = nucleus_struct(nc_id_vec == extra.ncID).time;
%         f1 = ismember(round(t1),round(base.time));
%         if sum(f1)==length(base.time)             
%             r1 = sqrt((xp1(f1)-base.xPos).^2+(yp1(f1)-base.yPos).^2);
%         else
%             r1 = Inf;
%         end
%         f2 = ismember(round(t2),round(base.time));
%         if sum(f2) < length(base.time)   
%             error('asfa')
%         end
        % Add particle info to nucleus struct
        fn = fieldnames(trace_struct);
        for f = 1:length(fn)
            trace_struct(ID1).(fn{f}) = base.(fn{f});
        end        
        remove_indices = [remove_indices ID2];
        match_indices = [match_indices ID1];
    end        
end

% remove extra entries
index_vector = 1:length(trace_struct);
nc_particles = [nucleus_struct.ParticleID];
tr_particles = [trace_struct.ParticleID];
trace_struct = trace_struct(~ismember(index_vector,[remove_indices dupe_indices]));
rm_particles = tr_particles([remove_indices dupe_indices]);
for i = 1:length(rm_particles)
    nucleus_struct(nc_particles==rm_particles(i)).ParticleID = NaN;
    nucleus_struct(i).xMean = mean(nucleus_struct(i).xPos);
    nucleus_struct(i).yMean = mean(nucleus_struct(i).yPos);
    nucleus_struct(i).xPosParticle = [];
    nucleus_struct(i).yPosParticle = [];
end
%%
%%% ------------------------ Find Stripe Centers ---------------------- %%%
%%% Smoothing Kernel
kernel_radius = 30; % radius of gauss kernel...nucleus diameter ~= 20-25
kernel_sigma = 15; 
[x_ref_kernel, y_ref_kernel] = meshgrid(1:2*kernel_radius+1,1:2*kernel_radius+1);
x_ref_kernel = x_ref_kernel - kernel_radius - 1;
y_ref_kernel = y_ref_kernel - kernel_radius - 1;
r_mat = sqrt(x_ref_kernel.^2 + y_ref_kernel.^2);
g_kernel = exp(-(r_mat/(2*kernel_sigma))); % gauss kernel
g_kernel(r_mat>kernel_radius) = 0;
%%% fitting parameters
min_mat_time = 30*60; % Minimum time for appearance of mature stripe
xDim = 512; % FOV x size in pixels
yDim = 256;
search_kernel_width = 2.0; % Size of summation bin to use for optimization
new_trace_struct = []; % Store trace results 
new_nucleus_struct = []; % Store nucleus results
search_swath = 15; % Search space (degrees). Size of deviation from AP orthogonal permitted
search_struct = struct; % Structure to store results
for i = find(include_vec==20)%1:length(include_vec)   
    setID = include_vec(i);
    set_trace_struct = trace_struct([trace_struct.setID] == include_vec(i));
    set_nucleus_struct = nucleus_struct([nucleus_struct.setID] == include_vec(i));
    CenterAngle = -round(set_trace_struct(1).APAngle); 
    
    % Allow for inferred center line to deviate by at most "search_swath"
    % degrees from perpendicular to AP axis
    theta_vec = (CenterAngle-search_swath):(CenterAngle+search_swath);    
    
    % Get position and production info
    xp_all = [set_trace_struct.xPos]; % X positions    
    yp_all = [set_trace_struct.yPos]; % Y positions    
    fluo_all = [set_trace_struct.fluo]; % Fluorescence
    time_all = [set_trace_struct.time]; % Time
    time_all = time_all(~isnan(fluo_all));
    fluo_all = fluo_all(~isnan(fluo_all));
    ap_all = [set_trace_struct.ap_vector]; % AP postions    
    
    % Filter for times later than specified maturation time. Also remove
    % obs with nonpositive fluorescence
    filter = (time_all>min_mat_time)&(fluo_all>0);
    fluo_all = fluo_all(filter);      
    time_all = time_all(filter);
    xp_all = xp_all(filter);
    yp_all = yp_all(filter);    
    ap_all = ap_all(filter);
    
    % Array to store Fluorescence values per pixel
    set_frame = zeros(yDim,xDim);   
    for j = 1:length(fluo_all) % Could vectorize this...
        row = round(yDim - yp_all(j) + 1);
        set_frame(row,ceil(xp_all(j))) ...
            = set_frame(row,ceil(xp_all(j))) + fluo_all(j);  
    end    
    norm_ref_array = ones(size(set_frame));
    norm_ref_array = conv2(norm_ref_array,g_kernel,'same');
    gauss_array = conv2(set_frame,g_kernel,'same');
    gauss_array = gauss_array./norm_ref_array;    
    % Calculate AP-per-pixel calibration
    xMat = abs(repmat(xp_all,length(xp_all),1)' - repmat(xp_all,length(xp_all),1)).^2;
    yMat = abs(repmat(yp_all,length(yp_all),1)' - repmat(yp_all,length(yp_all),1)).^2;
    rMat = sqrt(yMat + xMat);
    apMat = abs(repmat(ap_all,length(ap_all),1)-repmat(ap_all,length(ap_all),1)');
    apMat = apMat(rMat>100); % Guard against potential anomalies with smal p separations
    rMat = rMat(rMat>100);
    APperPixel = max(apMat./rMat)*100; % Max value should occur when spot separation is orthogonal to AP axis
    % Calculate Fluorecence Profile for each search angle
    f_sums = []; % Track fluorescence
    t_vec = []; % Track Angles
    p_vec = []; % Projected Position
    index_vec = []; % Convenience vector for indexing
    for t = 1:length(theta_vec)
        theta = theta_vec(t);
        projection_mat = zeros(yDim,xDim); % Store projected positions
        for m = 1:size(projection_mat,1)
            for n = 1:size(projection_mat,2)                
                row = size(projection_mat,1) - m + 1; % Flip Y direction
                projection_mat(row,n) = round(cosd(atand(row/n)-theta)*(sqrt(row^2 + n^2)));
            end
        end        
        search_struct(t).p_mat = projection_mat;
        search_struct(t).theta = theta;
        % Get unique projection values 
        unique_p = unique(projection_mat);        
        % We will use mean fluo per pixel for each position along
        % projection axis
        projection_means = zeros(1,length(unique(projection_mat)));        
        for p = 1:length(unique_p)
            projection_means(p) = mean(gauss_array(projection_mat==unique_p(p)));            
        end
        projection_means = projection_means/sum(projection_means); % Normalize
        search_struct(t).p_means = projection_means;
        search_struct(t).p_index = unique_p;
  
        % Find total share of fluorescence captured by prescribed window size 
        % centered at each point along projection axis
        for o = 1:length(unique_p)                        
            center = unique_p(o);
            w = search_kernel_width/APperPixel;
            f_sums = [f_sums sum(projection_means((unique_p>=(center-w))&(unique_p<=(center+w))))];            
            t_vec = [t_vec theta];
            p_vec = [p_vec center];
            index_vec = [index_vec t];
        end       
    end    
    %find best radius from mean set after using screening for well
    %populated orientations
    best_f = max(f_sums);  
    if sum(f_sums==best_f) > 1
        warning('Degenerate Solutions...Taking Median Index')
    end
    candidate_t = t_vec(f_sums==best_f);
    candidate_p = p_vec(f_sums==best_f);
    min_t = min(candidate_t);
    max_t = max(candidate_t);    
    candidate_t = sort(candidate_t+90)-90;
    %If there are degenerate solutions, take median angle
    med_index = ceil(length(candidate_t)/2);
    best_angle = candidate_t(med_index);
    best_center = candidate_p(candidate_t==best_angle);
    best_index = unique(index_vec(t_vec==best_angle));
    if length(best_center) > 1 || length(best_index) > 1
        error('Degenerate Centers or Indices')
    end       
    best_projection_mat = search_struct(best_index).p_mat-best_center;
    plot_mat = abs(best_projection_mat*APperPixel);
    plot_mat(plot_mat>8) = Inf;
    PixelperAP = APperPixel.^-1; 
    unique_p = sort(p_vec(t_vec==best_angle));
    %%% Make Figure      
    stripe_fig = figure('Position',[100 100 1024 512]);%,'Visible','off');        
    map = flipud(jet(64));
    colormap(map)    
    hold on    
%     norm_array = 8 - 4*gauss_array/max(gauss_array(:));     
    im = imagesc(plot_mat); 
    h = colorbar;
    set(im,'AlphaData',.5);  
%     im2 = imagesc(Rnorm_array);
%     set(im2,'AlphaData',.8);  
%     colormap(parula(64))
    plot_times = time_all - min_mat_time;
    ss = scatter(xp_all,yDim - yp_all + 1,(fluo_all+15)/15,8-plot_times/max(plot_times)*4,'filled');    
    set(ss,'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
    axis([0 xDim 0 yDim]);
    grid on
    title(['Estimated Stripe Position with Fluorescence Emissions, Set ' num2str(include_vec(i))]);
    xlabel('X Dim (Pixels)')
    ylabel('Y Dim (Pixels)')    
    ylabel(h,'distance from stripe (AP)')    
    saveas(stripe_fig, [ap_pos_path, '/stripe_fits/set' num2str(include_vec(i)) '_stripe_pos.png'],'png');
    save([ap_pos_path, '/stripe_fits/set' num2str(include_vec(i)) '_s_id_mat.mat'],'best_projection_mat');
    
    %%% Make illustrative fitting figures
    if setID == 20
        %%
        stripe_fig_1 = figure('Position',[100 100 1024/2 512/2]);
        map = jet(128);
        colormap(map)    
        hold on            
        plot_times = time_all - min_mat_time;
        ss = scatter(xp_all,yp_all ,(fluo_all+15)/10,plot_times,'filled');    
        set(ss,'MarkerFaceAlpha',.15,'MarkerEdgeAlpha',0);
        axis([0 xDim 0 yDim]);
        box on        
        xlabel('x dim (pixels)')
        ylabel('y dim (pixels)')            
        saveas(stripe_fig_1, ['../../fig/paper_figures/Fig2/set' num2str(include_vec(i)) '_raw_activity.png'],'png');        
        saveas(stripe_fig_1, ['../../fig/paper_figures/Fig2/set' num2str(include_vec(i)) '_raw_activity.pdf'],'pdf');        

        stripe_fig_2 = figure('Position',[100 100 1024/2 512/2]);        
        hold on
        map = jet(128);
        colormap(map)
        imagesc(flipud(gauss_array))        
        axis([0 xDim 0 yDim]);
        box on        
        xlabel('x dim (pixels)')
        ylabel('y dim (pixels)')  
        for j = [1,length(search_struct)]
            p_mat = search_struct(j).p_mat;
            v = p_mat(70,250);
            ind = find(p_mat==v);
            x_ind = floor(ind/256);
            y_ind = 256-mod(ind,256);
            plot(x_ind,y_ind,'Color','black','LineWidth',1);
        end
        p_mat = search_struct(best_index).p_mat;
        v = p_mat(70,250);
        ind = find(p_mat==v);
        x_ind = floor(ind/256);
        y_ind = 256-mod(ind,256);
        plot(x_ind,y_ind,'Color',map(30,:),'LineWidth',2)
        saveas(stripe_fig_2, ['../../fig/paper_figures/Fig2/set' num2str(include_vec(i)) '_conv_activity.png'],'png');        
        saveas(stripe_fig_2, ['../../fig/paper_figures/Fig2/set' num2str(include_vec(i)) '_conv_activity.pdf'],'pdf');        

        stripe_fig_3 = figure('Position',[100 100 1024/2 512/2]);  
        hold on
        for j = 1:length(search_struct)
            if j == best_index
                plot(search_struct(j).p_means,'Color',map(30,:),'LineWidth',2)
            else
                plot(search_struct(j).p_means,'Color',[1 1 1]/2,'LineWidth',.75)
            end
        end
        box on
        xlabel('position along projection axis')
        ylabel('average fluorescence per pixel')
        saveas(stripe_fig_3, ['../../fig/paper_figures/Fig2/set' num2str(include_vec(i)) '_best_profile.png'],'png');        
        saveas(stripe_fig_3, ['../../fig/paper_figures/Fig2/set' num2str(include_vec(i)) '_best_profile.pdf'],'pdf');        

        stripe_fig_4 = figure('Position',[100 100 568 256]);
        hold on
        map = flipud(jet(128));
        colormap(map)                   
        im = imagesc(flipud(plot_mat)); 
        h = colorbar;
        ss = scatter(xp_all,yp_all + 1,(fluo_all+15)/10,[1 1 1]/2,'filled');    
        set(ss,'MarkerFaceAlpha',.15,'MarkerEdgeAlpha',.15,'MarkerEdgeColor','black');
        xlabel('x dim (pixels)')
        ylabel('y dim (pixels)')    
        ylabel(h,'distance from stripe center (%AP)')    
        axis([0 512 0 256])
        box on
        saveas(stripe_fig_4, ['../../fig/paper_figures/Fig2/set' num2str(include_vec(i)) '_final_stripe.png'],'png');        
        saveas(stripe_fig_4, ['../../fig/paper_figures/Fig2/set' num2str(include_vec(i)) '_final_stripe.pdf'],'pdf');            
        %%
    end
    

    %Assign Relative AP Position to Each Particle
    for j = 1:length(set_trace_struct)
        x_vec = set_trace_struct(j).xPos;
        y_vec = yDim - set_trace_struct(j).yPos + 1;
        ap_pos_vec = [];
        for k = 1:length(x_vec)
            ap_pos_vec = [ap_pos_vec APperPixel*best_projection_mat(y_vec(k),x_vec(k))];
        end
        set_trace_struct(j).rel_ap_vector = ap_pos_vec;
        set_trace_struct(j).MeanAP = mean(ap_pos_vec);
    end    
    for j = 1:length(set_nucleus_struct)
        x_vec = round(set_nucleus_struct(j).xPos);
        y_vec = yDim - round(set_nucleus_struct(j).yPos) + 1;
        ap_pos_vec = [];
        for k = 1:length(x_vec)
            ap_pos_vec = [ap_pos_vec APperPixel*best_projection_mat(y_vec(k),x_vec(k))];
        end
        set_nucleus_struct(j).rel_ap_vector = ap_pos_vec;
        set_nucleus_struct(j).MeanAP = mean(ap_pos_vec);        
    end
    % Update set struct with search info
    for j = 1:length(set_trace_struct)                
        set_trace_struct(j).PixelperAP = PixelperAP;        
        set_trace_struct(j).stripe_angle = best_angle + 90;
        set_trace_struct(j).deviation = CenterAngle-best_angle;
        set_trace_struct(j).projection_center = best_center;
        set_trace_struct(j).search_radius = search_kernel_width;        
%         set_trace_struct(j).MeanAP = mean(set_trace_struct(j).ap_vector);               
        set_trace_struct(j).search_swath = search_swath;
    end
    new_trace_struct = [new_trace_struct set_trace_struct]; 
    new_nucleus_struct = [new_nucleus_struct set_nucleus_struct]; 
    disp(['Completed ' num2str(i) ' of ' num2str(length(include_vec))])
end         

% save([trace_name '.mat'],'new_trace_struct') 
% save([nucleus_name '.mat'],'new_nucleus_struct') 
%%
%%%% ------------------------- Make QC Plots ----------------------------%%%

%Make Titles for Plots (This will only work for eve2 format set titles
%currently)
sets = unique({trace_struct.source_path});
n_sets = length(sets);
set_titles = {};
for i = 1:length(sets)
    start = strfind(sets{i},'sec_') + 4;
    stop = strfind(sets{i},'/Compiled') - 1;
    string = sets{i};
    set_titles = {set_titles{:}  string(start:stop)};        
end

%Set dimensions for figs
xDim = 1;
yDim = 1;
toggle = 1;
while xDim*yDim < length(include_vec) 
    if toggle == 1
        xDim = xDim + 1;
    else 
        yDim = yDim + 1;
    end
    toggle = toggle * (toggle == 0) + toggle == 0;
end

precision = 0.01; % Set granularity (.01 = 1% AP precision)

cm = jet(64); % Make Color Palettes for use in figures
%Array to store color mappings
increment = floor(60 / n_sets);
set_colors = zeros(n_sets, 3);
for i = 1:n_sets
    set_colors(i,:) = cm(1+(i-1)*increment,:);
end

%------------------------- AP Averages -----------------------------------%
close all
ap_vec = round([new_trace_struct.MeanAP]);
set_vec = [new_trace_struct.setID];
ap_index = unique(ap_vec);
mean_fluo_ap_set = NaN(length(ap_index),length(include_vec));
cum_fluo_ap_set = NaN(length(ap_index),length(include_vec));
mean_fluo_ap = NaN(1,length(ap_index));
cum_fluo_ap = NaN(1,length(ap_index));
for a = 1:length(ap_index)
    ap = ap_index(a);
    for s = 1:length(include_vec)
        setID = include_vec(s);
        nan_fluo_vec_set_ap = [new_trace_struct((ap_vec==ap)&...
            (set_vec==setID)).fluo];
        z_fluo_vec_set_ap = nan_fluo_vec_set_ap; 
        z_fluo_vec_set_ap(isnan(z_fluo_vec_set_ap)) = 0;
        mean_fluo_ap_set(a,s) = nanmean(nan_fluo_vec_set_ap);
        cum_fluo_ap_set(a,s) = mean(z_fluo_vec_set_ap);
    end   
    nan_fluo_vec_ap = [new_trace_struct(ap_vec==ap).fluo];
    z_fluo_vec_ap = nan_fluo_vec_ap;
    z_fluo_vec_ap(isnan(z_fluo_vec_ap)) = 0;
    mean_fluo_ap(a) = nanmean(nan_fluo_vec_ap);
    cum_fluo_ap(a) = mean(z_fluo_vec_ap);
end
mean_fluo_fig = figure('Position',[0 0 1024 1024]);
for j = 1:n_sets        
    subplot(yDim,xDim,j)
    hold on
    bar(ap_index, mean_fluo_ap  , 'FaceColor','black',...
        'FaceAlpha', .3,'EdgeColor','black','EdgeAlpha',.3,'BarWidth', 1);
    bar(ap_index, mean_fluo_ap_set(:,j) , 'FaceColor',set_colors(j,:),...
        'FaceAlpha', .5,'EdgeColor',set_colors(j,:),'BarWidth', 1);
    set(gca,'fontsize',4)
    title(['Mean Fluorescence per Time Step (active tp) NC 14, Set: ' set_titles{j}]); %' Set:' sets{j}])
    axis([min(ap_index),max(ap_index),0 , 1.2*max(mean_fluo_ap_set(:))])    
    grid on
end
saveas(mean_fluo_fig, [ap_pos_path, 'mean_fluo_ap.png'],'png');

%%% ----------- Integrated Fluorescence With Stripe Centers ----------- %%%
cumulative_fluo_fig = figure('Position',[0 0 1024 1024]);
for j = 1:n_sets        
    subplot(yDim,xDim,j)
    hold on
    bar(ap_index, cum_fluo_ap  , 'FaceColor','black',...
        'FaceAlpha', .3,'EdgeColor','black','EdgeAlpha',.3,'BarWidth', 1);
    bar(ap_index, cum_fluo_ap_set(:,j) , 'FaceColor',set_colors(j,:),...
        'FaceAlpha', .5,'EdgeColor',set_colors(j,:),'BarWidth', 1);
    set(gca,'fontsize',4)
    title(['Mean Fluorescence per Time Step (all tp) NC 14, Set: ' set_titles{j}]); %' Set:' sets{j}])
    axis([min(ap_index),max(ap_index),0 , 1.2*max(cum_fluo_ap_set(:))])    
    grid on
end
saveas(cumulative_fluo_fig, [ap_pos_path, 'cumulative_fluo_ap.png'],'png');

%%%%%%%%%% ---------------- Temporal Trends -------------------- %%%%%%%%%%
close all
ap_vec_2 = round([new_trace_struct.MeanAP]/2)*2;
ap_2_index = unique(ap_vec_2);
ap_2_index = ap_2_index(abs(ap_2_index) < 8);
n_boots = 100; % number of boostraps to use for avg trends

line_types = {'-','-o','--'};
set_line_types = {};
for i = 1:length(include_vec)
    set_line_types = [set_line_types{:} {line_types{1+mod(i,3)}}];
end    

% set_fluo_mean_array = NaN(length(include_vec),45);
% set_fluo_se_array = NaN(length(include_vec),45);
for i = 1:length(ap_2_index)
    p = [];
    ap = ap_2_index(i);    
    sets = unique(set_vec(ap_vec_2==ap));    
    temporal_fig = figure('Position',[0 0 1024 512]);
    hold on
    time_vec = 1:45;
    f_avg = zeros(1,length(time_vec));
    legend_string = {};
    for j = 1:length(sets)
        setID = sets(j);
        set_trace_struct = new_trace_struct((ap_vec_2==ap)&(set_vec==setID));
        if length([set_trace_struct.fluo]) < 100 % skip poorly represented sets
            continue
        end
        time_list = ceil([set_trace_struct.time]/60);
        fluo_set = [set_trace_struct.fluo];
        n_dp = length(fluo_set);
        sample_vec = 1:n_dp;
        f_array = NaN(n_boots,length(time_vec));
        for n = 1:n_boots
            s_ids = randsample(sample_vec,n_dp,true);
            fluo_boot = fluo_set(s_ids);
            time_boot = time_list(s_ids);
            for t = time_vec
                f_array(n,t) = nanmean(fluo_boot(ismember(time_boot,t-2:t+2)));
            end
        end
        sfm = nanmean(f_array);
        sfe = nanstd(f_array);
%         set_fluo_mean_array(j,:) = sfm;
%         set_fluo_se_array(j,:) = sfe;
        plot([time_vec ; time_vec],[sfm-sfe ; sfm + sfe],'Color',set_colors(include_vec==sets(j),:),'LineWidth',1.5)
        p = [p plot(time_vec, sfm, set_line_types{include_vec==sets(j)},'Color',set_colors(include_vec==sets(j),:),'LineWidth',1.5)];                
        legend_string = {legend_string{:} ['Set ' num2str(sets(j))]};
    end
    title(['Average Fluorescence Over Time (5 min moving avg), AP ' num2str(ap)]);
    legend(p,legend_string{:},'Location','northwest')    
    saveas(temporal_fig,[fluo_path 'time_trends_bin' num2str(ap) '.png'],'png')
end
