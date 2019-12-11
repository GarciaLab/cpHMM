% Experiment with sochastic simulations in which rates vary in time
function temp_gillespie = synthetic_rate_gillespie_temporal(seq_length,...
            alpha, w, k_on_lower,k_on_upper,k_off_lower,k_off_upper, deltaT, r_emission, noise, ...
            pi0,t_shift,shift_width,granularity,f_per_mRNA,load_lag)
%%% Simulation Parameters
% t_process = 3000; % seconds (say)
t_process = deltaT*seq_length;
t_MS2 = alpha*deltaT;
% arrival_rate = r_emission(2)./f_per_mRNA;
t_ref = 0:granularity:t_process; % assuming that jump lengths are much shorter than obs time
beta = -log(.16)/shift_width/deltaT;
k_on_ref = k_on_lower + (k_on_upper-k_on_lower)*(1 + exp((t_ref+.5*granularity-t_shift*deltaT)*beta)).^-1;
k_off_ref = k_off_upper - (k_off_upper-k_off_lower)*(1 + exp((t_ref+.5*granularity-t_shift*deltaT)*beta)).^-1;

%%% Make Rate Array (3D rate matrix)
R_array = zeros(3,3,length(t_ref));
% from 0
R_array(1,1,:) = -2*k_on_ref;
R_array(2,1,:) = 2*k_on_ref;
R_array(3,1,:) = zeros(1,length(t_ref));
% from 1
R_array(1,2,:) = k_off_ref;
R_array(2,2,:) = -(k_on_ref+k_off_ref);
R_array(3,2,:) = k_on_ref;
% from 3
R_array(1,3,:) = zeros(1,length(t_ref));
R_array(2,3,:) = 2*k_off_ref;
R_array(3,3,:) = -2*k_off_ref;

%%% Simulations
jump_times = [0];
promoter_states = [randsample(1:3,1,true,pi0)];
obs_time_grid = deltaT:deltaT:t_process;
state_vec = 1:3;
for ts = 1:length(t_ref)
    T = t_ref(ts);
    t_step = 0; % process time within discrete step
    R = R_array(:,:,ts); % extract rate matrix for time step    
    cs = promoter_states(end); % get current state
    options = state_vec(state_vec~=cs);
    while t_step < granularity
        tau = -1/R(cs,cs);
        dt = exprnd(tau); % select jump time
        t_step = t_step + dt;
        options = state_vec(state_vec~=cs);        
        % if jump time is within discrete window, record
        if t_step < granularity 
            jump_times = [jump_times T + t_step];
            cs = randsample(options,1,true,R(options(:),cs)');
            promoter_states = [promoter_states cs];            
        end
    end
end
%%% Generate trace from promoter trajectory
e_states = r_emission(promoter_states);
% jump_times = [0 jump_times];
fluo_grid = zeros(1,length(obs_time_grid));
fluo_grid_MS2 = zeros(1,length(obs_time_grid));
for window = 1:length(fluo_grid)
    t_end = obs_time_grid(window);  
    t_start = max(0,t_end - w*deltaT);
    s_ind = find(jump_times<t_start,1,'last'); % index of first relevant promoter state
    if isempty(s_ind)
        s_ind = 1;
    end
    e_ind = find(jump_times<t_end,1,'last'); % last relevant state
    times = jump_times(s_ind:e_ind);
    times(1) = t_start;
    times = [times t_end];
    states = e_states(s_ind:e_ind);
    F = 0;
    F_alpha = 0;
    for t = 1:(length(times)-1)
        t1 = t_end - times(t);
        t2 = t_end - times(t+1);
%         if t1-t2 < 10 % assume that minimum of 1 sec window needed for loading
%             continue
%         end
        ta1 = min(t1,t_MS2);
        ta2 = min(t2,t_MS2);
        F = F + states(t)*(t1 - t2);
        F_alpha = F_alpha + states(t)*(ta1-ta2)*(ta1+ta2)/t_MS2*.5 + ...
            ((t1 - t2)-(ta1-ta2))*states(t);
    end    
    fluo_grid(window) = F;    
    fluo_grid_MS2(window) = F_alpha;
end
noise_vec = normrnd(0,noise,1,length(fluo_grid));
fluo_grid_noise = fluo_grid + noise_vec;
fluo_grid_noise(fluo_grid_noise<0) = 0;
fluo_grid_MS2_noise = fluo_grid_MS2 + noise_vec;
fluo_grid_MS2_noise(fluo_grid_MS2_noise<0) = 0;
% make clipped version with everything after last real nonzer obs == 0
last_time = find(round(fluo_grid),1,'last');
fluo_grid_clipped = fluo_grid_noise;
fluo_grid_clipped = fluo_grid_clipped(1:last_time);
fluo_grid_MS2_clipped = fluo_grid_MS2_noise;
fluo_grid_MS2_clipped = fluo_grid_MS2_clipped(1:last_time);
%%% Simulate Poisson Initiation Events
jump_times = [jump_times t_process];
arrival_times = [];
% calculate adjustment parameter for lag time
R_init = R_array(:,:,1);
[v,d] = eig(R_init);
min_index = find(diag(d)==0);
steady = v(:,min_index)/sum(v(:,min_index));
avg_dwell = -steady(2)*(R_init(2,2)^-1) + -steady(3)*(R_init(3,3)^-1);
shift_factor = avg_dwell/(avg_dwell-load_lag);

for s = 1:length(promoter_states)
    state = promoter_states(s);
    ar = r_emission(state)/f_per_mRNA*shift_factor;
    if ar > 0
        start_time = jump_times(s);
        stop_time = jump_times(s+1);
        T = stop_time - start_time;
        if s > 1
            if ismember(promoter_states(s-1),[1,2])
                t = load_lag;   
            else
                t = 0;
            end
        else
            t = 0;
        end                        
        while t < T
            dt = exprnd(1/ar);
            t = t + dt;
            if t < T
                arrival_times = [arrival_times start_time + t];
            end
        end
    end
end
% Generate poisson trace
fluo_grid_MS2_poisson = zeros(1,length(obs_time_grid));
for window = 1:length(fluo_grid_MS2_poisson)
    t_end = obs_time_grid(window);  
    t_start = max(0,t_end - w*deltaT);
    polII_arrivals = abs(arrival_times(arrival_times>t_start&arrival_times<=t_end) - t_end);
    for p = 1:length(polII_arrivals)
        at = polII_arrivals(p);
        if at < t_MS2
            fluo_grid_MS2_poisson(window) = fluo_grid_MS2_poisson(window) + at/t_MS2 * f_per_mRNA;
        else
            fluo_grid_MS2_poisson(window) = fluo_grid_MS2_poisson(window) + f_per_mRNA;
        end
    end
end

temp_gillespie.fluo_poisson_MS2_no_noise = fluo_grid_MS2_poisson;
temp_gillespie.fluo_no_noise = fluo_grid;
temp_gillespie.fluo = fluo_grid_noise;
temp_gillespie.fluo_clipped = fluo_grid_clipped;
temp_gillespie.fluo_MS2_no_noise = fluo_grid_MS2;
temp_gillespie.fluo_MS2 = fluo_grid_MS2_noise;
temp_gillespie.fluo_MS2_clipped = fluo_grid_MS2_clipped;
temp_gillespie.naive_states = promoter_states;
temp_gillespie.jump_times = jump_times;
temp_gillespie.t_ref = t_ref;
temp_gillespie.k_on_ref = k_on_ref;
temp_gillespie.k_off_ref = k_off_ref;
temp_gillespie.arrival_times = arrival_times;