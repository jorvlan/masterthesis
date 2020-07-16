function tremor_1st_level_TOS(subject)

subject = '001';
clc; close all; clearvars -except subject;

%% create correct path, add functions and scripts
restoredefaultpath;
addpath(genpath('/project/3024005.01/Analysis/SupportingScripts/'));
addpath('/project/3024005.01/Analysis/Tremor');
addpath /home/common/matlab/fieldtrip
ft_defaults %including all default Fieldtrip functions

%% Define input/output directories and files
conf.input.rawdir = '/project/3024005.01/DataVeni/';
conf.input.eegdir = fullfile(pf_findfile(conf.input.rawdir,['/' subject '/'],'fullfile'),[subject '_brainvision/raw']);
conf.input.eegfile = pf_findfile(conf.input.eegdir,['/' subject '/&/TOS/&/.eeg/'],'fullfile');

% conf.output.maindir = [ '/project/3024005.01/Analysis/Tremor/NewResults/' subject ];
% conf.output.figdir = [ '/project/3024005.01/Analysis/Tremor/NewResults/' subject '/Figures' ];
% conf.output.figsubdir = [ '/project/3024005.01/Analysis/Tremor/NewResults/' subject '/Figures/SingleTrials' ];
% conf.output.datadir = [ '/project/3024005.01/Analysis/Tremor/NewResults/' subject '/Data' ];
conf.output.maindir = [ '/project/3024005.01/Analysis/Tremor/NewResults3/' subject '/_0-60' ];
conf.output.figdir = [ '/project/3024005.01/Analysis/Tremor/NewResults3/' subject '/_0-60/Figures' ];
conf.output.figsubdir = [ '/project/3024005.01/Analysis/Tremor/NewResults3/' subject '/_0-60/Figures/SingleTrials' ];
conf.output.datadir = [ '/project/3024005.01/Analysis/Tremor/NewResults3/' subject '/_0-60/Data' ];

% Create all directories unless they exist already
temp = fieldnames(conf.output);
for freq=1:numel(temp)
    if ~exist( getfield(conf.output,temp{freq}) )
        mkdir( getfield(conf.output,temp{freq}) );
    end
end
conf.output.save = 'yes'; 

%% Create trial configuration and read raw data into trials
cfg = [];
cfg.dataset             = conf.input.eegfile;
[dir,name,ext]          = fileparts(cfg.dataset);
cfg.channel             = [2:8];                    %To discard channels that are not used here. Channel 2 to 5 is EMG, 6 to 8 is ACC.
cfg.trialfun            = 'ft_trialfun_general';
cfg.trialdef.eventtype  = 'Stimulus';
cfg.trialdef.eventvalue = {'S 70','S 71','S 80','S 81'};    
cfg.trialdef.prestim    = 6; % 7 for -6 %also read in 1s before and after trial for mtmconvol later (similar to padding around trial)
cfg.trialdef.poststim   = 64; %31 for 30 sec trials, 64 for + 3
trialcfg = ft_definetrial(cfg); 
%Read in all the raw data
cfg = trialcfg;
All_raw = ft_preprocessing(cfg);

%% Preprocessing EMG and ACC data
%EMG - step1 of preprocessing
cfg = trialcfg;
cfg.channel             = [2:5];    % Select the EMG channels
cfg.continuous          = 'yes';    % Load all data, select the trials later
cfg.preproc.detrend     = 'yes';    % Detrend to remove linear trend from the data (done per trial) (default = 'no')
cfg.preproc.demean      = 'yes';    % Demean  whether to apply baseline correction (default = 'no')
cfg.preproc.bsfilter    = 'yes';
cfg.preproc.bsfreq      = [2 8];    % frequency range
cfg.preproc.bsfiltord   = 1;        % Order of BS filter
cfg.preproc.bsfilttype  = 'but';    % Type of BS filter
cfg.preproc.bpfilter	= 'yes';    % BP filter
cfg.preproc.bpfreq	    = [20 250]; % BP frequency
cfg.preproc.bpfiltord   = 1;        % Order of BP filter
cfg.preproc.bpfilttype  = 'but';    % Type of BP filter
cfg.preproc.rectify     = 'yes';    % Rectify for tremor burst and to regain low frequencies
EMGdata_preproc1 = ft_preprocessing(cfg);

%Additional artifact identification to see if removal of data does improve data quality
cfg = [];
cfg.continuous = 'yes' ;
cfg.artfctdef.zvalue.channel = 'all';
cfg.artfctdef.zvalue.cutoff = 4;
% cfg.artfctdef.zvalue.interactive = 'yes'; %Only use if you visually want to check out the artifact rejection
[cfg, artifact_jump] = ft_artifact_zvalue(cfg,EMGdata_preproc1); 

%Adjusting the artifact length
midpoint = mean(artifact_jump,2);
artifact_jump = [midpoint-150 midpoint+150]; %originally -150 +150
cfg= [];
cfg.artfctdef.reject = 'nan';
cfg.artfctdef.jump.artifact = artifact_jump;
EMGdata_no_artifacts = ft_rejectartifact(cfg,EMGdata_preproc1);

%Replace nan with 0 -- an alternative could be to interpolate data
for freq=1:size(EMGdata_no_artifacts.trial,2)
    EMGdata_no_artifacts.trial{freq}(isnan(EMGdata_no_artifacts.trial{freq}))=0;
end

% EMG - step 2 of preprocessing
cfg = [];
cfg.preproc.bpfilter   = 'yes';
cfg.preproc.bpfreq     = [2 20];
cfg.preproc.bpfiltord  = 1;
EMGdata_preproc2 = ft_preprocessing(cfg,EMGdata_preproc1);

% ACC - only 1 step of preprocessing
cfg = trialcfg;
cfg.channel             = [6:8];
cfg.continuous          = 'yes';
cfg.preproc.detrend     = 'yes';
cfg.preproc.demean      = 'yes';
cfg.preproc.bpfreq      = [2 20];
cfg.preproc.bpfiltord   = 1;
ACC_preproc = ft_preprocessing(cfg);

% Combine preprocessed EMG and ACC data
cfg=[];
All_processed = ft_appenddata(cfg,EMGdata_preproc2,ACC_preproc);

%Remove some temporary variables
clearvars dir ext name temp midpoint artifact_jump;

%% Redefine trials into segment of 5 seconds & Calculate average power spectrum for EMG and ACC over all trial data (all frequencies and channels)
cfg = [];
cfg.toilim = [0 60]; 
segdata = ft_redefinetrial(cfg,All_processed); %for analysis using fft drop second before and after the trial
cfg = [];
cfg.length  = 5; 
cfg.overlap = 0;
cfg.trials  = 'all';
segdata = ft_redefinetrial(cfg,segdata); % cut data into 5s segments

cfg = [];
cfg.method      = 'mtmfft';     %mtmfft analyses an entire spectrum for the entire data length, implements multitaper frequency transformation
cfg.foi         = 2:0.2:16;     %frequency range you are interested in. Freq resolution depends on the length of the time window (1/T).
cfg.taper       = 'hanning';
cfg.t_ftimwin   = ones(length(cfg.foi),1)*5;
cfg.toi         = 'all';        %all trials (trial now =5s segment)
cfg.keeptrials  = 'no';         %return average (default = no)
fft_freq = ft_freqanalysis(cfg,segdata);

%% Plot the powerspectrum of all ACC channels & the EMG channels (channel 5 - 7) of the most affected arm (channel 1 & 2)
fig = figure('units','centimeters','outerposition',[0 0 20 20],'Color',[1 1 1]);
hold on;
plot(fft_freq.freq, fft_freq.powspctrm(5,:)) %acc_x
plot(fft_freq.freq, fft_freq.powspctrm(6,:)) %acc_y
plot(fft_freq.freq, fft_freq.powspctrm(7,:)) %acc_z
legend('ACC_x','ACC_y','ACC_z')
xlabel('Frequency (Hz)');
ylabel('Absolute tremor power (µV^2)');
title([ subject ' TOS - ACC powerspectrum']);
if strcmp(conf.output.save,'yes')
    fig_name = [ subject '_TOS_ACC_pspec' ];
    saveas(gcf, fullfile(conf.output.figdir,fig_name), 'jpg');
    saveas(gcf, fullfile(conf.output.figdir,fig_name), 'fig');
end

fig = figure('units','centimeters','outerposition',[0 0 20 20],'Color',[1 1 1]);
hold on;
plot(fft_freq.freq, fft_freq.powspctrm(1,:)) %ECR_ma
plot(fft_freq.freq, fft_freq.powspctrm(2,:)) %FCR_ma
legend('ECR_{ma}','FCR_{ma}')
xlabel('Frequency (Hz)');
ylabel('Absolute tremor power (µV^2)');
title([ subject ' TOS - EMG powerspectrum']);
if strcmp(conf.output.save,'yes')
    fig_name = [ subject '_TOS_EMG_pspec' ];
    saveas(gcf, fullfile(conf.output.figdir,fig_name), 'jpg');
    saveas(gcf, fullfile(conf.output.figdir,fig_name), 'fig');
end

%% Read out the tremor frequency and the peak channel - adjust manually if necessary  
[indexACC,indexFreq] = find(fft_freq.powspctrm(:,:) == max(max(fft_freq.powspctrm((5:7),:))));

%Manually adjust peak frequency and channel for subjects that have most power in the first harmonic (tremor frequency * 2)
if strcmp(subject, '007')
    indexACC = 5; 
    indexFreq = 14; %14 in 0.2 Hz resolution & 6 in 0.5 Hz resolution       %minor mistake in previous runs: it was 15 (0.2 Hz res) instead of 14 
elseif strcmp(subject, '014')
    indexACC = 5;
    indexFreq = 13; %13 in 0.2 Hz resolution & 6 in 0.5 Hz resolution
end

% optional: calculate principal component
% [~,pscore,~,~,pexp] = pca([acc_x acc_y acc_z]);           % calculate principle component analysis, rows are observations, columns variables

%% Calculate the powerspectrum for all 5s segments (=NumTrials * 12) and for all frequencies and channels 
power_array = [];
for segment = 1:length(segdata.trial)
    cfg=[];
    cfg.method      = 'mtmfft';
    cfg.output      = 'pow'; %=default
    cfg.trials      = segment;
    cfg.foi         = 2:0.2:16;
    cfg.taper       = 'hanning';
    cfg.t_ftimwin   = ones(length(cfg.foi),1)*5;
    cfg.toi         = 'all';
    cfg.keeptrials  = 'yes';
    fft_segments = ft_freqanalysis(cfg,segdata);
    
    format short g; %sets output format (short g  --> float or fixed (whatever is best) and 5 digits)
    power_array = cat(1,power_array, [fft_segments.powspctrm]); %add to array, which gets bigger every repetition
end

%% Calculate mean power in peak channel per 60s trials, plot and save data
%Calculate the mean power at the peak frequency (Hz) (bin) across the 12 segments of each trial from the fft_freq.powspctrm array
mean_power_per_trial_peakfreq = [];
for freq=1:17 %iterate over trials to get average power per trial
    
    total_power_per_seg_peakfreq = [];
    for time_range=1:12 %iterate over segements to get total power over peak frequency
        sum_freqbins_peakfreq = sum(power_array( ((freq-1)*12 + time_range), indexACC, indexFreq)); % select only the power at the peak frequency
        total_power_per_seg_peakfreq = [ total_power_per_seg_peakfreq sum_freqbins_peakfreq ];
    end
     %mean_power_peakfreq = mean(total_power_per_seg_peakfreq);
    mean_power_peakfreq = mean(log(total_power_per_seg_peakfreq));
    mean_power_per_trial_peakfreq = [ mean_power_per_trial_peakfreq mean_power_peakfreq ];

end
%log_mean_power_per_trial_peakfreq = log(mean_power_per_trial_peakfreq);
log_mean_power_per_trial_peakfreq = mean_power_per_trial_peakfreq;

%Plot power values unsorted (visualize trend over time) then sort the power values and save them
[sorted_trials, indeces_trials] = sort(All_processed.trialinfo);
sorted_trials = [sorted_trials(1:5); sorted_trials(7:17); sorted_trials(6)]; %Move 2nd Threat_odd block (=shock block) to the end
indeces_trials = [indeces_trials(1:5); indeces_trials(7:17); indeces_trials(6)]; %Use this array also later for sorting time courses trials

fig = figure('units','centimeters','outerposition',[0 0 25 20],'Color',[1 1 1]);
plot(1:17, log_mean_power_per_trial_peakfreq) %Plot the average log tremor power for the 17 trials (in sequential order) 
xlabel('Trial number');
ylabel('Average tremor power (log µV^2)');
title([ subject ' TOS - Average tremor power over trials']);
h = xline(indeces_trials(end),'--r','Shock trial');
%Save the plot
if strcmp(conf.output.save,'yes')
    fig_name = [ subject '_TOS_average_log_tremor_power_over_trials' ];
    saveas(gcf, fullfile(conf.output.figdir,fig_name), 'jpg');
    saveas(gcf, fullfile(conf.output.figdir,fig_name), 'fig');
end

%Resort trials (power values) and save the data
log_mean_power_per_trial_peakfreq = log_mean_power_per_trial_peakfreq(1,indeces_trials); %New order trials: thr_fix, thr_odd, safe_fix, safe_odd, shock

if strcmp(conf.output.save,'yes')
    power_values = log_mean_power_per_trial_peakfreq;
    save(fullfile(conf.output.datadir,[ subject '_TOS_average_log_power_per_trial']),'power_values');
end


%%%
% %threat-fix tremor power values average of the 4 trials
% Thr_fix_mean_tremor_power_values_n17(17,1) = mean(power_values(1,1:4));
% 
% %threat-odd tremor power values average of the 4 trials
% Thr_odd_mean_tremor_power_values_n17(17,1) = mean(power_values(1,5:8));
% 
% %safe-fix tremor power values average of the 4 trials
% Safe_fix_mean_tremor_power_values_n17(17,1) = mean(power_values(1,9:12));
% 
% %safe-odd tremor power values average of the 4 trials
% Safe_odd_mean_tremor_power_values_n17(17,1) = mean(power_values(1,13:16));
% 
% %shock
% shock_mean_tremor_power_values_n17(17,1) = power_values(:,17);


%% Calculate TFR for all trials and average trial per condition (i.e., create the 60s log power time course)
cfg = [];
cfg.method      = 'mtmconvol';                           % Select method (choose 'mtmconvol')
cfg.output      = 'pow';                                 % power  
cfg.taper       = 'hanning';                             % Windowing (because cut-off frequency), (Choose 'hanning' for low frequency)
cfg.foi         = 2:0.5:16; %2:0.2:16;                             % frequency range you are interested in (usually 1:0.5:20, make sure you at least include 3-8 Hz)   
nFoi            = length(cfg.foi);
cfg.t_ftimwin   = repmat(2,1,nFoi);                      % Wavelet length (seconds; 1 wavelet per frequency). This is important also for your NaN in the hanning taper (which is 0.5*this) a 2 seconds sliding window
cfg.toi         = All_processed.time{1}(1):0.001:All_processed.time{1}(end);% timeline the TFR (resolution in seconds) ('orig': original resolution; 'timedat': one step specified under conf.prepemg.timedat;)
cfg.pad         = 'maxperlen'; %'maxperlen'                          % Padding (use 'maxperlen')
cfg.keeptrials  = 'yes'; %no creates a 3D variable, yes a 4D variable
%cfg.trials       = [1:2];
mtmconvol_freq_trials = ft_freqanalysis(cfg, All_processed);  % The input data still contains 1s of data before and after the trial, to avoid NaNs at the beginning/end

%Match max freq from fft (0.2 Hz resolution) with frequency from mtmconvol (0.5 Hz resolution)
%     up = ceil(fft_freq.freq(indexFreq) * 2) / 2;
%     down = floor(fft_freq.freq(indexFreq) * 2) / 2;
%     if abs(fft_freq.freq(indexFreq) - up) < abs(fft_freq.freq(indexFreq) - down)
%         new_indexFreq = up;
%     else
%         new_indexFreq = down;
%     end


new_Freq=round(fft_freq.freq(indexFreq) * 2) / 2;
new_indexFreq = find(mtmconvol_freq_trials.freq == new_Freq);
%new_indexFreq = 8; %needed to be adjusted for sub 5, 11, 14, 15 , 16 when 0.5 Hz resolution.

%% Plot the powerspectrum MTMCONVOL of all ACC channels & the EMG channels (channel 5 - 7) of the most affected arm (channel 1 & 2)
all_mean_power_per_freq = [];
for trial=1:size(mtmconvol_freq_trials.powspctrm,1)
    
    trial_mean_power_per_freq = [];
    for freq = 1:size(mtmconvol_freq_trials.powspctrm,3)
        time_range = (6001:66000);
        mean_power_per_freq = mean(mtmconvol_freq_trials.powspctrm(trial,indexACC,freq,time_range));
        trial_mean_power_per_freq = [ trial_mean_power_per_freq mean_power_per_freq ];
    end
    
    all_mean_power_per_freq = [ all_mean_power_per_freq; trial_mean_power_per_freq ];
    
end

power_vals = mean(all_mean_power_per_freq,1);
plot(mtmconvol_freq_trials.freq, power_vals)
xlabel('Frequency (Hz)');
ylabel('Absolute tremor power (µV^2)');
title([ subject ' TOS - ACC powerspectrum']);
if strcmp(conf.output.save,'yes')
    fig_name = [ subject '_new_TOS_ACC_pspec' ];
    saveas(gcf, fullfile(conf.output.figdir,fig_name), 'jpg');
    saveas(gcf, fullfile(conf.output.figdir,fig_name), 'fig');
end




%Calculate and plot the (log transformed) tremor power time courses (TFR) of all 17 trials (trials of the same condition plotted in 1 figure)
power_over_time = {};
conditions = { 'threat fix', 'threat odd', 'safe fix', 'safe odd', 'shock' };
for condition=1:length(conditions)
 
   
    fig = figure('units','centimeters','outerposition',[0 0 25 20],'Color',[1 1 1]);
    hold on
    
    %Plot all trials of one condition in one plot
    if condition < length(conditions)
        for trial=1:4
            total_trial_num = indeces_trials( (condition-1)*4 + trial );
            power_over_time{condition}(trial,:) = log(squeeze(mtmconvol_freq_trials.powspctrm(total_trial_num,indexACC,new_indexFreq,1001:(length(mtmconvol_freq_trials.time)-1000))))';
            plot(mtmconvol_freq_trials.time(1001:(length(mtmconvol_freq_trials.time)-1000)), power_over_time{condition}(trial,:)) %Plot 60 seconds
        end
        %Add labels etc to the plot
         xlim([-6 63])%xlim([0 60000]) 
         xticks([-6 0 10 20 30 40 50 60 63]) %xticks([-10 0 10 20 30 40 50 60 70]) %
         xticklabels({'-6','onset','10','20','30','40','50','offset','+3'}) %xticklabels({'-10','onset','10','20','30','40','50','offset','+10'})
         h = xline(0,'r'); %this line is added withthe -10 + 10 range
         h1 = xline(60,'r'); %this line is added with the -10 + 10 range
         xlabel('Time (s)');
         ylabel('Tremor power (log µV^2)');
         title([subject ' TOS - ' conditions{condition} ' time courses']);
        
        %Save the plot
        if strcmp(conf.output.save,'yes')
            fig_name = [ subject '_TOS_' conditions{condition} '_time_courses' ];
            saveas(gcf, fullfile(conf.output.figsubdir,fig_name), 'jpg');
            saveas(gcf, fullfile(conf.output.figsubdir,fig_name), 'fig');
        end
        
    else
        trial = 1;
        total_trial_num = indeces_trials( (condition-1)*4 + trial );
        power_over_time{condition}(trial,:) = log(squeeze(mtmconvol_freq_trials.powspctrm(total_trial_num,indexACC,new_indexFreq,1001:(length(mtmconvol_freq_trials.time)-1000))))';
        plot(mtmconvol_freq_trials.time(1001:(length(mtmconvol_freq_trials.time)-1000)), power_over_time{condition}(trial,:)) %Plot 60 seconds
        %Add labels etc to the plot
        xlim([-6 63])%xlim([0 60000]) 
        xticks([-6 0 10 20 30 40 50 60 63]) %xticks([-10 0 10 20 30 40 50 60 70]) %
        xticklabels({'-6','onset','10','20','30','40','50','offset','+3'}) %xticklabels({'-10','onset','10','20','30','40','50','offset','+10'})
        h = xline(0,'r'); %this line is added withthe -10 + 10 range
        h1 = xline(60,'r'); %this line is added with the -10 + 10 range
        xlabel('Time (s)');
        ylabel('Tremor power (log µV^2)');
        title([subject ' TOS - ' conditions{condition} ' time course']);
        
        %Save the plot
        if strcmp(conf.output.save,'yes')
            fig_name = [ subject '_TOS_' conditions{condition} '_time_course' ];
            saveas(gcf, fullfile(conf.output.figdir,fig_name), 'jpg');
            saveas(gcf, fullfile(conf.output.figdir,fig_name), 'fig');
        end
    end
end

%Calculate and plot average power time courses (TFR) per condition
for condition=1:length(conditions)-1
%     power_over_time{condition}(5,1:60000) = mean(power_over_time{condition}(1:4,1:60000),1);
    
    cfg = [];
    cfg.method      = 'mtmconvol';                           % Select method (choose 'mtmconvol')
    cfg.output      = 'pow';                                 % power
    cfg.taper       = 'hanning';                             % Windowing (because cut-off frequency), (Choose 'hanning' for low frequency)
    cfg.foi         = 2:0.5:16;                              % frequency range you are interested in (usually 1:0.5:20, make sure you at least include 3-8 Hz)
    nFoi            = length(cfg.foi);
    cfg.t_ftimwin   = repmat(2,1,nFoi);                      % Wavelet length (seconds; 1 wavelet per frequency). This is important also for your NaN in the hanning taper (which is 0.5*this)
    cfg.toi         = All_processed.time{1}(1):0.001:All_processed.time{1}(end);% timeline the TFR (resolution in seconds) ('orig': original resolution; 'timedat': one step specified under conf.prepemg.timedat;)
    cfg.pad         = 'maxperlen';                           % Padding (use 'maxperlen')
    cfg.trials      = indeces_trials((condition-1)*4+1:(condition-1)*4+4);
    cfg.keeptrials  = 'no';
    mtmconvol_freq_conditions = ft_freqanalysis(cfg, All_processed);   % The input data still contains 1s of data before and after the trial, to avoid NaNs at the beginning/end
    
    power_over_time{condition}(5,:) = log(squeeze(mtmconvol_freq_conditions.powspctrm(indexACC,new_indexFreq,1001:(length(mtmconvol_freq_trials.time)-1000))))';
    
    fig = figure('units','centimeters','outerposition',[0 0 25 20],'Color',[1 1 1]);
    plot(mtmconvol_freq_trials.time(1001:(length(mtmconvol_freq_trials.time)-1000)), power_over_time{condition}(5,:)) %Plot 60 seconds
    xlim([-6 63])%xlim([0 60000]) 
    xticks([-6 0 10 20 30 40 50 60 63]) %xticks([-10 0 10 20 30 40 50 60 70]) %
    xticklabels({'-6','0','10','20','30','40','50','60','+3'}) %xticklabels({'-10','onset','10','20','30','40','50','offset','+10'})
    h = xline(0,'r','Trial onset'); %this line is added withthe -10 + 10 range
    h1 = xline(60,'r','Trial offset'); %this line is added with the -10 + 10 range
    xlabel('Time (s)');
    ylabel('Tremor power (log µV^2)');
    title([subject ' TOS - ' conditions{condition} ' average time course']);
    
    if strcmp(conf.output.save,'yes')
        fig_name = [subject '_TOS_' conditions{condition} '_average_time_course'];
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'jpg');
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'fig');
    end
end

%Save the tremor power time courses (data)
if strcmp(conf.output.save,'yes')
    power_time_courses = power_over_time;
    save(fullfile(conf.output.datadir,[ subject '_TOS_log_power_trial_time_courses']),'power_time_courses');
end
% 
%% TOS 2nd level - step 1
% Create variable with group average pupil values for the TOS conditions
subject = '018';

cd(['/project/3024005.01/Analysis/Tremor/Results_correct_time_course/' subject '/Data']);

%threat-rest
load([subject '_TOS_log_power_trial_time_courses.mat']);
tremor_time_courses_threat_rest_n17(17,:) = power_time_courses{1,1}(5,:);
%threat-odd
tremor_time_courses_threat_odd_n17(17,:) = power_time_courses{1,2}(5,:);
%safe-rest
tremor_time_courses_safe_rest_n17(17,:) = power_time_courses{1,3}(5,:);
%safe-odd
tremor_time_courses_safe_odd_n17(17,:) = power_time_courses{1,4}(5,:);
%shock
tremor_time_courses_shock_n17(17,:) = power_time_courses{1,5}(1,:);

cd('/project/3024005.01/Analysis/Tremor/Results_correct_time_course/Group60s');
save('tremor_time_courses_threat_rest_n17.mat','tremor_time_courses_threat_rest_n17') 
save('tremor_time_courses_threat_odd_n17.mat','tremor_time_courses_threat_odd_n17') 
save('tremor_time_courses_safe_rest_n17.mat','tremor_time_courses_safe_rest_n17') 
save('tremor_time_courses_safe_odd_n17.mat','tremor_time_courses_safe_odd_n17') 
save('tremor_time_courses_shock_n17.mat','tremor_time_courses_shock_n17') 


%Calculate the average across 17 subjects for threat-rest
average_tremor_time_courses_threat_rest_n17 = mean(tremor_time_courses_threat_rest_n17);
save('average_tremor_time_courses_threat_rest_n17.mat','average_tremor_time_courses_threat_rest_n17') 

%Calculate the average across 17 subjects for threat-odd
average_tremor_time_courses_threat_odd_n17 = mean(tremor_time_courses_threat_odd_n17);
save('average_tremor_time_courses_threat_odd_n17.mat','average_tremor_time_courses_threat_odd_n17') 

%Calculate the average across 17 subjects for safe-rest
average_tremor_time_courses_safe_rest_n17 = mean(tremor_time_courses_safe_rest_n17);
save('average_tremor_time_courses_safe_rest_n17.mat','average_tremor_time_courses_safe_rest_n17') 

%Calculate the average across 17 subjects for safe-odd
average_tremor_time_courses_safe_odd_n17 = mean(tremor_time_courses_safe_odd_n17);
save('average_tremor_time_courses_safe_odd_n17.mat','average_tremor_time_courses_safe_odd_n17') 

%Calculate the average across 17 subjects for shock
average_tremor_time_courses_shock_n17 = mean(tremor_time_courses_shock_n17);
save('average_tremor_time_courses_shock_n17.mat','average_tremor_time_courses_shock_n17') 

%% time-course correlations at the trial (1,2,3,4), condition (thr_fix,thr_odd,safe_fix,safe_odd) and group (N = 17) level



%% %%% threat-rest trials %%%%%
subject = '001';
clearvars -except subject;

%for the tremor data
cd(['/project/3024005.01/Analysis/Tremor/Results_correct_time_course','/' subject '/Data']);
load([subject '_TOS_log_power_trial_time_courses.mat']); %load the coco % rest , TOS file

%for the pupil data
cd(['/project/3024005.01/Analysis/Pupil/Results','/' subject '/Data/0-60 data']);
load([subject '_threat_rest_pupil_values_trial_time_courses.mat']); %load the threat-rest file
load([subject '_threat_odd_pupil_values_trial_time_courses.mat']); %load the threat-odd file
load([subject '_safe_rest_pupil_values_trial_time_courses.mat']); %load the safe-rest file
load([subject '_safe_odd_pupil_values_trial_time_courses.mat']); %load the safe-odd file
load([subject '_shock_pupil_values_trial_time_courses.mat']); %load the shock file 
%% Trial 1
power_time_courses{1,1}(1,5000:65000); %tremor time-course trial 1 threat-rest
pupil_threat_rest_time_courses(1,:); %pupil time-course trial 1 threat-rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,1}(1,5000:65000)',pupil_threat_rest_time_courses(1,:)','Type','Pearson');
%Put the R and P value in an array which is saved later in this script
threat_rest_trial_correlations_array(:,:) = [R P];

%% Trial 2
power_time_courses{1,1}(2,5000:65000); %tremor time-course trial 1 threat-rest
pupil_threat_rest_time_courses(2,:); %pupil time-course trial 1 threat-rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,1}(2,5000:65000)',pupil_threat_rest_time_courses(2,:)','Type','Pearson');
%Put the R and P value in an array which is saved later in this script
threat_rest_trial_correlations_array(2,:) = [R P];

%% Trial 3
power_time_courses{1,1}(3,5000:65000); %tremor time-course trial 1 threat-rest
pupil_threat_rest_time_courses(3,:); %pupil time-course trial 1 threat-rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,1}(3,5000:65000)',pupil_threat_rest_time_courses(3,:)','Type','Pearson');
%Put the R and P value in an array which is saved later in this script
threat_rest_trial_correlations_array(3,:) = [R P];

%% Trial 4
power_time_courses{1,1}(4,5000:65000); %tremor time-course trial 1 threat-rest
pupil_threat_rest_time_courses(4,:); %pupil time-course trial 1 threat-rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,1}(4,5000:65000)',pupil_threat_rest_time_courses(4,:)','Type','Pearson');
%Put the R and P value in an array which is saved later in this script
threat_rest_trial_correlations_array(4,:) = [R P];

%% Save the trial correlations in an array
cd(['/project/3024005.01/Analysis/Tremor/Results_correct_time_course','/' subject '/Data']);
save('threat_rest_trial_correlations_array.mat','threat_rest_trial_correlations_array')

%% Save the condition values per subject in 1 array
cd('/project/3024005.01/Analysis/Tremor/Results_correct_time_course/Group60s');
load('avg_corr_threat_rest.mat')
avg_corr_threat_rest(:,17) = threat_rest_trial_correlations_array(:,1);
save('avg_corr_threat_rest.mat','avg_corr_threat_rest')

%% calculate the mean per per subject for threat-rest
avg_corr_threat_rest_n17 = mean(avg_corr_threat_rest(:,:))';
save('avg_corr_threat_rest_n17.mat','avg_corr_threat_rest_n17')
%% %%% threat-odd trials %%%
subject = '018';
clearvars -except subject;

%for the tremor data
cd(['/project/3024005.01/Analysis/Tremor/Results_correct_time_course','/' subject '/Data']);
load([subject '_TOS_log_power_trial_time_courses.mat']); %load the coco % rest file

%for the pupil data
cd(['/project/3024005.01/Analysis/Pupil/Results','/' subject '/Data/0-60 data']);
load([subject '_threat_rest_pupil_values_trial_time_courses.mat']); %load the threat-rest file
load([subject '_threat_odd_pupil_values_trial_time_courses.mat']); %load the threat-odd file
load([subject '_safe_rest_pupil_values_trial_time_courses.mat']); %load the safe-rest file
load([subject '_safe_odd_pupil_values_trial_time_courses.mat']); %load the safe-odd file
load([subject '_shock_pupil_values_trial_time_courses.mat']); %load the shock file

%% Trial 1
power_time_courses{1,2}(1,5000:65000); %tremor time-course trial 1 threat-rest
pupil_threat_odd_time_courses(1,:); %pupil time-course trial 1 threat-rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,2}(1,5000:65000)',pupil_threat_odd_time_courses(1,:)','Type','Pearson');
%Put the R and P value in an array which is saved later in this script
threat_odd_trial_correlations_array(:,:) = [R P];

%% Trial 2
power_time_courses{1,2}(2,5000:65000); %tremor time-course trial 1 threat-rest
pupil_threat_odd_time_courses(2,:); %pupil time-course trial 1 threat-rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,2}(2,5000:65000)',pupil_threat_odd_time_courses(2,:)','Type','Pearson');
%Put the R and P value in an array which is saved later in this script
threat_odd_trial_correlations_array(2,:) = [R P];

%% Trial 3
power_time_courses{1,2}(3,5000:65000); %tremor time-course trial 1 threat-rest
pupil_threat_odd_time_courses(3,:); %pupil time-course trial 1 threat-rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,2}(3,5000:65000)',pupil_threat_odd_time_courses(3,:)','Type','Pearson');
%Put the R and P value in an array which is saved later in this script
threat_odd_trial_correlations_array(3,:) = [R P];

%% Trial 4
power_time_courses{1,2}(4,:); %tremor time-course trial 1 threat-rest
pupil_threat_odd_time_courses(4,1:60000); %pupil time-course trial 1 threat-rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,2}(4,5000:65000)',pupil_threat_odd_time_courses(4,:)','Type','Pearson');
%Put the R and P value in an array which is saved later in this script
threat_odd_trial_correlations_array(4,:) = [R P];

%% Save the trial correlations in an array
cd(['/project/3024005.01/Analysis/Tremor/Results_correct_time_course','/' subject '/Data']);
save('threat_odd_trial_correlations_array.mat','threat_odd_trial_correlations_array')

%% Save the condition values per subject in 1 array
cd('/project/3024005.01/Analysis/Tremor/Results_correct_time_course/Group60s');
load('avg_corr_threat_odd.mat')
avg_corr_threat_odd(:,17) = threat_odd_trial_correlations_array(:,1);
save('avg_corr_threat_odd.mat','avg_corr_threat_odd')

%% calculate the mean per per subject for threat-odd
avg_corr_threat_odd_n17 = mean(avg_corr_threat_odd(:,:))';
save('avg_corr_threat_odd_n17.mat','avg_corr_threat_odd_n17')

%% %%% safe-rest trials %%%
subject = '018';
clearvars -except subject;

%for the tremor data
cd(['/project/3024005.01/Analysis/Tremor/Results_correct_time_course','/' subject '/Data']);
load([subject '_TOS_log_power_trial_time_courses.mat']); %load the coco % rest file

%for the pupil data
cd(['/project/3024005.01/Analysis/Pupil/Results','/' subject '/Data/0-60 data']);
load([subject '_threat_rest_pupil_values_trial_time_courses.mat']); %load the threat-rest file
load([subject '_threat_odd_pupil_values_trial_time_courses.mat']); %load the threat-odd file
load([subject '_safe_rest_pupil_values_trial_time_courses.mat']); %load the safe-rest file
load([subject '_safe_odd_pupil_values_trial_time_courses.mat']); %load the safe-odd file
load([subject '_shock_pupil_values_trial_time_courses.mat']); %load the shock file
%% Trial 1
power_time_courses{1,3}(1,5000:65000); %tremor time-course trial 1 threat-rest
pupil_safe_rest_time_courses(1,1:60000); %pupil time-course trial 1 threat-rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,3}(1,5000:65000)',pupil_safe_rest_time_courses(1,:)','Type','Pearson');
%Put the R and P value in an array which is saved later in this script
safe_rest_trial_correlations_array(:,:) = [R P];

%% Trial 2
power_time_courses{1,3}(2,5000:65000); %tremor time-course trial 1 threat-rest
pupil_safe_rest_time_courses(2,:); %pupil time-course trial 1 threat-rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,3}(2,5000:65000)',pupil_safe_rest_time_courses(2,:)','Type','Pearson');
%Put the R and P value in an array which is saved later in this script
safe_rest_trial_correlations_array(2,:) = [R P];

%% Trial 3
power_time_courses{1,3}(3,5000:65000); %tremor time-course trial 1 threat-rest
pupil_safe_rest_time_courses(3,:); %pupil time-course trial 1 threat-rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,3}(3,5000:65000)',pupil_safe_rest_time_courses(3,:)','Type','Pearson');
%Put the R and P value in an array which is saved later in this script
safe_rest_trial_correlations_array(3,:) = [R P];

%% Trial 4
power_time_courses{1,3}(4,5000:65000); %tremor time-course trial 1 threat-rest
pupil_safe_rest_time_courses(4,:); %pupil time-course trial 1 threat-rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,3}(4,5000:65000)',pupil_safe_rest_time_courses(4,:)','Type','Pearson');
%Put the R and P value in an array which is saved later in this script
safe_rest_trial_correlations_array(4,:) = [R P];

%% Save the trial correlations in an array
cd(['/project/3024005.01/Analysis/Tremor/Results_correct_time_course','/' subject '/Data']);
save('safe_rest_trial_correlations_array.mat','safe_rest_trial_correlations_array')

%% Save the condition values per subject in 1 array
cd('/project/3024005.01/Analysis/Tremor/Results_correct_time_course/Group60s');
load('avg_corr_safe_rest.mat')
avg_corr_safe_rest(:,17) = safe_rest_trial_correlations_array(:,1);
save('avg_corr_safe_rest.mat','avg_corr_safe_rest')

%% calculate the mean per per subject for safe-rest
avg_corr_safe_rest_n17 = mean(avg_corr_safe_rest(:,:))';
save('avg_corr_safe_rest_n17.mat','avg_corr_safe_rest_n17')


%% %%% safe-odd trials %%%
subject = '018';
clearvars -except subject;

%for the tremor data
cd(['/project/3024005.01/Analysis/Tremor/Results_correct_time_course','/' subject '/Data']);
load([subject '_TOS_log_power_trial_time_courses.mat']); %load the coco % rest file

%for the pupil data
cd(['/project/3024005.01/Analysis/Pupil/Results','/' subject '/Data/0-60 data']);
load([subject '_threat_rest_pupil_values_trial_time_courses.mat']); %load the threat-rest file
load([subject '_threat_odd_pupil_values_trial_time_courses.mat']); %load the threat-odd file
load([subject '_safe_rest_pupil_values_trial_time_courses.mat']); %load the safe-rest file
load([subject '_safe_odd_pupil_values_trial_time_courses.mat']); %load the safe-odd file
load([subject '_shock_pupil_values_trial_time_courses.mat']); %load the shock file
%% Trial 1
power_time_courses{1,4}(1,5000:65000); %tremor time-course trial 1 threat-rest
pupil_safe_odd_time_courses(1,:); %pupil time-course trial 1 threat-rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,4}(1,5000:65000)',pupil_safe_odd_time_courses(1,:)','Type','Pearson');
%Put the R and P value in an array which is saved later in this script
safe_odd_trial_correlations_array(:,:) = [R P];

%% Trial 2
power_time_courses{1,4}(2,5000:65000); %tremor time-course trial 1 threat-rest
pupil_safe_odd_time_courses(2,:); %pupil time-course trial 1 threat-rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,4}(2,5000:65000)',pupil_safe_odd_time_courses(2,:)','Type','Pearson');
%Put the R and P value in an array which is saved later in this script
safe_odd_trial_correlations_array(2,:) = [R P];

%% Trial 3
power_time_courses{1,4}(3,5000:65000); %tremor time-course trial 1 threat-rest
pupil_safe_odd_time_courses(3,:); %pupil time-course trial 1 threat-rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,4}(3,5000:65000)',pupil_safe_odd_time_courses(3,:)','Type','Pearson');
%Put the R and P value in an array which is saved later in this script
safe_odd_trial_correlations_array(3,:) = [R P];

%% Trial 4
power_time_courses{1,4}(4,5000:65000); %tremor time-course trial 1 threat-rest
pupil_safe_odd_time_courses(4,:); %pupil time-course trial 1 threat-rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,4}(4,5000:65000)',pupil_safe_odd_time_courses(4,:)','Type','Pearson');
%Put the R and P value in an array which is saved later in this script
safe_odd_trial_correlations_array(4,:) = [R P];

%% Save the trial correlations in an array
cd(['/project/3024005.01/Analysis/Tremor/Results_correct_time_course','/' subject '/Data']);
save('safe_odd_trial_correlations_array.mat','safe_odd_trial_correlations_array')

%% Save the condition values per subject in 1 array
cd('/project/3024005.01/Analysis/Tremor/Results_correct_time_course/Group60s');
load('avg_corr_safe_odd.mat')
avg_corr_safe_odd(:,17) = safe_odd_trial_correlations_array(:,1);
save('avg_corr_safe_odd.mat','avg_corr_safe_odd')

%% calculate the mean per per subject for safe-odd
avg_corr_safe_odd_n17 = mean(avg_corr_safe_odd(:,:))';
save('avg_corr_safe_odd_n17.mat','avg_corr_safe_odd_n17')

%% %%% shock trials %%%
subject = '018';
clearvars -except subject;

%for the tremor data
cd(['/project/3024005.01/Analysis/Tremor/Results_correct_time_course','/' subject '/Data']);
load([subject '_TOS_log_power_trial_time_courses.mat']); %load the coco % rest file

%for the pupil data
cd(['/project/3024005.01/Analysis/Pupil/Results','/' subject '/Data/0-60 data']);
load([subject '_threat_rest_pupil_values_trial_time_courses.mat']); %load the threat-rest file
load([subject '_threat_odd_pupil_values_trial_time_courses.mat']); %load the threat-odd file
load([subject '_safe_rest_pupil_values_trial_time_courses.mat']); %load the safe-rest file
load([subject '_safe_odd_pupil_values_trial_time_courses.mat']); %load the safe-odd file
load([subject '_shock_pupil_values_trial_time_courses.mat']); %load the shock file
%% Trial 1
power_time_courses{1,5}(1,5000:65000); %tremor time-course trial 1 threat-rest
pupil_shock_time_courses(1,:); %pupil time-course trial 1 threat-rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,5}(1,5000:65000)',pupil_shock_time_courses(1,:)','Type','Pearson');

%Put the R and P value in an array which is saved later in this script
shock_trial_correlations_array(:,:) = [R P];

%% Save the trial correlations in an array
cd(['/project/3024005.01/Analysis/Tremor/Results_correct_time_course','/' subject '/Data']);
save('shock_trial_correlations_array.mat','shock_trial_correlations_array')

%% Save the condition values per subject in 1 array
cd('/project/3024005.01/Analysis/Tremor/Results_correct_time_course/Group60s');
load('avg_corr_shock.mat')
avg_corr_shock(:,17) = shock_trial_correlations_array(:,1);
save('avg_corr_shock.mat','avg_corr_shock')
avg_corr_shock_n17 = avg_corr_shock;
save('avg_corr_shock_n17.mat','avg_corr_shock_n17')

%% correlate the VAS self-report tremor to stress with delta coco-rest, delta threat-safe
%Perform the correlation
[R,P] = corr(VAS_stress_tremor,delta_coco_rest,'Type','Pearson'); %R = 0.0636, P = 0.8085
[R,P] = corr(VAS_stress_tremor,delta_threat_safe,'Type','Pearson'); 

h = figure ('Color', [1 1 1]);
s1 = plot(VAS_stress_tremor, delta_threat_safe, 'k+');
set(s1, 'MarkerSize', 8, 'LineWidth', 2);
%%% regression line
hold on
l = lsline ;
set(l,'LineWidth', 2)
%%% axis display
xlim([0 100])
xlabel('self-report tremor response to stress', 'FontSize', 20)
ylabel('delta threat-safe', 'FontSize', 20)
set(gca, 'FontSize', 20, 'YMinorTick','on','XMinorTick','on')
r = corrcoef(VAS_stress_tremor(:, 1), delta_threat_safe(:, 1));
disp(r(1,2));

str=['r= ',num2str(r(1,2))]
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
p = 0.8085;
p_val = legend([p],{'pval'});
  set(p, 'Position', [.65 .16 .11 .08], 'Color', [1 1 1],'FontSize',14);
%% for barplotting with avnBarplotscatter.m

%for the 4 TOS conditions
avg_corr_TOS_cell = num2cell(avg_corr_threat_rest,[1 2]);
avg_corr_TOS_cell{1,1}(:,2) = (avg_corr_safe_rest);

avg_corr_TOS_cell{1,2}(:,1) = (avg_corr_threat_odd);
avg_corr_TOS_cell{1,2}(:,2) = (avg_corr_safe_odd);
cd('/project/3024005.01/Analysis/Tremor/NewResults/Group60s/0-60 time-courses')
save('avg_corr_TOS_cell.mat', 'avg_corr_TOS_cell')

% 
% %for the shock condition
% %calculate the mean correlation for the 4 TOS conditions and put this in 1
% %column in 1 cell array
% avg_corr_allTOS_cell{1,1}(1,1) = mean([ avg_corr_TOS_cell{1,1}(1,:) avg_corr_TOS_cell{1,2}(1,:)]);
% avg_corr_allTOS_cell{1,1}(2,1) = mean([ avg_corr_TOS_cell{1,1}(2,:) avg_corr_TOS_cell{1,2}(2,:)]);
% avg_corr_allTOS_cell{1,1}(3,1) = mean([ avg_corr_TOS_cell{1,1}(3,:) avg_corr_TOS_cell{1,2}(3,:)]);
% avg_corr_allTOS_cell{1,1}(4,1) = mean([ avg_corr_TOS_cell{1,1}(4,:) avg_corr_TOS_cell{1,2}(4,:)]);
% avg_corr_allTOS_cell{1,1}(5,1) = mean([ avg_corr_TOS_cell{1,1}(5,:) avg_corr_TOS_cell{1,2}(5,:)]);
% avg_corr_allTOS_cell{1,1}(6,1) = mean([ avg_corr_TOS_cell{1,1}(6,:) avg_corr_TOS_cell{1,2}(6,:)]);
% avg_corr_allTOS_cell{1,1}(7,1) = mean([ avg_corr_TOS_cell{1,1}(7,:) avg_corr_TOS_cell{1,2}(7,:)]);
% avg_corr_allTOS_cell{1,1}(8,1) = mean([ avg_corr_TOS_cell{1,1}(8,:) avg_corr_TOS_cell{1,2}(8,:)]);
% avg_corr_allTOS_cell{1,1}(9,1) = mean([ avg_corr_TOS_cell{1,1}(9,:) avg_corr_TOS_cell{1,2}(9,:)]);
% avg_corr_allTOS_cell{1,1}(10,1) = mean([ avg_corr_TOS_cell{1,1}(10,:) avg_corr_TOS_cell{1,2}(10,:)]);
% avg_corr_allTOS_cell{1,1}(11,1) = mean([ avg_corr_TOS_cell{1,1}(11,:) avg_corr_TOS_cell{1,2}(11,:)]);
% avg_corr_allTOS_cell{1,1}(12,1) = mean([ avg_corr_TOS_cell{1,1}(12,:) avg_corr_TOS_cell{1,2}(12,:)]);
% avg_corr_allTOS_cell{1,1}(13,1) = mean([ avg_corr_TOS_cell{1,1}(13,:) avg_corr_TOS_cell{1,2}(13,:)]);
% avg_corr_allTOS_cell{1,1}(14,1) = mean([ avg_corr_TOS_cell{1,1}(14,:) avg_corr_TOS_cell{1,2}(14,:)]);
% avg_corr_allTOS_cell{1,1}(15,1) = mean([ avg_corr_TOS_cell{1,1}(15,:) avg_corr_TOS_cell{1,2}(15,:)]);
% avg_corr_allTOS_cell{1,1}(16,1) = mean([ avg_corr_TOS_cell{1,1}(16,:) avg_corr_TOS_cell{1,2}(16,:)]);
% avg_corr_allTOS_cell{1,1}(17,1) = mean([ avg_corr_TOS_cell{1,1}(17,:) avg_corr_TOS_cell{1,2}(17,:)]);
% 
% cd('/project/3024005.01/Analysis/Tremor/NewResults/Group60s/0-60 time-courses')
% save('avg_corr_allTOS_cell.mat', 'avg_corr_allTOS_cell')
% 
% %put the shock double in the avg_corr_allTOS_cell 
% avg_corr_allTOS_cell{1,2} = avg_corr_shock;
% save('avg_corr_allTOS_cell.mat', 'avg_corr_allTOS_cell')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Correlation plot %%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% Plot the correlation for threat-rest at the group level
% threat_rest_group_correlations_figure = figure('units','centimeters','outerposition',[0 0 25 20],'Color',[1 1 1]);
% s1 = plot(power_time_courses_60s_group_thr_rest(:,1),pupil_time_courses_60s_group_thr_rest(:,1), '+', 'MarkerFaceColor', 'k');
%     ylabel('Pupil diameter','FontSize', 12);
%     xlabel ('Tremor power','FontSize', 12);
%     title('Threat-rest correlation plot (N = 17)');
%     set(s1, 'MarkerSize', 8, 'LineWidth', 2);
% %Add the correlation coefficient
%     disp(RHO(1,1));
%     str=['RHO= ',num2str(RHO(1,1))]
%     RHOVAL = text(min(get(gca, 'xlim')+1), max(get(gca, 'ylim')-100), str);
%     set(RHOVAL, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left')
% %Add the p-value
%     disp(PVAL(1,1));
%     str=['P= ',num2str(PVAL(1,1))]
%     SIGVAL = text(min(get(gca, 'xlim')+1), max(get(gca, 'ylim')-150), str);
%     set(SIGVAL, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left')    
% %Add a regression line
%     hold on
%     l = lsline ;
%     set(l,'LineWidth', 2)
% %Save figures in .fig and .jpg format
%     saveas(threat_rest_group_correlations_figure,'threat_rest_group_correlations_figure.jpg')
%     savefig('threat_rest_group_correlations_figure.fig')
 
 %%%%%%%%%%%%%%%%%%%%%%%%%% OLD correlation scripts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

%subject = '001';
%clearvars -except subject;

%for the tremor data
%cd(['/project/3024005.01/Analysis/Tremor/NewResults','/' subject '/Data']);
%load([subject '_TOS_log_power_trial_time_courses.mat']); %load the all the TOS conditions file

%for the pupil data
%cd(['/project/3024005.01/Analysis/Pupil/Results','/' subject '/Data']);
%load([subject '_threat_rest_pupil_values_trial_time_courses.mat']); %load the threat-rest file
%load([subject '_threat_odd_pupil_values_trial_time_courses.mat']); %load the threat-odd file
%load([subject '_safe_rest_pupil_values_trial_time_courses.mat']); %load the safe-rest file
%load([subject '_safe_odd_pupil_values_trial_time_courses.mat']); %load the safe-odd file
%load([subject '_shock_pupil_values_trial_time_courses.mat']); %load the shock file


%%%%% threat-rest trials %%%%%
%power_time_courses{1,1}(1,:); %tremor time-course trial 1 threat-rest
%pupil_threat_rest_time_courses(1,1:60000); %pupil time-course trial 1 threat-rest
%[RHO,PVAL] = corr(power_time_courses{1,1}(1,:)',pupil_threat_rest_time_courses(1,1:60000)','Type','Spearman');

%threat_rest_trial_correlations_array(:,:) = [RHO PVAL];

%power_time_courses{1,1}(2,:); %tremor time-course trial 2 threat-rest
%pupil_threat_rest_time_courses(2,1:60000); %pupil time-course trial 2 threat-rest
%[RHO,PVAL] = corr(power_time_courses{1,1}(2,:)',pupil_threat_rest_time_courses(2,1:60000)','Type','Spearman');

%threat_rest_trial_correlations_array(2,:) = [RHO PVAL];

%power_time_courses{1,1}(3,:); %tremor time-course trial 3 threat-rest
%pupil_threat_rest_time_courses(3,1:60000); %pupil time-course trial 3 threat-rest
%[RHO,PVAL] = corr(power_time_courses{1,1}(3,:)',pupil_threat_rest_time_courses(3,1:60000)','Type','Spearman');

%threat_rest_trial_correlations_array(3,:) = [RHO PVAL];

%power_time_courses{1,1}(4,:); %tremor time-course trial 4 threat-rest
%pupil_threat_rest_time_courses(4,1:60000); %pupil time-course  trial 4 threat-rest
%[RHO,PVAL] = corr(power_time_courses{1,1}(4,:)',pupil_threat_rest_time_courses(4,1:60000)','Type','Spearman');

%threat_rest_trial_correlations_array(4,:) = [RHO PVAL];

%cd(['/project/3024005.01/Analysis/Tremor/NewResults','/' subject '/Data']);
%save('threat_rest_trial_correlations_array.mat','threat_rest_trial_correlations_array')

%%%%% threat-odd trials %%%%%
%power_time_courses{1,2}(1,:); %tremor time-course trial 1 threat-odd
%pupil_threat_odd_time_courses(1,1:60000); %pupil time-course trial 1 threat-odd
%[RHO,PVAL] = corr(power_time_courses{1,2}(1,:)',pupil_threat_odd_time_courses(1,1:60000)','Type','Spearman');

%threat_odd_trial_correlations_array(:,:) = [RHO PVAL];

%power_time_courses{1,2}(2,:); %tremor time-course trial 2 threat-rest
%pupil_threat_odd_time_courses(2,1:60000); %pupil time-course trial 2 threat-rest
%[RHO,PVAL] = corr(power_time_courses{1,2}(2,:)',pupil_threat_odd_time_courses(2,1:60000)','Type','Spearman');

%threat_odd_trial_correlations_array(2,:) = [RHO PVAL];

%power_time_courses{1,2}(3,:); %tremor time-course trial 3 threat-rest
%pupil_threat_odd_time_courses(3,1:60000); %pupil time-course trial 3 threat-rest
%[RHO,PVAL] = corr(power_time_courses{1,2}(3,:)',pupil_threat_odd_time_courses(3,1:60000)','Type','Spearman');

%threat_odd_trial_correlations_array(3,:) = [RHO PVAL];

%power_time_courses{1,2}(4,:); %tremor time-course trial 4 threat-rest
%pupil_threat_odd_time_courses(4,1:60000); %pupil time-course  trial 4 threat-rest
%[RHO,PVAL] = corr(power_time_courses{1,2}(4,:)',pupil_threat_odd_time_courses(4,1:60000)','Type','Spearman');

%threat_odd_trial_correlations_array(4,:) = [RHO PVAL];

%cd(['/project/3024005.01/Analysis/Tremor/NewResults','/' subject '/Data']);
%save('threat_odd_trial_correlations_array.mat','threat_odd_trial_correlations_array')

%%%%% safe-rest trials %%%%%
%power_time_courses{1,3}(1,:); %tremor time-course trial 1 threat-odd
%pupil_safe_rest_time_courses(1,1:60000); %pupil time-course trial 1 threat-odd
%[RHO,PVAL] = corr(power_time_courses{1,3}(1,:)',pupil_safe_rest_time_courses(1,1:60000)','Type','Spearman');

%safe_rest_trial_correlations_array(:,:) = [RHO PVAL];

%power_time_courses{1,3}(2,:); %tremor time-course trial 2 threat-rest
%pupil_safe_rest_time_courses(2,1:60000); %pupil time-course trial 2 threat-rest
%[RHO,PVAL] = corr(power_time_courses{1,3}(2,:)',pupil_safe_rest_time_courses(2,1:60000)','Type','Spearman');

%safe_rest_trial_correlations_array(2,:) = [RHO PVAL];

%power_time_courses{1,3}(3,:); %tremor time-course trial 3 threat-rest
%pupil_safe_rest_time_courses(3,1:60000); %pupil time-course trial 3 threat-rest
%[RHO,PVAL] = corr(power_time_courses{1,3}(3,:)',pupil_safe_rest_time_courses(3,1:60000)','Type','Spearman');

%safe_rest_trial_correlations_array(3,:) = [RHO PVAL];

%power_time_courses{1,3}(4,:); %tremor time-course trial 4 threat-rest
%pupil_safe_rest_time_courses(4,1:60000); %pupil time-course  trial 4 threat-rest
%[RHO,PVAL] = corr(power_time_courses{1,3}(4,:)',pupil_safe_rest_time_courses(4,1:60000)','Type','Spearman');

%safe_rest_trial_correlations_array(4,:) = [RHO PVAL];

%cd(['/project/3024005.01/Analysis/Tremor/NewResults','/' subject '/Data']);
%save('safe_rest_trial_correlations_array.mat','safe_rest_trial_correlations_array')

%%%%% safe-odd trials %%%%%
%power_time_courses{1,4}(1,:); %tremor time-course trial 1 threat-odd
%pupil_safe_odd_time_courses(1,1:60000); %pupil time-course trial 1 threat-odd
%[RHO,PVAL] = corr(power_time_courses{1,4}(1,:)',pupil_safe_odd_time_courses(1,1:60000)','Type','Spearman');

%safe_odd_trial_correlations_array(:,:) = [RHO PVAL];

%power_time_courses{1,4}(2,:); %tremor time-course trial 2 threat-rest
%pupil_safe_odd_time_courses(2,1:60000); %pupil time-course trial 2 threat-rest
%[RHO,PVAL] = corr(power_time_courses{1,4}(2,:)',pupil_safe_odd_time_courses(2,1:60000)','Type','Spearman');

%safe_odd_trial_correlations_array(2,:) = [RHO PVAL];

%power_time_courses{1,4}(3,:); %tremor time-course trial 3 threat-rest
%pupil_safe_odd_time_courses(3,1:60000); %pupil time-course trial 3 threat-rest
%[RHO,PVAL] = corr(power_time_courses{1,4}(3,:)',pupil_safe_odd_time_courses(3,1:60000)','Type','Spearman');

%safe_odd_trial_correlations_array(3,:) = [RHO PVAL];

%power_time_courses{1,4}(4,:); %tremor time-course trial 4 threat-rest
%pupil_safe_odd_time_courses(4,1:60000); %pupil time-course  trial 4 threat-rest
%[RHO,PVAL] = corr(power_time_courses{1,4}(4,:)',pupil_safe_odd_time_courses(4,1:60000)','Type','Spearman');

%safe_odd_trial_correlations_array(4,:) = [RHO PVAL];

%cd(['/project/3024005.01/Analysis/Tremor/NewResults','/' subject '/Data']);
%save('safe_odd_trial_correlations_array.mat','safe_odd_trial_correlations_array')

%%%%% shock trial %%%%%
%power_time_courses{1,5}(1,:); %tremor time-course shock trial
%pupil_shock_time_courses(1,1:60000); %pupil time-course shock trial
%[RHO,PVAL] = corr(power_time_courses{1,5}(1,:)',pupil_shock_time_courses(1,1:60000)','Type','Spearman');

%shock_trial_correlations_array(:,:) = [RHO PVAL];

%cd(['/project/3024005.01/Analysis/Tremor/NewResults','/' subject '/Data']);
%save('shock_trial_correlations_array.mat','shock_trial_correlations_array')

%%% CONDITIONS %%%%
%cd(['/project/3024005.01/Analysis/Pupil/Results','/' subject '/Data']);
%load([subject '_average_threat_rest_pupil_values_trial_time_courses.mat']); %load the threat-rest file
%load([subject '_average_threat_odd_pupil_values_trial_time_courses.mat']); %load the threat-odd file
%load([subject '_average_safe_rest_pupil_values_trial_time_courses.mat']); %load the safe-rest file
%load([subject '_average_safe_odd_pupil_values_trial_time_courses.mat']); %load the safe-odd file
%load([subject '_average_shock_pupil_values_trial_time_courses.mat']); %load the shock file

%%% threat-rest condition %%%%
%power_time_courses{1,1}(5,:); %average tremor time-course threat-rest condition
%average_pupil_threat_rest_time_courses(1,1:60000); %average pupil time-course threat-rest condition
%[RHO,PVAL] = corr(power_time_courses{1,1}(5,:)',average_pupil_threat_rest_time_courses(1,1:60000)','Type','Spearman');

%threat_rest_condition_correlations_array(:,:) = [RHO PVAL];

%cd(['/project/3024005.01/Analysis/Tremor/NewResults','/' subject '/Data']');
%save('threat_rest_condition_correlations_array.mat','threat_rest_condition_correlations_array')

%%% threat-odd condition %%%%
%power_time_courses{1,2}(5,:); %average tremor time-course threat-rest condition
%average_pupil_threat_odd_time_courses(1,1:60000); %average pupil time-course threat-rest condition
%[RHO,PVAL] = corr(power_time_courses{1,2}(5,:)',average_pupil_threat_odd_time_courses(1,1:60000)','Type','Spearman');

%threat_odd_condition_correlations_array(:,:) = [RHO PVAL];

%cd(['/project/3024005.01/Analysis/Tremor/NewResults','/' subject '/Data']');
%save('threat_odd_condition_correlations_array.mat','threat_odd_condition_correlations_array')

%%% safe-rest condition %%%%
%power_time_courses{1,3}(5,:); %average tremor time-course threat-rest condition
%average_pupil_safe_rest_time_courses(1,1:60000); %average pupil time-course threat-rest condition
%[RHO,PVAL] = corr(power_time_courses{1,3}(5,:)',average_pupil_safe_rest_time_courses(1,1:60000)','Type','Spearman');

%safe_rest_condition_correlations_array(:,:) = [RHO PVAL];

%cd(['/project/3024005.01/Analysis/Tremor/NewResults','/' subject '/Data']');
%save('safe_rest_condition_correlations_array.mat','safe_rest_condition_correlations_array')

%%% safe-odd condition %%%%
%power_time_courses{1,4}(5,:); %average tremor time-course threat-rest condition
%average_pupil_safe_odd_time_courses(1,1:60000); %average pupil time-course threat-rest condition
%[RHO,PVAL] = corr(power_time_courses{1,4}(5,:)',average_pupil_safe_odd_time_courses(1,1:60000)','Type','Spearman');

%safe_odd_condition_correlations_array(:,:) = [RHO PVAL];

%cd(['/project/3024005.01/Analysis/Tremor/NewResults','/' subject '/Data']');
%save('safe_odd_condition_correlations_array.mat','safe_odd_condition_correlations_array')

%%% shock condition %%%%
%power_time_courses{1,5}(1,:); %average tremor time-course threat-rest condition
%average_pupil_shock_time_courses(1,1:60000); %average pupil time-course threat-rest condition
%[RHO,PVAL] = corr(power_time_courses{1,5}(1,:)',average_pupil_shock_time_courses(1,1:60000)','Type','Spearman');

%shock_condition_correlations_array(:,:) = [RHO PVAL];

%cd(['/project/3024005.01/Analysis/Tremor/NewResults','/' subject '/Data']');
%save('shock_condition_correlations_array.mat','shock_condition_correlations_array')

end
