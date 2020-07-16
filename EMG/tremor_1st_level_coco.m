function tremor_1st_level_coco(subject)

subject = '018';
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
conf.input.eegfile = pf_findfile(conf.input.eegdir,['/' subject '/&/coco/&/.eeg/'],'fullfile');      

% conf.output.maindir = [ '/project/3024005.01/Analysis/Tremor/NewResults/' subject ];
% conf.output.figdir = [ '/project/3024005.01/Analysis/Tremor/NewResults/' subject '/Figures' ];
% conf.output.figsubdir = [ '/project/3024005.01/Analysis/Tremor/NewResults/' subject '/Figures/SingleTrials' ];
% conf.output.datadir = [ '/project/3024005.01/Analysis/Tremor/NewResults/' subject '/Data' ];
conf.output.maindir = [ '/project/3024005.01/Analysis/Tremor/Results_correct_time_course/' subject ];
conf.output.figdir = [ '/project/3024005.01/Analysis/Tremor/Results_correct_time_course/' subject '/Figures' ];
conf.output.figsubdir = [ '/project/3024005.01/Analysis/Tremor/Results_correct_time_course/' subject '/Figures/SingleTrials' ];
conf.output.datadir = [ '/project/3024005.01/Analysis/Tremor/Results_correct_time_course/' subject '/Data' ];

temp = fieldnames(conf.output);
for i=1:numel(temp)
    if ~exist( getfield(conf.output,temp{i}) )
        mkdir( getfield(conf.output,temp{i}) );
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
cfg.trialdef.eventvalue = {'S 11','S 12'};    
cfg.trialdef.prestim    = 6; %7 for -6 %also read in 1s before and after trial for mtmconvol later (similar to padding around trial)
cfg.trialdef.poststim   = 64; %64 for +3 %31 for 30 sec trials
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
for i=1:size(EMGdata_no_artifacts.trial,2)
    EMGdata_no_artifacts.trial{i}(isnan(EMGdata_no_artifacts.trial{i}))=0;
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
cfg.toilim = [0 60]; %Length of the data for fast fourier transform
segdata = ft_redefinetrial(cfg,All_processed); %for analysis using fft drop second before and after the trial
cfg = [];
cfg.length  = 5; 
cfg.overlap = 0;
cfg.trials  = 'all';
segdata = ft_redefinetrial(cfg,segdata); % cut data into 5s segments

cfg = [];
cfg.method      = 'mtmfft'; % mtmfft analyses an entire spectrum for the entire data length, implements multitaper frequency transformation
cfg.foi         = 2:0.2:16;     % frequency range you are interested in. Freq resolution depends on the length of the time window (1/T).
cfg.taper       = 'hanning';
cfg.t_ftimwin   = ones(length(cfg.foi),1)*5; %*5
cfg.toi         = 'all';  %all trials (trial now =5s segment)
cfg.keeptrials  = 'no'; %return average (default = no)
fft_freq        = ft_freqanalysis(cfg,segdata);

%% Plot the powerspectrum of all ACC channels & the EMG channels (channel 5 - 7) of the most affected arm (channel 1 & 2)
fig = figure('units','centimeters','outerposition',[0 0 20 20],'Color',[1 1 1]);
hold on;
plot(fft_freq.freq, fft_freq.powspctrm(5,:)) %acc_x
plot(fft_freq.freq, fft_freq.powspctrm(6,:)) %acc_y
plot(fft_freq.freq, fft_freq.powspctrm(7,:)) %acc_z
legend('ACC_x','ACC_y','ACC_z')
xlabel('Frequency (Hz)');
ylabel('Absolute tremor power (µV^2)');
title([ subject ' coco - ACC powerspectrum']);
if strcmp(conf.output.save,'yes')
    fig_name = [ subject '_coco_ACC_pspec' ];
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
title([ subject ' coco - EMG powerspectrum']);
if strcmp(conf.output.save,'yes')
    fig_name = [ subject '_coco_EMG_pspec' ];
    saveas(gcf, fullfile(conf.output.figdir,fig_name), 'jpg');
    saveas(gcf, fullfile(conf.output.figdir,fig_name), 'fig');
end

%% Read out the tremor frequency and the peak channel - adjust manually if necessary  
[indexACC,indexFreq] = find(fft_freq.powspctrm(:,:) == max(max(fft_freq.powspctrm((5:7),:))));

%Manually adjust peak frequency and channel for subjects that have most power in the first harmonic (tremor frequency * 2)
if strcmp(subject, '017')
    indexACC = 5;
    indexFreq = 16; %Tremor peak is at 5 Hz
end

% optional: calculate principal component
% [~,pscore,~,~,pexp] = pca([acc_x acc_y acc_z]);           % calculate principle component analysis, rows are observations, columns variables

%% Calculate the powerspectrum for all 5s segments (=NumTrials * 12) and for all frequencies and channels
power_array = [];
for segment = 1:length(segdata.trial)
    cfg=[];
    cfg.method      = 'mtmfft';
    cfg.output      = 'pow';
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
%Calculate the mean power at the peak frequency (bin) across the 12 segments of each trial from the fft_freq.powspctrm array

mean_power_per_trial_peakfreq = [];
for freq=1:8 %iterate over trials to get average power per trial
    
    total_power_per_seg_peakfreq = [];
    for time_range=1:12 %iterate over segements to get total power over peak frequency
        sum_freqbins_peakfreq = sum(power_array( ((freq-1)*12 + time_range) ,indexACC, indexFreq)); % select only the power at the peak frequency from the dominant axis
        total_power_per_seg_peakfreq = [ total_power_per_seg_peakfreq sum_freqbins_peakfreq ];
    end
   % mean_power_peakfreq = mean(total_power_per_seg_peakfreq);
    mean_power_peakfreq = mean(log(total_power_per_seg_peakfreq));
    mean_power_per_trial_peakfreq = [ mean_power_per_trial_peakfreq mean_power_peakfreq ];

end
%log_mean_power_per_trial_peakfreq = log(mean_power_per_trial_peakfreq);
log_mean_power_per_trial_peakfreq = mean_power_per_trial_peakfreq;

%Plot power values unsorted (visualize trend over time) then sort the power values and save them
[sorted_trials, indeces_trials] = sort(All_processed.trialinfo);

fig = figure('units','centimeters','outerposition',[0 0 25 20],'Color',[1 1 1]);
plot(1:8, log_mean_power_per_trial_peakfreq) %Plot the average log tremor power for the 8 trials (in sequential order) 
xlabel('Trial number');
ylabel('Average tremor power (log µV^2)');
title([ subject ' coco - Average tremor power over trials']);
%Save the plot
if strcmp(conf.output.save,'yes')
    fig_name = [ subject '_coco_average_log_tremor_power_over_trials' ];
    saveas(gcf, fullfile(conf.output.figdir,fig_name), 'jpg');
    saveas(gcf, fullfile(conf.output.figdir,fig_name), 'fig');
end

%Resort trials (power values) and save the data
log_mean_power_per_trial_peakfreq = log_mean_power_per_trial_peakfreq(1,indeces_trials); %New order trials: coco, rest

if strcmp(conf.output.save,'yes')
    power_values = log_mean_power_per_trial_peakfreq;
    save(fullfile(conf.output.datadir,[ subject '_coco_average_log_power_per_trial']),'power_values');
end

% %%%
% %threat-fix tremor power values average of the 4 trials
% coco_mean_tremor_power_values_n17_(17,1) = mean(power_values(1,1:4));
% 
% %threat-odd tremor power values average of the 4 trials
% rest_mean_tremor_power_values_n17_(17,1) = mean(power_values(1,5:8));

%% Calculate TFR for all trials and average trial per condition (i.e., create the 60s log power time course)
cfg = [];
cfg.method      = 'mtmconvol';                           % Select method (choose 'mtmconvol')
cfg.output      = 'pow';                                 % power  
cfg.taper       = 'hanning';                             % Windowing (because cut-off frequency), (Choose 'hanning' for low frequency)
cfg.foi         = 2:0.5:16;    %2:0.2:16;                          % frequency range you are interested in (usually 1:0.5:20, make sure you at least include 3-8 Hz)   
nFoi            = length(cfg.foi);
cfg.t_ftimwin   = repmat(2,1,nFoi);                      % Wavelet length (seconds; 1 wavelet per frequency). This is important also for your NaN in the hanning taper (which is 0.5*this)
cfg.toi         = All_processed.time{1}(1):0.001:All_processed.time{1}(end);% timeline the TFR (resolution in seconds) ('orig': original resolution; 'timedat': one step specified under conf.prepemg.timedat;)
cfg.pad         = 'maxperlen';                           % Padding (use 'maxperlen')
cfg.keeptrials  = 'yes';
mtmconvol_freq_trials = ft_freqanalysis(cfg, All_processed);   % The input data still contains 1s of data before and after the trial, to avoid NaNs at the beginning/end


new_Freq=round(fft_freq.freq(indexFreq) * 2) / 2;
new_indexFreq = find(mtmconvol_freq_trials.freq == new_Freq);

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
conditions = { 'coco', 'rest'};
for condition=1:length(conditions)
    
    fig = figure('units','centimeters','outerposition',[0 0 25 20],'Color',[1 1 1]);
    hold on
    
    %Plot all trials of one condition in one plot
    for trial=1:4
        total_trial_num = indeces_trials( (condition-1)*4 + trial );
        power_over_time{condition}(trial,:) = log(squeeze(mtmconvol_freq_trials.powspctrm(total_trial_num,indexACC,new_indexFreq,1001:(length(mtmconvol_freq_trials.time)-1000))))';
        plot(mtmconvol_freq_trials.time(1001:(length(mtmconvol_freq_trials.time)-1000)), power_over_time{condition}(trial,:)) %Plot 60 seconds
    end
    title([subject ' coco - ' conditions{condition} ' time courses']);
    
    %Add labels etc to the plot
    xlim([-5 63])%xlim([-8 68])%xlim([0 60000]) 
    xticks([-5 0 10 20 30 40 50 60 63]) %xticks([-8 0 10 20 30 40 50 60 68]) %xticks([-10 0 10 20 30 40 50 60 70]) %
    xticklabels({'-5','onset','10','20','30','40','50','offset','+3'})%xticklabels({'-8','onset','10','20','30','40','50','offset','+8'}) %xticklabels({'-10','onset','10','20','30','40','50','offset','+10'})
    h = xline(0,'r'); %this line is added withthe -10 + 10 range
    h1 = xline(60,'r'); %this line is added with the -10 + 10 range
    xlabel('Time (s)');
    ylabel('Tremor power (log µV^2)');
    title([subject ' coco - ' conditions{condition} ' time courses']);

    %Save the plot
    if strcmp(conf.output.save,'yes')
        fig_name = [ subject '_coco_' conditions{condition} '_time_courses' ];
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'jpg');
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'fig');
    end
    
end

%Calculate and plot average power time courses (TFR) per condition
for condition=1:length(conditions)
     %power_over_time{condition}(5,1:80000) = mean(power_over_time{condition}(1:4,1:80000),1);
    
    cfg = [];
    cfg.method      = 'mtmconvol';                           % Select method (choose 'mtmconvol')
    cfg.output      = 'pow';                                 % power
    cfg.taper       = 'hanning';                             % Windowing (because cut-off frequency), (Choose 'hanning' for low frequency)
    cfg.foi         = 2:0.5:16; %2:0.2:16;                             % frequency range you are interested in (usually 1:0.5:20, make sure you at least include 3-8 Hz)
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
    xlim([-5 63])%xlim([-8 68])%xlim([0 60000]) 
    xticks([-5 0 10 20 30 40 50 60 63]) %xticks([-8 0 10 20 30 40 50 60 68]) %xticks([-10 0 10 20 30 40 50 60 70]) %
    xticklabels({'-5','0','10','20','30','40','50','60','+3'}) %xticklabels({'-10','onset','10','20','30','40','50','offset','+10'})
    h = xline(0,'r','Trial onset'); %this line is added withthe -10 + 10 range
    h1 = xline(60,'r','Trial offset'); %this line is added with the -10 + 10 range
    xlabel('Time (s)');
    ylabel('Tremor power (log µV^2)');
    title([subject ' coco - ' conditions{condition} ' average time course']);
    
    if strcmp(conf.output.save,'yes')
        fig_name = [subject '_coco_' conditions{condition} '_average_time_course'];
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'jpg');
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'fig');
    end

end
%Save the tremor power time courses (data)

    
if strcmp(conf.output.save,'yes')
    power_time_courses = power_over_time;
    save(fullfile(conf.output.datadir,[ subject '_coco_log_power_trial_time_courses']),'power_time_courses');
end


%% Coco 2nd level - step 1.
% Create variable with group average pupil values for coco & rest
subject = '001';

cd(['/project/3024005.01/Analysis/Tremor/Results_correct_time_course/' subject '/Data']);

%coco
load([subject '_coco_log_power_trial_time_courses.mat']);
tremor_time_courses_coco_n17(17,:) = power_time_courses{1,1}(5,:);
%rest
tremor_time_courses_rest_n17(17,:) = power_time_courses{1,2}(5,:);

cd('/project/3024005.01/Analysis/Tremor/Results_correct_time_course/Group60s');
save('tremor_time_courses_coco_n17.mat','tremor_time_courses_coco_n17') 
save('tremor_time_courses_rest_n17.mat','tremor_time_courses_rest_n17') 

%Calculate the average across 17 subjects for coco
average_tremor_time_courses_coco_n17 = mean(tremor_time_courses_coco_n17);
save('average_tremor_time_courses_coco_n17.mat','average_tremor_time_courses_coco_n17') 

%Calculate the average across 17 subjects for rest
average_tremor_time_courses_rest_n17 = mean(tremor_time_courses_rest_n17);
save('average_tremor_time_courses_rest_n17.mat','average_tremor_time_courses_rest_n17') 
% 
%% time-course correlations at the trial (1,2,3,4) & condition (coco,rest) level
subject = '018';
clearvars -except subject;

%for the tremor data
cd(['/project/3024005.01/Analysis/Tremor/Results_correct_time_course','/' subject '/Data']);
load([subject '_coco_log_power_trial_time_courses.mat']); %load the coco % rest file

%for the pupil data
cd(['/project/3024005.01/Analysis/Pupil/Results','/' subject '/Data/0-60 data']);
load([subject '_coco_pupil_values_trial_time_courses.mat']); %load the coco file
load([subject '_rest_pupil_values_trial_time_courses.mat']); %load the rest file

%%%%% COCO trials %%%%%
%% Trial 1
power_time_courses{1,1}(1,5000:65000); %tremor time-course trial 1 coco
pupil_coco_time_courses(1,:); %pupil time-course trial 1 coco

%Perform the correlation
[R,P] = corr(power_time_courses{1,1}(1,5000:65000)',pupil_coco_time_courses(1,:)','Type','Pearson');

%Put the R and P value in an array which is saved later in this script
coco_trial_correlations_array(:,:) = [R P];

%% Trial 2
power_time_courses{1,1}(2,5000:65000); %tremor time-course trial 2 coco
pupil_coco_time_courses(2,:); %pupil time-course trial 2 coco

%Perform the correlation
[R,P] = corr(power_time_courses{1,1}(2,5000:65000)',pupil_coco_time_courses(2,:)','Type','Pearson');

%Put the R and P value in an array which is saved later in this script
coco_trial_correlations_array(2,:) = [R P];

%% Trial 3
power_time_courses{1,1}(3,5000:65000); %tremor time-course trial 3 coco
pupil_coco_time_courses(3,:); %pupil time-course trial 3 coco

%Perform the correlation
[R,P] = corr(power_time_courses{1,1}(3,5000:65000)',pupil_coco_time_courses(3,:)','Type','Pearson');

%Put the R and P value in an array which is saved later in this script
coco_trial_correlations_array(3,:) = [R P];

%% Trial 4
power_time_courses{1,1}(4,5000:65000); %tremor time-course trial 4 coco
pupil_coco_time_courses(4,:); %pupil time-course  trial 4 coco

%Perform the correlation
[R,P] = corr(power_time_courses{1,1}(4,5000:65000)',pupil_coco_time_courses(4,:)','Type','Pearson');

%Put the R and P value in an array which is saved later in this script
coco_trial_correlations_array(4,:) = [R P];

%% Save the trial correlations in an array
cd(['/project/3024005.01/Analysis/Tremor/Results_correct_time_course','/' subject '/Data']);
save('coco_trial_correlations_array.mat','coco_trial_correlations_array') %save the array with the R and P values 

%% Save the condition values per subject in 1 array
cd('/project/3024005.01/Analysis/Tremor/Results_correct_time_course/Group60s');
load('avg_corr_coco.mat')
avg_corr_coco(:,17) = coco_trial_correlations_array(:,1);
save('avg_corr_coco.mat','avg_corr_coco')

%% calculate the mean per per subject for coco
avg_corr_coco_n17 = mean(avg_corr_coco(:,:))';
save('avg_corr_coco_n17.mat','avg_corr_coco_n17')
%%
subject = '018';
clearvars -except subject;

%for the tremor data
cd(['/project/3024005.01/Analysis/Tremor/Results_correct_time_course','/' subject '/Data']);
load([subject '_coco_log_power_trial_time_courses.mat']); %load the coco % rest file

%for the pupil data
cd(['/project/3024005.01/Analysis/Pupil/Results','/' subject '/Data/0-60 data']);
load([subject '_coco_pupil_values_trial_time_courses.mat']); %load the coco file
load([subject '_rest_pupil_values_trial_time_courses.mat']); %load the rest file

%%%% REST trials %%%%%
%% Trial 1
power_time_courses{1,2}(1,:); %tremor time-course trial 1 rest
pupil_rest_time_courses(1,1:60000); %pupil time-course trial 1 rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,2}(1,5000:65000)',pupil_rest_time_courses(1,:)','Type','Pearson');

%Put the R and P value in an array which is saved later in this script
rest_trial_correlations_array(:,:) = [R P];

%% Trial 2
power_time_courses{1,2}(2,:); %tremor time-course trial 2 rest
pupil_rest_time_courses(2,1:60000); %pupil time-course trial 2 rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,2}(2,5000:65000)',pupil_rest_time_courses(2,:)','Type','Pearson');

%Put the R and P value in an array which is saved later in this script
rest_trial_correlations_array(2,:) = [R P];

%% Trial 3
power_time_courses{1,2}(3,:); %tremor time-course trial 3 rest
pupil_rest_time_courses(3,1:60000); %pupil time-course  trial 3 rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,2}(3,5000:65000)',pupil_rest_time_courses(3,:)','Type','Pearson');

%Put the R and P value in an array which is saved later in this script
rest_trial_correlations_array(3,:) = [R P];

%% Trial 4
power_time_courses{1,2}(4,:); %tremor time-course trial 4 rest
pupil_rest_time_courses(4,1:60000); %pupil time-course trial 4 rest

%Perform the correlation
[R,P] = corr(power_time_courses{1,2}(4,5000:65000)',pupil_rest_time_courses(4,:)','Type','Pearson');

%Put the R and P value in an array which is saved later in this script
rest_trial_correlations_array(4,:) = [R P];

%% Save the trial correlations in an array
cd(['/project/3024005.01/Analysis/Tremor/Results_correct_time_course','/' subject '/Data']);
save('rest_trial_correlations_array.mat','rest_trial_correlations_array')

%% Save the condition values per subject in 1 array
cd('/project/3024005.01/Analysis/Tremor/Results_correct_time_course/Group60s');
load('avg_corr_rest.mat')
avg_corr_rest(:,17) = rest_trial_correlations_array(:,1);
save('avg_corr_rest.mat','avg_corr_rest')

%% calculate the mean per per subject for rest
avg_corr_rest_n17 = mean(avg_corr_rest(:,:))';
save('avg_corr_rest_n17.mat','avg_corr_rest_n17')






%% for barplotting with avnBarplotscatter.m

avg_corr_coco_cell = num2cell(avg_corr_coco,[1 2]);
avg_corr_coco_cell{1,2} = (avg_corr_rest);

save('avg_corr_coco_cell.mat', 'avg_corr_coco_cell')
end
% 
% 
% 
% 
% 







%%%%%%%%%%%%%%%%%%% OLD scripts for correlations at group level %%%%%%%%%%%%%%%%%%%%%%%%%
%% Correlating coco time-courses on a group level for tremor & pupil

%%% COCO & REST files %%%
%cd('/project/3024005.01/Analysis/Tremor/NewResults/Group60s/0-60 time-courses');
%load('average_power_time_courses_coco_n17.mat');
%load('average_power_time_courses_rest_n17.mat');

%cd('/project/3024005.01/Analysis/Pupil/Results/group-level (60 sec)/coco/0-60 data'); %load directory with coco group level time course values
%load('average_pupil_time_courses_coco_n17.mat');
%load('average_pupil_time_courses_rest_n17.mat');


%% COCO average %%%%
%average_power_time_courses_coco_n17(1,1:60000); %average tremor time-course coco condition
%average_pupil_time_courses_coco_n17(1,1:60000); %average pupil time-course  coco condition

% Down-sampling of the tremor time-course to 1 sec segments
%power_time_courses_group_coco = average_power_time_courses_coco_n17(1,:);
%n = 1000; % average every n values
%power_time_courses_60s_group_coco = arrayfun(@(i) mean(power_time_courses_group_coco(i:i+n-1)),1:n:length(power_time_courses_group_coco)-n+1)'; % the averaged vector

% Down-sampling of the pupil time-course to 1 sec segments
%pupil_time_courses_group_coco = average_pupil_time_courses_coco_n17(1,1:60000);
%n = 1000; % average every n values
%pupil_time_courses_60s_group_coco = arrayfun(@(i) mean(pupil_time_courses_group_coco(i:i+n-1)),1:n:length(pupil_time_courses_group_coco)-n+1)'; % the averaged vector
%
%[RHO,PVAL] = corr(power_time_courses_60s_group_coco,pupil_time_courses_60s_group_coco,'Type','Spearman');

%Put the R and P value in an array which is saved later in this script
%coco_group_correlations_array(:,:) = [RHO PVAL];

%cd('/project/3024005.01/Analysis/Tremor/NewResults/Group60s/0-60 time-courses');
%save('coco_group_correlations_array.mat','coco_group_correlations_array')

% Plot the correlation for coco at the group level
%coco_group_correlations_figure = figure('units','centimeters','outerposition',[0 0 25 20],'Color',[1 1 1]);
%s1 = plot(power_time_courses_60s_group_coco(:,1),pupil_time_courses_60s_group_coco(:,1), '+', 'MarkerFaceColor', 'k');
%%    ylabel('Pupil diamter','FontSize', 12);
%%    xlabel ('Tremor power','FontSize', 12);
%    title('coco correlation plot (N = 17)');
%    set(s1, 'MarkerSize', 8, 'LineWidth', 2);
%Add the correlation coefficient
%    disp(RHO(1,1));
%    str=['RHO= ',num2str(RHO(1,1))]
%    RHOVAL = text(min(get(gca, 'xlim')+1), max(get(gca, 'ylim')-15), str);
%    set(RHOVAL, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left')
%Add the p-value
%    disp(PVAL(1,1));
%    str=['P= ',num2str(PVAL(1,1))]
%%    SIGVAL = text(min(get(gca, 'xlim')+1), max(get(gca, 'ylim')-30), str);
%    set(SIGVAL, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left')    
%Add a regression line
%%    hold on
%    l = lsline ;
%    set(l,'LineWidth', 2)
%Save figures in .fig and .jpg format
%    saveas(coco_group_correlations_figure,'coco_group_correlations_figure.jpg')
%    savefig('coco_group_correlations_figure.fig')
    


