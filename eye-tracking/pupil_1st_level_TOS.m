%% Pupil data - TOS

subject = '001'; %%%%% specify subject number

clc; close all; clearvars -except subject; %%%%% make sure to have a clean workspace

%% create correct path, add functions and scripts
restoredefaultpath;
addpath(genpath('/project/3024005.01/Analysis/SupportingScripts/'));
addpath('/project/3024005.01/Analysis/Pupil');

%% Define input directory and files
conf.input.rawdir = '/project/3024005.01/Analysis/Pupil/';
conf.input.matdir = fullfile(pf_findfile(conf.input.rawdir,['/' subject '/'],'fullfile'));
conf.input.matfile = pf_findfile(conf.input.matdir,['/' subject '/&/TOS/&/.mat/'],'fullfile');

%% Define where saved data and figures are stored 
conf.output.dir = [ '/project/3024005.01/Analysis/Pupil/Results/' subject];
conf.output.figdir = [ '/project/3024005.01/Analysis/Pupil/Results/' subject '/Figures' ];
conf.output.datadir = [ '/project/3024005.01/Analysis/Pupil/Results/' subject '/Data' ];
conf.output.save = 'yes';
if ~exist(conf.output.dir)
    mkdir(conf.output.dir);
end

%% Load the the file
filename = conf.input.matfile;
load (filename);

%% Define what you want to analyse from variables pupdat and interpoldat 
conf.ana.var     =  {
%                      'Raw X';
%                      'Raw Y';
                       'Pupil';}; % pupdat.names you want to analyse (e.g. Pupil)


conf.ana.type    =  'fullinterpoldat'; % field saved in interpoldat that you want to plot
conf.ana.check   =  {
%                       'all';
                        'usable';
%                       'non-usable';
                        }; % Indicate if you want to use data that has been marked usable, non-usable, indeterminate or all            

%% Select the data to analyse 
nVar = length(conf.ana.var); %the amount of variables analyzed
CurVar = conf.ana.var; %the name of the current variable that is analyzed
fullinterpoldat  = interpoldat.(conf.ana.type); %indicates the fullinterpoldat structure that will be used for plotting    
%% Get the onset
onset  = pupdat.trialonsets{1, 3}; %onset of trials for pupil as defined in pupdat.trialonset
offset = pupdat.trialonsets{1, 3} + 60*pupdat.samplerate; % offset of trials for pupil 
nOns    = length(onset); %get the number of onsets

%% %% TOS %%


%defining the individual TOS trials according to the onsets and offsets
trial1TOS = fullinterpoldat{1, 3}(onset(1,1)-6000:offset(1,1)+3000); %(onset(1,1):offset(1,1)); %%fill in onset(1,1)-9792 for Subject, because its 10 max sec before 1st trial onset is 9.792 sec -> see onset 17x1 double. Also for subject 11 %fill in onset(1,1)-7972 for Subject 11, because its 10 max sec before 1st trial onset is 7.972 sec -> see onset 17x1 double. Also for subject 15 %fill in onset(1,1)-9360 for Subject 15, because its 10 max sec before 1st trial onset is 9.361 sec -> see onset 17x1 double. Also for subject 16 %fill in onset(1,1)-7956 for Subject 16, because its 10 max sec before 1st trial onset is 7.956 sec -> see onset 17x1 double.
trial2TOS = fullinterpoldat{1, 3}(onset(2,1)-6000:offset(2,1)+3000);
trial3TOS = fullinterpoldat{1, 3}(onset(3,1)-6000:offset(3,1)+3000);
trial4TOS = fullinterpoldat{1, 3}(onset(4,1)-6000:offset(4,1)+3000);
trial5TOS = fullinterpoldat{1, 3}(onset(5,1)-6000:offset(5,1)+3000);
trial6TOS = fullinterpoldat{1, 3}(onset(6,1)-6000:offset(6,1)+3000);
trial7TOS = fullinterpoldat{1, 3}(onset(7,1)-6000:offset(7,1)+3000);
trial8TOS = fullinterpoldat{1, 3}(onset(8,1)-6000:offset(8,1)+3000);
trial9TOS = fullinterpoldat{1, 3}(onset(9,1)-6000:offset(9,1)+3000);
trial10TOS = fullinterpoldat{1, 3}(onset(10,1)-6000:offset(10,1)+3000);
trial11TOS = fullinterpoldat{1, 3}(onset(11,1)-6000:offset(11,1)+3000);
trial12TOS = fullinterpoldat{1, 3}(onset(12,1)-6000:offset(12,1)+3000);
trial13TOS = fullinterpoldat{1, 3}(onset(13,1)-6000:offset(13,1)+3000);
trial14TOS = fullinterpoldat{1, 3}(onset(14,1)-6000:offset(14,1)+3000);
trial15TOS = fullinterpoldat{1, 3}(onset(15,1)-6000:offset(15,1)+3000);
trial16TOS = fullinterpoldat{1, 3}(onset(16,1)-6000:offset(16,1)+3000);
trial17TOS = fullinterpoldat{1, 3}(onset(17,1)-6000:offset(17,1)+3000);
alltrialsTOS = fullinterpoldat{1, 3}(:); % all TOS trials together
alldataTOS = fullinterpoldat{1, 3}; %all data from begin (0) to end

%for subject 10 add 208 datapoints (10000-9792(onset)) (start number at trial onset is 5370) in trial1TOS
%extradata_trial1TOS(1:208,1) = 5370;%first value at onset
%trial1TOS(209:80001,1) = trial1TOS;
%trial1TOS(1:208,1) = extradata_trial1TOS;

%for subject 11 add 2025 datapoints (10000-7975 (onset))(start number at trial onset is 4826) in trial1TOS 
%extradata_trial1TOS(1:2025,1) = 4826;%first value at onset
%trial1TOS(2026:80001,1) = trial1TOS;
%trial1TOS(1:2025,1) = extradata_trial1TOS;

%for subject 15 add 639 datapoints (10000-9361 (onset))(start number at trial onset is 4083) in trial1TOS 
%extradata_trial1TOS(1:640,1) = 4083;%first value at onset
%trial1TOS(641:80001,1) = trial1TOS;
%trial1TOS(1:639,1) = extradata_trial1TOS;

%for subject 16 add 2044 datapoints (10000-7956 (onset)) (start number at trial onset is 4288) in trial1TOS 
%extradata_trial1TOS(1:2044,1) = 4288; %first value at onset
%trial1TOS(2045:80001,1) = trial1TOS;
%trial1TOS(1:2044,1) = extradata_trial1TOS;

%for subject 17 add 220 datapoints (10000-9780 (onset)) (start number at trial onset is 4288) in trial1TOS 
%extradata_trial1TOS(1:220,1) = 3787; %first value at onset
%trial1TOS(221:80001,1) = trial1TOS;
%trial1TOS(1:220,1) = extradata_trial1TOS;
%% For all trials calculate mean and SD
%%%%%% TOS %%%%%%

%Trial 1 TOS 
mean_trial1TOS = mean(trial1TOS); %
stdev_trial1TOS = std(trial1TOS); %

%Trial 2 TOS 
mean_trial2TOS = mean(trial2TOS); %
stdev_trial2TOS = std(trial2TOS); %

%Trial 3 TOS 
mean_trial3TOS = mean(trial3TOS); %
stdev_trial3TOS = std(trial3TOS); %

%Trial 4 TOS 
mean_trial4TOS = mean(trial4TOS); %
stdev_trial4TOS = std(trial4TOS); %

%Trial 5 TOS 
mean_trial5TOS = mean(trial5TOS); %
stdev_trial5TOS = std(trial5TOS); %

%Trial 6 TOS 
mean_trial6TOS = mean(trial6TOS); %
stdev_trial6TOS = std(trial6TOS); %

%Trial 7 TOS 
mean_trial7TOS = mean(trial7TOS); %
stdev_trial7TOS = std(trial7TOS); %

%Trial 8 TOS 
mean_trial8TOS = mean(trial8TOS); %
stdev_trial8TOS = std(trial8TOS); %

%Trial 9 TOS 
mean_trial9TOS = mean(trial9TOS); %
stdev_trial9TOS = std(trial9TOS); %

%Trial 10 TOS 
mean_trial10TOS = mean(trial10TOS); %
stdev_trial10TOS = std(trial10TOS); %

%Trial 11 TOS 
mean_trial11TOS = mean(trial11TOS); %
stdev_trial11TOS = std(trial11TOS); %

%Trial 12 TOS 
mean_trial12TOS = mean(trial12TOS); %
stdev_trial12TOS = std(trial12TOS); %

%Trial 13 TOS 
mean_trial13TOS = mean(trial13TOS); %
stdev_trial13TOS = std(trial13TOS); %

%Trial 14 TOS 
mean_trial14TOS = mean(trial14TOS); %
stdev_trial14TOS = std(trial14TOS); %

%Trial 15 TOS 
mean_trial15TOS = mean(trial15TOS); %
stdev_trial15TOS = std(trial15TOS); %

%Trial 16 TOS 
mean_trial16TOS = mean(trial16TOS); %
stdev_trial16TOS = std(trial16TOS); %

%Trial 17 TOS 
mean_trial17TOS = mean(trial17TOS); %
stdev_trial17TOS = std(trial17TOS); %

%Put all the values together to calculate the mean of all trials in the original sequence order
TOS_mean_vector = [mean_trial1TOS;mean_trial2TOS;mean_trial3TOS;mean_trial4TOS;mean_trial5TOS;mean_trial6TOS;mean_trial7TOS;mean_trial8TOS;mean_trial9TOS;mean_trial10TOS;mean_trial11TOS;mean_trial12TOS;mean_trial13TOS;mean_trial14TOS;mean_trial15TOS;mean_trial16TOS;mean_trial17TOS];
TOS_sd_vector = [stdev_trial1TOS;stdev_trial2TOS;stdev_trial3TOS;stdev_trial4TOS;stdev_trial5TOS;stdev_trial6TOS;stdev_trial7TOS;stdev_trial8TOS;stdev_trial9TOS;stdev_trial10TOS;stdev_trial11TOS;stdev_trial12TOS;stdev_trial13TOS;stdev_trial14TOS;stdev_trial15TOS;stdev_trial16TOS;stdev_trial17TOS];


%% TOS %% Time-course plotting of pupil data - Step 1
%Plot pupil values unsorted (visualize trend over time) then sort the pupil values and save them
[sorted_trials, indeces_trials] = sort(pupdat.trialcodes{1,3});
sorted_trials = [sorted_trials(1:5); sorted_trials(7:17); sorted_trials(6)]; %Move 2nd Threat_odd block (=shock block) to the end
indeces_trials = [indeces_trials(1:5); indeces_trials(7:17); indeces_trials(6)]; %Use this array also later for sorting time courses trials


fig = figure('units','centimeters','outerposition',[0 0 25 20],'Color',[1 1 1]);
plot(1:17, TOS_mean_vector) %Plot the average pupil values for the 8 trials (in sequential order) 
xlabel('Trial number');
ylabel('Average pupil dilation');
title([ subject ' TOS - Average pupil dilation over trials']);
h = xline(indeces_trials(end),'--r','Shock trial');
%Save the plot
if strcmp(conf.output.save,'yes')
    fig_name = [ subject '_TOS_average_pupil_dilation_over_trials' ]; %  fig_name = [ subject '_TOS_average_pupil_dilation_over_trials' ];
    saveas(gcf, fullfile(conf.output.figdir,fig_name), 'jpg');
    saveas(gcf, fullfile(conf.output.figdir,fig_name), 'fig');
end

%Resort trials (power values) and save the data
TOS_mean_vector = TOS_mean_vector'; %First reshape the data structure
TOS_mean_vector = TOS_mean_vector(1,indeces_trials); %New order trials: thr_fix, thr_odd, safe_fix, safe_odd, shock

if strcmp(conf.output.save,'yes')
    pupil_values= TOS_mean_vector; %pupil_values = TOS_mean_vector;
    save(fullfile(conf.output.datadir,[ subject '_TOS_average_pupil_dilation_per_trial']),'pupil_values'); % save(fullfile(conf.output.datadir,[ subject '_TOS_average_pupil_dilation_per_trial']),'pupil_values');
end

%% TOS %% Time-course plotting of pupil data - Step 2
% Put all the trial values together and resort the trials
TOS_alltrials = [trial1TOS,trial2TOS,trial3TOS,trial4TOS,trial5TOS,trial6TOS,trial7TOS,trial8TOS,trial9TOS,trial10TOS,trial11TOS,trial12TOS,trial13TOS,trial14TOS,trial15TOS,trial16TOS,trial17TOS];

%Resort trials (power values) and save the data
TOS_alltrials = TOS_alltrials(:,indeces_trials); %New order trials: thr_fix, thr_odd, safe_fix, safe_odd, shock
TOS_alltrials = TOS_alltrials'; %Reshape the data structure

if strcmp(conf.output.save,'yes')
    pupil_values_TOS = TOS_alltrials; %pupil_values_TOS = TOS_alltrials;
    save(fullfile(conf.output.datadir,[ subject '_TOS_pupil_values_per_trial']),'pupil_values_TOS'); %save(fullfile(conf.output.datadir,[ subject '_TOS_pupil_values_per_trial']),'pupil_values_TOS');
end

%Calculate and plot the (absolute) pupil values time all 8 trials (trials of the same condition plotted in 1 figure)
TOS_thr_fix_average = mean(TOS_alltrials(1:4,:));
TOS_thr_odd_average = mean(TOS_alltrials(5:8,:));
TOS_safe_fix_average = mean(TOS_alltrials(9:12,:));
TOS_safe_odd_average = mean(TOS_alltrials(13:16,:));
TOS_shock_average = (TOS_alltrials(17,:));

TOS_conditions = {TOS_thr_fix_average,TOS_thr_odd_average,TOS_safe_fix_average,TOS_safe_odd_average,TOS_shock_average };

conditions = { 'threat-rest'};
for condition=1:length(conditions)
    
    fig = figure('units','centimeters','outerposition',[0 0 25 20],'Color',[1 1 1]);
    hold on
    plot(TOS_conditions{1,1}(:,:)); %plot(TOS_conditions{1,1}(1,1:60001));
   h = xline(6000,'r','Trial onset'); %this line is added with the -10 + 10 range
   h1 = xline(66000,'r','Trial offset'); %this line is added with the -10 + 10 range
    title([subject ' threat-rest ' ' average time course']);
end 

%Add labels etc to the plot
    xlim([0 69000]) %xlim([0 60000])
    xticks([0 6000 16000 26000 36000 46000 56000 66000 69000]) %xticks([0 10000 20000 30000 40000 50000 60000])
    xticklabels({'-6','0','10','20','30','40','50','60','+3'}) %xticklabels(}'0','10','20','30','40','50','60'})
    xlabel('Time (s)');
    ylabel('Pupil dilation');

 %Save the plot
    if strcmp(conf.output.save,'yes')
        fig_name = [ subject '_threat_rest_' conditions{condition} '_time_courses' ]; %  fig_name = [ subject '_threat_rest_' conditions{condition} '_time_courses' ];
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'jpg');
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'fig');
    end
        
%Save the pupil values time courses threat-rest (data)
if strcmp(conf.output.save,'yes')
    pupil_threat_rest_time_courses = TOS_alltrials(1:4,:); %pupil_threat_rest_time_courses = TOS_alltrials(1:4,:);
    save(fullfile(conf.output.datadir,[ subject '_threat_rest_pupil_values_trial_time_courses']),'pupil_threat_rest_time_courses'); %save(fullfile(conf.output.datadir,[ subject '_threat_rest_pupil_values_trial_time_courses']),'pupil_threat_rest_time_courses');
end

%Save the average pupil values time course threat-rest (data)
if strcmp(conf.output.save,'yes')
    average_pupil_threat_rest_time_courses = TOS_conditions{1,1}(:,:); %pupil_threat_rest_time_courses = TOS_alltrials(1:4,:);
    save(fullfile(conf.output.datadir,[ subject '_average_threat_rest_pupil_values_trial_time_courses']),'average_pupil_threat_rest_time_courses'); %save(fullfile(conf.output.datadir,[ subject '_threat_rest_pupil_values_trial_time_courses']),'pupil_threat_rest_time_courses');
end



conditions = { 'threat-odd'};
for condition=1:length(conditions)
    
    fig = figure('units','centimeters','outerposition',[0 0 25 20],'Color',[1 1 1]);
    hold on
    plot(TOS_conditions{1,2}(:,:)); %plot(TOS_conditions{1,2}(1,1:60001));
    h = xline(6000,'r','Trial onset'); %this line is added with the -10 + 10 range
    h1 = xline(66000,'r','Trial offset'); %this line is added with the -10 + 10 range
    title([subject ' threat-odd ' ' average time course']);
end 

%Add labels etc to the plot
    xlim([0 69000]) %xlim([0 60000])
    xticks([0 6000 16000 26000 36000 46000 56000 66000 69000]) %xticks([0 10000 20000 30000 40000 50000 60000])
    xticklabels({'-6','0','10','20','30','40','50','60','+3'}) %xticklabels(}'0','10','20','30','40','50','60'})
    xlabel('Time (s)');
    ylabel('Pupil dilation');
    
 %Save the plot
    if strcmp(conf.output.save,'yes')
        fig_name = [ subject '_threat_odd_' conditions{condition} '_time_courses' ]; %  fig_name = [ subject '_threat_odd_' conditions{condition} '_time_courses' ];
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'jpg');
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'fig');
    end
        
%Save the pupil values time courses threat-odd (data)
if strcmp(conf.output.save,'yes')
    pupil_threat_odd_time_courses= TOS_alltrials(5:8,:); %pupil_threat_odd_time_courses = TOS_alltrials(5:8,:);
    save(fullfile(conf.output.datadir,[ subject '_threat_odd_pupil_values_trial_time_courses']),'pupil_threat_odd_time_courses');%save(fullfile(conf.output.datadir,[ subject '_threat_odd_pupil_values_trial_time_courses']),'pupil_threat_odd_time_courses');
end

%Save the average pupil values time course threat-rest (data)
if strcmp(conf.output.save,'yes')
    average_pupil_threat_odd_time_courses = TOS_conditions{1,2}(:,:); %pupil_threat_rest_time_courses = TOS_alltrials(1:4,:);
    save(fullfile(conf.output.datadir,[ subject '_average_threat_odd_pupil_values_trial_time_courses']),'average_pupil_threat_odd_time_courses'); %save(fullfile(conf.output.datadir,[ subject '_threat_rest_pupil_values_trial_time_courses']),'pupil_threat_rest_time_courses');
end

conditions = { 'safe-rest'};
for condition=1:length(conditions)
    
    fig = figure('units','centimeters','outerposition',[0 0 25 20],'Color',[1 1 1]);
    hold on
    plot(TOS_conditions{1,3}(:,:)); %plot(TOS_conditions{1,3}(1,1:60001));
    h = xline(6000,'r','Trial onset'); %this line is added with the -10 + 10 range
    h1 = xline(66000,'r','Trial offset'); %this line is added with the -10 + 10 range
    title([subject ' safe-rest ' ' average time course']);
end 

%Add labels etc to the plot
    xlim([0 69000]) %xlim([0 60000])
    xticks([0 6000 16000 26000 36000 46000 56000 66000 69000]) %xticks([0 10000 20000 30000 40000 50000 60000])
    xticklabels({'-6','0','10','20','30','40','50','60','+3'}) %xticklabels(}'0','10','20','30','40','50','60'})
    xlabel('Time (s)');
    ylabel('Pupil dilation');

 %Save the plot
    if strcmp(conf.output.save,'yes')
        fig_name = [ subject '_safe_rest_' conditions{condition} '_time_courses' ]; %  fig_name = [ subject '_safe_rest_' conditions{condition} '_time_courses' ];
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'jpg');
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'fig');
    end
        
%Save the pupil values time courses threat-odd (data)
if strcmp(conf.output.save,'yes')
    pupil_safe_rest_time_courses = TOS_alltrials(9:12,:); %pupil_safe_rest_time_courses = TOS_alltrials(9:12,:);
    save(fullfile(conf.output.datadir,[ subject '_safe_rest_pupil_values_trial_time_courses']),'pupil_safe_rest_time_courses'); %    save(fullfile(conf.output.datadir,[ subject '_safe_rest_pupil_values_trial_time_courses']),'pupil_safe_rest_time_courses');
end

%Save the average pupil values time course threat-rest (data)
if strcmp(conf.output.save,'yes')
    average_pupil_safe_rest_time_courses = TOS_conditions{1,3}(:,:); %pupil_threat_rest_time_courses = TOS_alltrials(1:4,:);
    save(fullfile(conf.output.datadir,[ subject '_average_safe_rest_pupil_values_trial_time_courses']),'average_pupil_safe_rest_time_courses'); %save(fullfile(conf.output.datadir,[ subject '_threat_rest_pupil_values_trial_time_courses']),'pupil_threat_rest_time_courses');
end

conditions = { 'safe-odd'};
for condition=1:length(conditions)
    
    fig = figure('units','centimeters','outerposition',[0 0 25 20],'Color',[1 1 1]);
    hold on
    plot(TOS_conditions{1,4}(:,:)); %plot(TOS_conditions{1,4}(1,1:60001));
    h = xline(6000,'r','Trial onset'); %this line is added with the -10 + 10 range
    h1 = xline(66000,'r','Trial offset'); %this line is added with the -10 + 10 range
    title([subject ' safe-odd ' ' average time course']);
end 

%Add labels etc to the plot
    xlim([0 69000]) %xlim([0 60000])
    xticks([0 6000 16000 26000 36000 46000 56000 66000 69000]) %xticks([0 10000 20000 30000 40000 50000 60000])
    xticklabels({'-6','0','10','20','30','40','50','60','+3'}) %xticklabels(}'0','10','20','30','40','50','60'})
    xlabel('Time (s)');
    ylabel('Pupil dilation');
    
 %Save the plot
    if strcmp(conf.output.save,'yes')
        fig_name = [ subject '_safe_odd_' conditions{condition} '_time_courses' ]; %  fig_name = [ subject '_safe_odd_' conditions{condition} '_time_courses' ];
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'jpg');
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'fig');
    end
        
%Save the pupil values time courses threat-odd (data)
if strcmp(conf.output.save,'yes')
    pupil_safe_odd_time_courses = TOS_alltrials(13:16,:); %pupil_safe_odd_time_courses = TOS_alltrials(13:16,:);
    save(fullfile(conf.output.datadir,[ subject '_safe_odd_pupil_values_trial_time_courses']),'pupil_safe_odd_time_courses'); %save(fullfile(conf.output.datadir,[ subject '_safe_odd_pupil_values_trial_time_courses']),'pupil_safe_odd_time_courses');
end

%Save the average pupil values time course threat-rest (data)
if strcmp(conf.output.save,'yes')
    average_pupil_safe_odd_time_courses = TOS_conditions{1,4}(:,:); %pupil_threat_rest_time_courses = TOS_alltrials(1:4,:);
    save(fullfile(conf.output.datadir,[ subject '_average_safe_odd_pupil_values_trial_time_courses']),'average_pupil_safe_odd_time_courses'); %save(fullfile(conf.output.datadir,[ subject '_threat_rest_pupil_values_trial_time_courses']),'pupil_threat_rest_time_courses');
end

conditions = { 'shock'};
for condition=1:length(conditions)
    
    fig = figure('units','centimeters','outerposition',[0 0 25 20],'Color',[1 1 1]);
    hold on
    plot(TOS_conditions{1,5}(:,:)); %plot(TOS_conditions{1,5}(1,1:60001));
    h = xline(6000,'r','Trial onset'); %this line is added with the -10 + 10 range
    h1 = xline(66000,'r','Trial offset'); %this line is added with the -10 + 10 range
    title([subject ' shock ' ' average time course']);
end 

%Add labels etc to the plot
    xlim([0 69000]) %xlim([0 60000])
    xticks([0 6000 16000 26000 36000 46000 56000 66000 69000]) %xticks([0 10000 20000 30000 40000 50000 60000])
    xticklabels({'-6','0','10','20','30','40','50','60','+3'}) %xticklabels(}'0','10','20','30','40','50','60'})
    xlabel('Time (s)');
    ylabel('Pupil dilation');
    
 %Save the plot
    if strcmp(conf.output.save,'yes')
        fig_name = [ subject '_shock_' conditions{condition} '_time_courses' ]; %  fig_name = [ subject '_shock_' conditions{condition} '_time_courses' ];
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'jpg');
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'fig');
    end
        
%Save the pupil values time courses threat-odd (data)
if strcmp(conf.output.save,'yes')
    pupil_shock_time_courses = TOS_alltrials(17,:); %pupil_shock_time_courses = TOS_alltrials(17,:);
    save(fullfile(conf.output.datadir,[ subject '_shock_pupil_values_trial_time_courses']),'pupil_shock_time_courses');  %save(fullfile(conf.output.datadir,[ subject '_shock_pupil_values_trial_time_courses']),'pupil_shock_time_courses');
end

%Save the average pupil values time course threat-rest (data)
if strcmp(conf.output.save,'yes')
    average_pupil_shock_time_courses = TOS_conditions{1,5}(:,:); %pupil_threat_rest_time_courses = TOS_alltrials(1:4,:);
    save(fullfile(conf.output.datadir,[ subject '_average_shock_pupil_values_trial_time_courses']),'average_pupil_shock_time_courses'); %save(fullfile(conf.output.datadir,[ subject '_threat_rest_pupil_values_trial_time_courses']),'pupil_threat_rest_time_courses');
end

%% Create variable with group average pupil values for coco & rest

subject = '018';

cd(['/project/3024005.01/Analysis/Pupil/Results/' subject '/Data']);

%threat-rest
load([subject '_threat_rest_pupil_values_trial_time_courses.mat']);
pupil_time_courses_threat_rest_n17(17,:) = mean(pupil_threat_rest_time_courses(1:4,:));

%threat-odd
load([subject '_threat_odd_pupil_values_trial_time_courses.mat']);
pupil_time_courses_threat_odd_n17(17,:) = mean(pupil_threat_odd_time_courses(1:4,:));

%safe-rest
load([subject '_safe_rest_pupil_values_trial_time_courses.mat']);
pupil_time_courses_safe_rest_n17(17,:) = mean(pupil_safe_rest_time_courses(1:4,:));

%safe-odd
load([subject '_safe_odd_pupil_values_trial_time_courses.mat']);
pupil_time_courses_safe_odd_n17(17,:) = mean(pupil_safe_odd_time_courses(1:4,:));

%shock
load([subject '_shock_pupil_values_trial_time_courses.mat']);
pupil_time_courses_shock_n17(17,:) = pupil_shock_time_courses(1,:);

cd('/project/3024005.01/Analysis/Pupil/Results/group-level (60 sec)/TOS');
save('pupil_time_courses_threat_rest_n17.mat','pupil_time_courses_threat_rest_n17') 
save('pupil_time_courses_threat_odd_n17.mat','pupil_time_courses_threat_odd_n17') 
save('pupil_time_courses_safe_rest_n17.mat','pupil_time_courses_safe_rest_n17') 
save('pupil_time_courses_safe_odd_n17.mat','pupil_time_courses_safe_odd_n17') 
save('pupil_time_courses_shock_n17.mat','pupil_time_courses_shock_n17') 

%Calculate the average across 17 subjects for threat-rest
average_pupil_time_courses_threat_rest_n17 = mean(pupil_time_courses_threat_rest_n17);
save('average_pupil_time_courses_threat_rest_n17.mat','average_pupil_time_courses_threat_rest_n17') 

%Calculate the average across 17 subjects for threat-rest
average_pupil_time_courses_threat_odd_n17 = mean(pupil_time_courses_threat_odd_n17);
save('average_pupil_time_courses_threat_odd_n17.mat','average_pupil_time_courses_threat_odd_n17') 

%Calculate the average across 17 subjects for threat-rest
average_pupil_time_courses_safe_rest_n17 = mean(pupil_time_courses_safe_rest_n17);
save('average_pupil_time_courses_safe_rest_n17.mat','average_pupil_time_courses_safe_rest_n17') 

%Calculate the average across 17 subjects for threat-rest
average_pupil_time_courses_safe_odd_n17 = mean(pupil_time_courses_safe_odd_n17);
save('average_pupil_time_courses_safe_odd_n17.mat','average_pupil_time_courses_safe_odd_n17') 

%Calculate the average across 17 subjects for threat-rest
average_pupil_time_courses_shock_n17 = mean(pupil_time_courses_shock_n17);
save('average_pupil_time_courses_shock_n17.mat','average_pupil_time_courses_shock_n17') 
