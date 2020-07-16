%% Pupil data coco

subject = '018'; %%%%% specify subject number

clc; close all; clearvars -except subject; %%%%% make sure to have a clean workspace

%% create correct path, add functions and scripts
restoredefaultpath;
addpath(genpath('/project/3024005.01/Analysis/SupportingScripts/'));
addpath('/project/3024005.01/Analysis/Pupil');

%% Define input directory and files
conf.input.rawdir = '/project/3024005.01/Analysis/Pupil/';
conf.input.matdir = fullfile(pf_findfile(conf.input.rawdir,['/' subject '/'],'fullfile'));
conf.input.matfile = pf_findfile(conf.input.matdir,['/' subject '/&/coco/&/.mat/'],'fullfile');


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

%% %%%% Coco %%%%%% ADJUSTED TO 10 sec before trial onset and 10 sec after trial offset
%defining the individual coco & rest trials according to the onsets and offsets
trial1coco = fullinterpoldat{1, 3}(onset(1,1)-6000:offset(1,1)+3000); %trial1coco = fullinterpoldat{1, 3}(onset(1,1):offset(1,1));
trial1rest = fullinterpoldat{1, 3}(onset(2,1)-6000:offset(2,1)+3000);
trial2coco = fullinterpoldat{1, 3}(onset(3,1)-6000:offset(3,1)+3000);
trial2rest = fullinterpoldat{1, 3}(onset(4,1)-6000:offset(4,1)+3000);
trial3coco = fullinterpoldat{1, 3}(onset(5,1)-6000:offset(5,1)+3000);
trial3rest = fullinterpoldat{1, 3}(onset(6,1)-6000:offset(6,1)+3000);
trial4coco = fullinterpoldat{1, 3}(onset(7,1)-6000:offset(7,1)+3000);
trial4rest = fullinterpoldat{1, 3}(onset(8,1)-6000:offset(8,1)+3000);
alltrialsCOCO = fullinterpoldat{1, 3}(:); %all coco & rest trials together
alldataCOCO = fullinterpoldat{1, 3};%all data from begin to end

%%%%%% Coco %%%%%%

%Trial 1 coco
mean_trial1coco = mean(trial1coco); %
stdev_trial1coco = std(trial1coco); %

%Trial 1 rest
mean_trial1rest = mean(trial1rest); %
stdev_trial1rest = std(trial1rest); %

%Trial 2 coco
mean_trial2coco = mean(trial2coco); %
stdev_trial2coco = std(trial2coco); %

%Trial 2 rest
mean_trial2rest = mean(trial2rest); %
stdev_trial2rest = std(trial2rest); %

%Trial 3 coco
mean_trial3coco = mean(trial3coco); %
stdev_trial3coco = std(trial3coco); %

%Trial 3 rest
mean_trial3rest = mean(trial3rest); %
stdev_trial3rest = std(trial3rest); %

%Trial 4 coco
mean_trial4coco = mean(trial4coco); %
stdev_trial4coco = std(trial4coco); %

%Trial 4 rest
mean_trial4rest = mean(trial4rest); %
stdev_trial4rest = std(trial4rest); %

%Put all the mean + SD values together
coco_mean_vector = [mean_trial1coco;mean_trial1rest;mean_trial2coco;mean_trial2rest;mean_trial3coco;mean_trial3rest;mean_trial4coco;mean_trial4rest];
coco_sd_vector = [stdev_trial1coco;stdev_trial1rest;stdev_trial2coco;stdev_trial2rest;stdev_trial3coco;stdev_trial3rest;stdev_trial4coco;stdev_trial4rest];

%% COCO %% Time-course plotting of pupil data 
%Plot pupil values unsorted (visualize trend over time) then sort the pupil values and save them
[sorted_trials, indeces_trials] = sort(pupdat.trialcodes{1,3});

fig = figure('units','centimeters','outerposition',[0 0 25 20],'Color',[1 1 1]);
plot(1:8, coco_mean_vector) %Plot the average pupil values for the 8 trials (in sequential order) 
xlabel('Trial number');
ylabel('Average pupil dilation');
title([ subject ' coco - Average pupil dilation over trials']);

%Save the plot
if strcmp(conf.output.save,'yes')
    fig_name = [ subject '_coco_average_pupil_dilation_over_trials' ]; %fig_name = [ subject '_coco_average_pupil_dilation_over_trials' ]; 
    saveas(gcf, fullfile(conf.output.figdir,fig_name), 'jpg');
    saveas(gcf, fullfile(conf.output.figdir,fig_name), 'fig');
end

%Resort trials (power values) and save the data
coco_mean_vector = coco_mean_vector'; %First reshape the data structure
coco_mean_vector = coco_mean_vector(1,indeces_trials); %New order trials: coco, rest

if strcmp(conf.output.save,'yes')
    pupil_values = coco_mean_vector;
    save(fullfile(conf.output.datadir,[ subject '_coco_average_pupil_dilation_per_trial']),'pupil_values'); % save(fullfile(conf.output.datadir,[ subject '_coco_average_pupil_dilation_per_trial']),'pupil_values');
end


%Calculate and plot the (absolute) pupil values time all 8 trials (trials of the same condition plotted in 1 figure)

trial1coco = trial1coco';
trial2coco = trial2coco';
trial3coco = trial3coco';
trial4coco = trial4coco';
trialscoco = [trial1coco;trial2coco;trial3coco;trial4coco];

trial1rest = trial1rest';
trial2rest = trial2rest';
trial3rest = trial3rest';
trial4rest = trial4rest';
trialsrest = [trial1rest;trial2rest;trial3rest;trial4rest];

%create an additional 5th pupil channel (extra row) that takes the mean of the 4 trials
trialscoco(5,:) = mean(trialscoco(1:4,:));
trialsrest(5,:) = mean(trialsrest(1:4,:));

conditions = { 'coco'};
for condition=1:length(conditions)
    
    fig = figure('units','centimeters','outerposition',[0 0 25 20],'Color',[1 1 1]);
    hold on
    plot(trialscoco(5,:));
    h = xline(6000,'r','Trial onset'); %this line is added with the -10 + 10 range
    h1 = xline(66000,'r','Trial offset'); %this line is added with the -10 + 10 range
    title([subject ' coco - ' ' average time course']);
end 

%Add labels etc to the plot
    xlim([0 69000]) %xlim([0 60000])
    xticks([0 6000 16000 26000 36000 46000 56000 66000 69000]) %xticks([0 10000 20000 30000 40000 50000 60000])
    xticklabels({'-6','0','10','20','30','40','50','60','+3'}) %xticklabels(}'0','10','20','30','40','50','60'})
    xlabel('Time (s)');
    ylabel('Pupil dilation');

 %Save the plot
    if strcmp(conf.output.save,'yes')
        fig_name = [ subject '_coco_' conditions{condition} '_time_courses' ]; % fig_name = [ subject '_coco_' conditions{condition} '_time_courses' ];
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'jpg');
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'fig');
    end
        
%Save the pupil values time courses coco (data)
if strcmp(conf.output.save,'yes')
    pupil_coco_time_courses = trialscoco; %  pupil_coco_time_courses = trialscoco;
    save(fullfile(conf.output.datadir,[ subject '_coco_pupil_values_trial_time_courses']),'pupil_coco_time_courses'); % save(fullfile(conf.output.datadir,[ subject '_coco_pupil_values_trial_time_courses']),'pupil_coco_time_courses');
end
   
conditions = { 'rest'};
for condition=1:length(conditions)
    
    fig = figure('units','centimeters','outerposition',[0 0 25 20],'Color',[1 1 1]);
    hold on
    plot(trialsrest(5,:));
    h = xline(6000,'r','Trial onset'); %this line is added with the -10 + 10 range
    h1 = xline(66000,'r','Trial offset');%this line is added with the -10 + 10 range
    title([subject ' rest - ' ' average time course']);
end 

%Add labels etc to the plot
    xlim([0 69000]) %xlim([0 60000])
    xticks([0 6000 16000 26000 36000 46000 56000 66000 69000]) %xticks([0 10000 20000 30000 40000 50000 60000])
    xticklabels({'-6','0','10','20','30','40','50','60','+3'}) %xticklabels(}'0','10','20','30','40','50','60'})
    xlabel('Time (s)');
    ylabel('Pupil dilation');
   
    %Save the plot
    if strcmp(conf.output.save,'yes')
        fig_name = [ subject '_rest_' conditions{condition} '_time_courses' ];   %fig_name = [ subject '_rest_' conditions{condition} '_time_courses' ];
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'jpg');
        saveas(gcf, fullfile(conf.output.figdir,fig_name), 'fig');
    end
    
%Save the pupil values time courses rest (data)
if strcmp(conf.output.save,'yes')
    pupil_rest_time_courses = trialsrest; % pupil_rest_time_courses = trialsrest;
    save(fullfile(conf.output.datadir,[ subject '_rest_pupil_values_trial_time_courses']),'pupil_rest_time_courses'); %save(fullfile(conf.output.datadir,[ subject '_rest_pupil_values_trial_time_courses']),'pupil_rest_time_courses');
end



%% Create variable with group average pupil values for coco & rest

subject = '018';

cd(['/project/3024005.01/Analysis/Pupil/Results/' subject '/Data']);

%coco
load([subject '_coco_pupil_values_trial_time_courses.mat']);
pupil_time_courses_coco_n17(17,:) = pupil_coco_time_courses(5,:);

%rest
load([subject '_rest_pupil_values_trial_time_courses.mat']);
pupil_time_courses_rest_n17(17,:) = pupil_rest_time_courses(5,:);

cd('/project/3024005.01/Analysis/Pupil/Results/group-level (60 sec)/coco');
save('pupil_time_courses_coco_n17.mat','pupil_time_courses_coco_n17') 
save('pupil_time_courses_rest_n17.mat','pupil_time_courses_rest_n17') 



%% Plot coco & rest average time-course N = 17
cd('/project/3024005.01/Analysis/Pupil/Results/group-level (60 sec)/coco');

%Calculate the average across 17 subjects for coco
average_pupil_time_courses_coco_n17 = mean(pupil_time_courses_coco_n17);
save('average_pupil_time_courses_coco_n17.mat','average_pupil_time_courses_coco_n17') 

%Calculate the average across 17 subjects for rest
average_pupil_time_courses_rest_n17 = mean(pupil_time_courses_rest_n17);
save('average_pupil_time_courses_rest_n17.mat','average_pupil_time_courses_rest_n17') 

fig = figure('units','centimeters','outerposition',[0 0 25 20],'Color',[1 1 1]);
        
plot(average_pupil_time_courses_coco_n17(1,1:69000)); %Plot average coco power across 17 subjects in 1 line in 1 figure.
        xlim([0 69000])   %x = linspace(0,60000);
        xticks([0 6000 16000 26000 36000 46000 56000 66000 69000]) %xticks([0 10000 20000 30000 40000 50000 60000])
        xticklabels({'-6','0','10','20','30','40','50','60','+3'})
        h = xline(6000,'r','Trial onset'); %this line is added with the -10 + 10 range
        h1 = xline(66000,'r','Trial offset'); %this line is added with the -10 + 10 range
        xlabel('Time (s)');
        ylabel('Pupil diameter');
        title(['coco average time course (N = 17)']); 
    

        fig_name = ['coco_average_time_course_n17'];
        saveas(gcf, fullfile('/project/3024005.01/Analysis/Pupil/Results/group-level (60 sec)/coco',fig_name), 'jpg');
        saveas(gcf, fullfile('/project/3024005.01/Analysis/Pupil/Results/group-level (60 sec)/coco',fig_name), 'fig');

hold on
plot(average_pupil_time_courses_rest_n17(1,1:69000)); %Plot average rest power across 17 subjects in 1 line in 1 figure.
leg2 = legend('coco','rest');
set(leg2, 'Position', [.6 .8 .11 .08], 'Color', [1 1 1],'FontSize',10);
 
 
        xlim([0 69000])   %x = linspace(0,60000);
        xticks([0 6000 16000 26000 36000 46000 56000 66000 69000]) %xticks([0 10000 20000 30000 40000 50000 60000])
        xticklabels({'-6','0','10','20','30','40','50','60','+3'})
        h = xline(6000,'r','Trial onset'); %this line is added with the -10 + 10 range
        h1 = xline(66000,'r','Trial offset'); %this line is added with the -10 + 10 range
        xlabel('Time (s)');
        ylabel('Pupil diameter');
        title(['coco & rest average time course (N = 17)']); 
    

        fig_name = ['rest_average_time_course_n17'];
        saveas(gcf, fullfile('/project/3024005.01/Analysis/Pupil/Results/group-level (60 sec)/coco',fig_name), 'jpg');
        saveas(gcf, fullfile('/project/3024005.01/Analysis/Pupil/Results/group-level (60 sec)/coco',fig_name), 'fig');
 
 