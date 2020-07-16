 function tremor_2nd_level(subject_list)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Adjust this part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subject_list = {'001','002','003','004','005','006','007','008','010','011','012','013','014','015','016','017','018'};
tasks = { 'coco', 'TOS' };
conditions = { { 'Coco', 'Rest' }, { 'Threat-fix', 'Threat-odd', 'Safe-fix', 'Safe-odd', 'Shock' } };%Conditions of each task - ordered
trials = { { 4, 4 }, { 4, 4, 4, 4, 1 } }; %Number of trials for each condition of each task


%% Define directories and filenames
outputdir = '/project/3024005.01/Analysis/Tremor/Results_correct_time_course/Group60s';
inputbasedir = '/project/3024005.01/Analysis/Tremor/Results_correct_time_course/'; %Folder in which all subject subfolder with 1st level results are in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist(outputdir)
    mkdir(outputdir);
end

for task=1:length(tasks)
    for sb=1:length(subject_list)
        
        %% Load average power values and combine data from all subjects
        inputsubjdir = [ inputbasedir subject_list{sb} '/Data/' ];
        
        filename = [ subject_list{sb} '_' tasks{task} '_average_log_power_per_trial' ];
        load( fullfile(inputsubjdir,filename) );
        all_power_vals(sb,:) = power_values;

        %% Calculate individual averages for all conditions
        start_trial = 1;
        for condition=1:length(conditions{task})
            avg_power_vals(sb,condition) = mean(all_power_vals(sb,start_trial:start_trial+(trials{task}{condition}-1)));
            start_trial = start_trial + trials{task}{condition};
        end
        
        %% Load trial power values (time courses) and combine data from all subjects
        filename = [ subject_list{sb} '_' tasks{task} '_log_power_trial_time_courses' ];
        load( fullfile(inputsubjdir,filename) );
        
        for i=1:size(power_time_courses,2)
            all_cond_power_time_courses{i}(sb,:) = power_time_courses{i}(end,:);
        end
        
    end

    %% Plot the group results (average log tremor power per condition)
    Ylim = [ floor(min(avg_power_vals(:)))-2 ceil(max(avg_power_vals(:)))+2 ];
    figure;
    bar(1:size(avg_power_vals,2),mean(avg_power_vals,1));
    hold on
    errorbar( mean(avg_power_vals,1), std(avg_power_vals,1)/sqrt(size(avg_power_vals,1)), '.b' )
    % grid on
    ylabel('Average tremor power (log µV^2)')
    title([ tasks{task} ' - Average log tremor power per condition (N = ' num2str(length(subject_list)) ')'])
    XTickLabel = conditions{task};
    set(gca, 'XTick',1:size(avg_power_vals,2));
    set(gca, 'XTickLabel', XTickLabel);
    ylim( Ylim )
    
    %Add individual data points
    for sb=1:size(avg_power_vals,1)
        if rem(sb,2) == 0
            plot(1:size(avg_power_vals,2),avg_power_vals(sb,:),'-o','Color',[ (1-0.08*(sb-9))*floor(sb/9) (1-0.08*sb)*rem(floor(sb/9)+1,2) 0.5 ])
        else
            plot(1:size(avg_power_vals,2),avg_power_vals(sb,:),'--o','Color',[ (1-0.08*(sb-9))*floor(sb/9) (1-0.08*sb)*rem(floor(sb/9)+1,2) 0.5 ])
        end
    end
    
    fig_name = [ tasks{task} '_average_log_tremor_power_per_condition_N' num2str(length(subject_list)) ];
    saveas(gcf, fullfile(outputdir,fig_name), 'jpg');
    saveas(gcf, fullfile(outputdir,fig_name), 'fig');
    
    %% Plot Coco and TOS tremor time courses per condition
    for condition=1:size(all_cond_power_time_courses,2)
        
        avg_cond_power_time_courses(condition,:) = mean(all_cond_power_time_courses{condition},1);
        
        fig = figure('units','centimeters','outerposition',[0 0 25 20],'Color',[1 1 1]);
        plot(avg_cond_power_time_courses(condition,:));
        xlim([0 69000])
        xticks([0 6000 16000 26000 36000 46000 56000 66000 69000])
        xticklabels({'-6','0','10','20','30','40','50','60','+3'})
        h = xline(6000,'r','Trial onset'); %this line is added withthe -10 + 10 range
        h1 = xline(66000,'r','Trial offset'); %this line is added with the -10 + 10 range
        xlabel('Time (s)');
        ylabel('Tremor power (log µV^2)');
        title([ conditions{task}{condition} ' average group time course (N = ' num2str(length(subject_list)) ')' ]);
        
        fig_name = [ conditions{task}{condition} ' average group time course (N = ' num2str(length(subject_list)) ')' ];
        saveas(gcf, fullfile(outputdir,fig_name), 'jpg');
        saveas(gcf, fullfile(outputdir,fig_name), 'fig');
        
    end
    
end


%% further plotting of the time-courses for COCO & TOS
%plot coco & rest together
plot(avg_cond_power_time_courses(1,1:76000));
hold on
plot(avg_cond_power_time_courses(2,1:76000)); %Plot average rest power across 17 subjects in 1 line in 1 figure.

plot(average_tremor_time_courses_threat_rest_n17_min8plus8(:,:));
hold on
plot(average_tremor_time_courses_threat_odd_n17_min8plus8(:,:));
plot(average_tremor_time_courses_safe_rest_n17_min8plus8(:,:));
plot(average_tremor_time_courses_safe_odd_n17_min8plus8(:,:));
plot(average_tremor_time_courses_shock_n17_min8plus8(:,:));
        xlim([0 76000])
        xticks([0 8000 18000 28000 38000 48000 58000 68000 76000])
        xticklabels({'-8','0','10','20','30','40','50','60','+8'})
        h = xline(8000,'r','Trial onset'); %this line is added withthe -10 + 10 range
        h1 = xline(68000,'r','Trial offset'); %this line is added with the -10 + 10 range
        xlabel('Time (s)');
        ylabel('Tremor power (log µV^2)');
        leg2 = legend('safe-rest','safe-odd');%leg2 = legend('coco','rest');
        set(leg2, 'Position', [.5 .8 .11 .08], 'Color', [1 1 1],'FontSize',10);
        title('TOS average group time course (N = 17)');

%% boundedline toolbox for plotting the time-course COCO & TOS

%Reshape the data for plotting purposes

%coco
tremor_time_courses_coco_n17 = tremor_time_courses_coco_n17';

%rest
tremor_time_courses_rest_n17 = tremor_time_courses_rest_n17';


%Reshape the data for plotting purposes

%threat-rest
tremor_time_courses_threat_rest_n17 = tremor_time_courses_threat_rest_n17';

%threat-odd
tremor_time_courses_threat_odd_n17 = tremor_time_courses_threat_odd_n17';

%safe-rest
tremor_time_courses_safe_rest_n17 = tremor_time_courses_safe_rest_n17';

%safe-odd
tremor_time_courses_safe_odd_n17 = tremor_time_courses_safe_odd_n17';

%shock
tremor_time_courses_shock_n17 = tremor_time_courses_shock_n17';

%pupil
pupil_time_courses_coco_n17 = pupil_time_courses_coco_n17';
pupil_time_courses_rest_n17 = pupil_time_courses_rest_n17';

pupil_time_courses_threat_rest_n17 = pupil_time_courses_threat_rest_n17';
pupil_time_courses_threat_odd_n17 = pupil_time_courses_threat_odd_n17';
pupil_time_courses_safe_rest_n17 = pupil_time_courses_safe_rest_n17';
pupil_time_courses_safe_odd_n17 = pupil_time_courses_safe_odd_n17';
pupil_time_courses_shock_n17 = pupil_time_courses_shock_n17';


%shock figure with 17 individual lines

clc; close all;
scsz = get(0,'ScreenSize');
pos1 = [scsz(3)/10  scsz(4)/10  scsz(3)/1.5  scsz(4)/1.5];
fig55 = figure(55);
set(fig55,'Renderer','OpenGL','Units','pixels','OuterPosition',pos1,'Color',[1,1,1])
%set(figure(1),'Color',[1 1 1]);
%-------------------------------------------------------------
 
 
%-------------------------------------------------------------
%coco
c1 = [0.55 0 0]; c2 = [0 0.55 0];
c11 = [0.55 0.1 0.1]; c22 = [0.1 0.55 0.1];
%TOS
c1 = [0.55 0 0]; c2 = [0.85 0 0]; c3 = [0 0.55 0]; c4 = [0 0.85 0]; c5 = [.8 .1 .8];
c11 = [0.55 0.1 0.1]; c22 = [0.85 0.1 0.1]; c33 = [0.1 0.55 0.1]; c44 = [0.1 0.85 0.1]; c55 = [.8 .2 .8];
%c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6]; c5= [.01 .9 .01];
%c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7]; c55=[.01 .9 .01];
applered= [.9 .2 .2]; oceanblue= [.2 .4 .6]; neongreen = [.1 .9 .1];
liteblue = [.2 .9 .9]; hotpink=[.9 .1 .9]; c11 = 'none';
MSz = 7;
ax = [.10 .10 .85 .85];
%-------------------------------------------------------------
 
 
 
% Assuming each line represents the average of 3 columns of data...
 
% dataset1 : 55x3 (row x col) dataset
 
% dataset2 : 55x3 (row x col) dataset
 
 
%===========================================================%
 
% Massage Data
 
%===========================================================%
 
 
nSETS = 2; %5
rNcN = size(pupil_time_courses_coco_n17);
DP_REP = size(pupil_time_courses_rest_n17);
safe_rest = size(tremor_time_courses_safe_rest_n17 );
safe_odd = size(tremor_time_courses_safe_odd_n17);
shock = size(tremor_time_courses_shock_n17);
AveOver = 1;
DATARATE = 1;
t = 1;
 
 
	%==============================================%
	MuDATA=pupil_time_courses_coco_n17; repDATA=rNcN(2);
	%------------------------------
	Mu = mean(MuDATA,2)';       Sd = std(MuDATA,0,2)';      Se = Sd./sqrt(repDATA);
	y_Mu = Mu;                  x_Mu = 1:(size(y_Mu,2));    e_Mu = Se;
	xx_Mu = 1:0.1:max(x_Mu);
	yy_Mu = interp1(x_Mu,y_Mu,xx_Mu,'pchip');
	ee_Mu = interp1(x_Mu,e_Mu,xx_Mu,'pchip');
	p_Mu = polyfit(x_Mu,Mu,1);
	x2_Mu = 1:0.1:max(x_Mu);	y2_Mu = polyval(p_Mu,x2_Mu);
	XT_Mu = xx_Mu';				YT_Mu = yy_Mu';		ET_Mu = ee_Mu';
	%==============================================%
 
 
	%hax = axes('Position',ax);
 
[ph1, po1] = boundedline(XT_Mu,YT_Mu, ET_Mu,'cmap',c1,'alpha','transparency', 0.2);
	hold on
 
	%==============================================%
	MuDATA=pupil_time_courses_rest_n17; repDATA=DP_REP(2);
	%------------------------------
	Mu = mean(MuDATA,2)';		Sd = std(MuDATA,0,2)';		Se = Sd./sqrt(repDATA);
	y_Mu = Mu;				x_Mu = 1:(size(y_Mu,2));	e_Mu = Se;
	xx_Mu = 1:0.1:max(x_Mu);
	 yy_Mu = spline(x_Mu,y_Mu,xx_Mu);	% ee_Mu = spline(x_Mu,e_Mu,xx_Mu);
	yy_Mu = interp1(x_Mu,y_Mu,xx_Mu,'pchip');
	ee_Mu = interp1(x_Mu,e_Mu,xx_Mu,'pchip');
	p_Mu = polyfit(x_Mu,Mu,1);
	x2_Mu = 1:0.1:max(x_Mu);	y2_Mu = polyval(p_Mu,x2_Mu);
	XT_Mu = xx_Mu';				YT_Mu = yy_Mu';		ET_Mu = ee_Mu';
	%==============================================%
 
[ph2, po2] = boundedline(XT_Mu,YT_Mu, ET_Mu,'cmap',c2,'alpha','transparency', 0.2);

%==============================================%
	MuDATA=tremor_time_courses_safe_rest_n17 ; repDATA=safe_rest(2);
	%------------------------------
	Mu = mean(MuDATA,2)';		Sd = std(MuDATA,0,2)';		Se = Sd./sqrt(repDATA);
	y_Mu = Mu;				x_Mu = 1:(size(y_Mu,2));	e_Mu = Se;
	xx_Mu = 1:0.1:max(x_Mu);
	 yy_Mu = spline(x_Mu,y_Mu,xx_Mu);	% ee_Mu = spline(x_Mu,e_Mu,xx_Mu);
	yy_Mu = interp1(x_Mu,y_Mu,xx_Mu,'pchip');
	ee_Mu = interp1(x_Mu,e_Mu,xx_Mu,'pchip');
	p_Mu = polyfit(x_Mu,Mu,1);
	x2_Mu = 1:0.1:max(x_Mu);	y2_Mu = polyval(p_Mu,x2_Mu);
	XT_Mu = xx_Mu';				YT_Mu = yy_Mu';		ET_Mu = ee_Mu';
	%==============================================%
 
[ph3, po3] = boundedline(XT_Mu,YT_Mu, ET_Mu,'cmap',c3,'alpha','transparency', 0.1);

%==============================================%
	MuDATA=tremor_time_courses_safe_odd_n17; repDATA=safe_odd(2);
	%------------------------------
	Mu = mean(MuDATA,2)';		Sd = std(MuDATA,0,2)';		Se = Sd./sqrt(repDATA);
	y_Mu = Mu;				x_Mu = 1:(size(y_Mu,2));	e_Mu = Se;
	xx_Mu = 1:0.1:max(x_Mu);
	 yy_Mu = spline(x_Mu,y_Mu,xx_Mu);	% ee_Mu = spline(x_Mu,e_Mu,xx_Mu);
	yy_Mu = interp1(x_Mu,y_Mu,xx_Mu,'pchip');
	ee_Mu = interp1(x_Mu,e_Mu,xx_Mu,'pchip');
	p_Mu = polyfit(x_Mu,Mu,1);
	x2_Mu = 1:0.1:max(x_Mu);	y2_Mu = polyval(p_Mu,x2_Mu);
	XT_Mu = xx_Mu';				YT_Mu = yy_Mu';		ET_Mu = ee_Mu';
	%==============================================%
 
[ph4, po4] = boundedline(XT_Mu,YT_Mu, ET_Mu,'cmap',c4,'alpha','transparency', 0.1);

%==============================================%
	MuDATA=tremor_time_courses_shock_n17; repDATA=shock(2);
	%------------------------------
	Mu = mean(MuDATA,2)';		Sd = std(MuDATA,0,2)';		Se = Sd./sqrt(repDATA);
	y_Mu = Mu;				x_Mu = 1:(size(y_Mu,2));	e_Mu = Se;
	xx_Mu = 1:0.1:max(x_Mu);
	 yy_Mu = spline(x_Mu,y_Mu,xx_Mu);	% ee_Mu = spline(x_Mu,e_Mu,xx_Mu);
	yy_Mu = interp1(x_Mu,y_Mu,xx_Mu,'pchip');
	ee_Mu = interp1(x_Mu,e_Mu,xx_Mu,'pchip');
	p_Mu = polyfit(x_Mu,Mu,1);
	x2_Mu = 1:0.1:max(x_Mu);	y2_Mu = polyval(p_Mu,x2_Mu);
	XT_Mu = xx_Mu';				YT_Mu = yy_Mu';		ET_Mu = ee_Mu';
	%==============================================%
 
[ph5, po5] = boundedline(XT_Mu,YT_Mu, ET_Mu,'cmap',c5,'alpha','transparency', 0.1);
%	    axis tight; hold on;
	
    %leg2 = legend([ph1,ph2,ph3,ph4,ph5],{'threat','threat-odd','safe','safe-odd','shock'});
    leg2 = legend([ph1,ph2],{'coco','rest'});
    set(leg2, 'Position', [.65 .16 .11 .08], 'Color', [1 1 1],'FontSize',14);
	
%leg1 = legend([ph1],{'avg pupil diameter'});
%    set(leg1, 'Position', [.15 .85 .11 .08], 'Color', [1 1 1],'FontSize',10);
 
	%------ Legend &amp; Tick Labels -------
	%if verLessThan('matlab', '8.3.1');
%		xt = roundn((get(gca,'XTick')).*AveOver*DATARATE.*(t)./(60),0);
%		set(gca,'XTickLabel', sprintf('%.0f|',xt))
%	else
%        hax2 = (get(gca));
%        hax2.XTick = ([0,10000,20000,30000,40000,50000,60000]);%
%		xt = hax2.XTick;%
%		xtt = roundn(xt*AveOver*DATARATE*(t)/(60),0);
%		hax2.XTickLabel = xtt;
%        hax2.XLim = ([0 60000]);%
%	end
	%------
 
 
        MS1 = 5; MS2 = 2; %MS3 = 1; MS4 = 4;
   set(ph1,'LineStyle','-','Color',c1,'LineWidth',3,...
        'Marker','none','MarkerSize',MS2,'MarkerEdgeColor',c1);
   set(ph2,'LineStyle','-','Color',c2,'LineWidth',3,...
      'Marker','none','MarkerSize',MS2,'MarkerEdgeColor',c2);
   set(ph3,'LineStyle','-','Color',c3,'LineWidth',3,...
      'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c3);
   set(ph4,'LineStyle','-','Color',c4,'LineWidth',3,...
      'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c4);
    set(ph5,'LineStyle','-','Color',c5,'LineWidth',3,...
      'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c5);
   
 
    %hTitle  = title ('\fontsize{12} threat-rest & threat-oddball average time-courses N = 17');
    %hXLabel = xlabel('\fontsize{12} Time (s)');
    %hYLabel = ylabel('\fontsize{12} Tremor power (log µV^2)');
    %set(gca,'FontName','Helvetica','FontSize',12);
    %set([hTitle, hXLabel, hYLabel],'FontName','Century Gothic');
    %set(gca,'Box','off','TickDir','out','TickLength',[.01 .01], ...
    %'XMinorTick','off','YMinorTick','on','YGrid','on', ...
    %'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',2);
 
    %------
    % Extra axis for boxing
    %haxes1 = gca; % handle to axes
	%haxes1_pos = get(haxes1,'Position'); % store position of first axes
	%haxes2 = axes('Position',haxes1_pos,'Color','none',...
	%	  'XAxisLocation','top','YAxisLocation','right');
    %set(gca,'Box','off','TickDir','out','TickLength',[.01 .01], ...
	%'XMinorTick','off','YMinorTick','off','XGrid','off','YGrid','off', ...
	%'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',2, ...
    %'XTick', [], 'YTick', []);
    %------
 
% for COCO    
        xlim([0 69000])
        xticks([0 6000 16000 26000 36000 46000 56000 66000 69000])
        xticklabels({'-6','0','10','20','30','40','50','60','+3'})
        h = xline(6000,'r','Trial onset','LineWidth',2); %this line is added withthe -10 + 10 range
        h1 = xline(66000,'r','Trial offset','LineWidth',2); %this line is added with the -10 + 10 range
        xlabel('Time (s)');
        ylabel('Pupil diameter (+/- SEM)');%ylabel('Tremor power (log µV^2)');
       % title(['coco & rest average time courses (N = 17)']); 
        
        title('');
        set(gca,'FontWeight', 'bold','FontName','Helvetica','FontSize',18);
        set(gca,'Box','off','TickLength',[.01 .01], ...
        'XMinorTick','off','YMinorTick','off', ...
        'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',2);
        ylabel('\fontsize{20} mean pupil diameter (+/- SEM) ','Color', [0 0 0]);
        xlabel('\fontsize{18} Time (s)','Color', [0 0 0]);
        %fig4_v2 = set(figure(1),'units','centimeters','outerposition',[0 0 25 20],'Color',[1 1 1]);
 % for TOS
        xlim([0 68000])  %x = linspace(0,72000); %x = linspace(0,60000);
        xticks([0 5000 15000 25000 35000 45000 55000 65000 68000]) %xticks([0 10000 20000 30000 40000 50000 60000])
        xticklabels({'-5','0','10','20','30','40','50','60','+3'})
        h = xline(5000,'r','Trial onset','LineWidth',2); %this line is added with the -10 + 10 range
        h1 = xline(65000,'r','Trial offset','LineWidth',2); %this line is added with the -10 + 10 range
        xlabel('Time (s)');
        ylabel('Pupil diameter (+/- SEM)');%ylabel('Tremor power (log µV^2)');
        title(['TOS average time courses (N = 17)']);
%===========================================================%
 
%% comparing pre vs. post-shock effect on tremor (-10 sec, +10 sec)

cd('/project/3024005.01/Analysis/Tremor/Results_correct_time_course/Group60s');
load('tremor_time_courses_shock_n17');

shock_admin_n17(17,:) = tremor_time_courses_shock_n17(17,[mean(30000:40000) mean(42000:52000)]);

shock_admin_n17_cell_v2{1,1} = shock_admin_n17(:,1);
shock_admin_n17_cell_v2{1,2} = shock_admin_n17(:,2);

%% comparing pre vs. post-shock effect on pupil diameter (-10 sec, +10 sec)

cd('/project/3024005.01/Analysis/Pupil/Results/group-level (60 sec)/TOS/0-60 data');
load('Pupil_time_courses_shock_n17');

shock_admin_n17(17,:) = Pupil_time_courses_shock_n17(17,[mean(25000:35000) mean(37000:47000)]);

shock_admin_n17_cell_pupil{1,1} = shock_admin_n17(:,1);
shock_admin_n17_cell_pupil{1,2} = shock_admin_n17(:,2);


%% fisher's z transformation on the averaged correlations
fishers_z(:,:) = atanh(all_correlations(:,:));
fishers_z_corr{1,1} = fishers_z(:,6);
fishers_z_corr{1,1}(:,2) = fishers_z(:,7);
fishers_z_corr{1,2} = fishers_z(:,1);
fishers_z_corr{1,2}(:,2) = fishers_z(:,3);
fishers_z_corr{1,3} = fishers_z(:,2);
fishers_z_corr{1,3}(:,2) = fishers_z(:,4);
fishers_z_corr{1,4} = fishers_z(:,5);

end




