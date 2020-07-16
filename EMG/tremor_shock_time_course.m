%plot time-course individual subjects for shock

%transpose data
sub_1 = sub_1';
sub_2 = sub_2';
sub_3 = sub_3';
sub_4 = sub_4';
sub_5 = sub_5';
sub_6 = sub_6';
sub_7 = sub_7';
sub_8 = sub_8';
sub_10 = sub_10';
sub_11 = sub_11';
sub_12 = sub_12';
sub_13 = sub_13';
sub_14 = sub_14';
sub_15 = sub_15';
sub_16 = sub_16';
sub_17 = sub_17';
sub_18 = sub_18';


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
c6 = [0.40 0 0]; c7 = [0.75 0 0]; c8 = [0 0.40 0]; c9 = [0 0.75 0]; c10 = [.5 .1 .5];
c66 = [0.40 0.1 0.1]; c77 = [0.75 0.1 0.1]; c88 = [0.1 0.40 0.1]; c99 = [0.1 0.75 0.1]; c100 = [.5 .2 .5];

c11 = [0.2 0 0]; c12 = [0.9 0 0]; c13 = [0 0.3 0]; c14 = [0 0.65 0]; c15 = [.4 .1 .4];
c111 = [0.2 0.1 0.1]; c122 = [0.9 0.1 0.1]; c133 = [0.1 0.3 0.1]; c144 = [0.1 0.65 0.1]; c155 = [.4 .2 .4];

c16 = [0.15 0 0]; c17 = [0.3 0 0]; 
c166 = [0.15 0.1 0.1]; c177 = [0.3 0.1 0.1];

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
 
 
nSETS = 18; %5
rNcN = size(sub_1);
DP_REP = size(sub_2);
safe_rest = size(sub_3);
safe_odd = size(sub_4);
shock = size(sub_5);
sub6 = size(sub_6);
sub7 = size(sub_7);
sub8 = size(sub_8);
sub10 = size(sub_10);
sub11 = size(sub_11);
sub12 = size(sub_12);
sub13 = size(sub_13);
sub14 = size(sub_14);
sub15 = size(sub_15);
sub16 = size(sub_16);
sub17 = size(sub_17);
sub18 = size(sub_18);
AveOver = 1;
DATARATE = 1;
t = 1;
 
 
	%==============================================%
	MuDATA=sub_1; repDATA=rNcN(2);
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
 
[ph1, po1] = boundedline(XT_Mu,YT_Mu, ET_Mu,'cmap',c1,'alpha','transparency', 0.1);
	hold on
 
	%==============================================%
	MuDATA=sub_2; repDATA=DP_REP(2);
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
 
[ph2, po2] = boundedline(XT_Mu,YT_Mu, ET_Mu,'cmap',c2,'alpha','transparency', 0.1);

%==============================================%
	MuDATA=sub_3 ; repDATA=safe_rest(2);
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
	MuDATA=sub_4; repDATA=safe_odd(2);
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
	MuDATA=sub_5; repDATA=shock(2);
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
	
   %==============================================%
	MuDATA=sub_6; repDATA=sub_6(2);
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
 
[ph6, po6] = boundedline(XT_Mu,YT_Mu, ET_Mu,'cmap',c6,'alpha','transparency', 0.1);
%	    axis tight; hold on;


%==============================================%
	MuDATA=sub_7; repDATA=sub_7(2);
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
 
[ph7, po7] = boundedline(XT_Mu,YT_Mu, ET_Mu,'cmap',c7,'alpha','transparency', 0.1);
%	    axis tight; hold on;


%==============================================%
	MuDATA=sub_8; repDATA=sub_8(2);
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
 
[ph8, po8] = boundedline(XT_Mu,YT_Mu, ET_Mu,'cmap',c8,'alpha','transparency', 0.1);
%	    axis tight; hold on;

%==============================================%
	MuDATA=sub_10; repDATA=sub_10(2);
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
 
[ph10, po10] = boundedline(XT_Mu,YT_Mu, ET_Mu,'cmap',c9,'alpha','transparency', 0.1);
%	    axis tight; hold on;


%==============================================%
	MuDATA=sub_11; repDATA=sub_11(2);
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
 
[ph11, po11] = boundedline(XT_Mu,YT_Mu, ET_Mu,'cmap',c10,'alpha','transparency', 0.1);
%	    axis tight; hold on;


%==============================================%
	MuDATA=sub_12; repDATA=sub_12(2);
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
 
[ph12, po12] = boundedline(XT_Mu,YT_Mu, ET_Mu,'cmap',c11,'alpha','transparency', 0.1);
%	    axis tight; hold on;

%==============================================%
	MuDATA=sub_13; repDATA=sub_13(2);
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
 
[ph13, po13] = boundedline(XT_Mu,YT_Mu, ET_Mu,'cmap',c12,'alpha','transparency', 0.1);
%	    axis tight; hold on;


%==============================================%
	MuDATA=sub_14; repDATA=sub_14(2);
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
 
[ph14, po14] = boundedline(XT_Mu,YT_Mu, ET_Mu,'cmap',c13,'alpha','transparency', 0.1);
%	    axis tight; hold on;


%==============================================%
	MuDATA=sub_15; repDATA=sub_15(2);
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
 
[ph15, po15] = boundedline(XT_Mu,YT_Mu, ET_Mu,'cmap',c14,'alpha','transparency', 0.1);
%	    axis tight; hold on;

%==============================================%
	MuDATA=sub_16; repDATA=sub_16(2);
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
 
[ph16, po16] = boundedline(XT_Mu,YT_Mu, ET_Mu,'cmap',c15,'alpha','transparency', 0.1);
%	    axis tight; hold on;

%==============================================%
	MuDATA=sub_17; repDATA=sub_17(2);
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
 
[ph17, po17] = boundedline(XT_Mu,YT_Mu, ET_Mu,'cmap',c16,'alpha','transparency', 0.1);
%	    axis tight; hold on;

%==============================================%
	MuDATA=sub_18; repDATA=sub_18(2);
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
 
[ph18, po18] = boundedline(XT_Mu,YT_Mu, ET_Mu,'cmap',c17,'alpha','transparency', 0.1);
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
   set(ph1,'LineStyle','-','Color',c1,'LineWidth',2,...
        'Marker','none','MarkerSize',MS2,'MarkerEdgeColor',c1);
   set(ph2,'LineStyle','-','Color',c2,'LineWidth',2,...
      'Marker','none','MarkerSize',MS2,'MarkerEdgeColor',c2);
   set(ph3,'LineStyle','-','Color',c3,'LineWidth',2,...
      'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c3);
   set(ph4,'LineStyle','-','Color',c4,'LineWidth',2,...
      'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c4);
    set(ph5,'LineStyle','-','Color',c5,'LineWidth',2,...
      'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c5);
   set(ph6,'LineStyle','-','Color',c6,'LineWidth',2,...
        'Marker','none','MarkerSize',MS2,'MarkerEdgeColor',c6);
   set(ph7,'LineStyle','-','Color',c7,'LineWidth',2,...
      'Marker','none','MarkerSize',MS2,'MarkerEdgeColor',c7);
   set(ph8,'LineStyle','-','Color',c8,'LineWidth',2,...
      'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c8);
    set(ph10,'LineStyle','-','Color',c9,'LineWidth',2,...
      'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c9);
   set(ph11,'LineStyle','-','Color',c10,'LineWidth',2,...
        'Marker','none','MarkerSize',MS2,'MarkerEdgeColor',c10);
   set(ph12,'LineStyle','-','Color',c11,'LineWidth',2,...
      'Marker','none','MarkerSize',MS2,'MarkerEdgeColor',c11);
   set(ph13,'LineStyle','-','Color',c12,'LineWidth',2,...
      'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c12);
   set(ph14,'LineStyle','-','Color',c13,'LineWidth',2,...
      'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c13);
    set(ph15,'LineStyle','-','Color',c14,'LineWidth',2,...
      'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c14);
   set(ph16,'LineStyle','-','Color',c15,'LineWidth',2,...
        'Marker','none','MarkerSize',MS2,'MarkerEdgeColor',c15);
   set(ph17,'LineStyle','-','Color',c16,'LineWidth',2,...
      'Marker','none','MarkerSize',MS2,'MarkerEdgeColor',c16);
   set(ph18,'LineStyle','-','Color',c17,'LineWidth',2,...
      'Marker','none','MarkerSize',MS2,'MarkerEdgeColor',c17);
   
 
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
        ylabel('tremor power');%ylabel('Tremor power (log µV^2)');
        title(['Tremor power time-course during shock N = 17']); 
        
        title('');
        set(gca,'Box','off','TickLength',[.01 .01], ...
        'XMinorTick','off','YMinorTick','off', ...
        'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',2);
        ylabel('\fontsize{20} tremor power ','Color', [0 0 0]);
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