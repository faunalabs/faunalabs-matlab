%% ftChapterFourStats.m
%
%   Dave Haas, 13 October 2021
%
%

clc;

global TAG_PATHS;

% named colors
blue = [0 0.4470 0.7410];           % '#0072BD';
red = [0.8500 0.3250 0.0980];       % '#D95319'	
goldenrod = [0.9290 0.6940 0.1250];	% '#EDB120'	
purple = [0.4940 0.1840 0.5560];    % '#7E2F8E'	
green = [0.4660 0.6740 0.1880];     % '#77AC30'	
cyan = [0.3010 0.7450 0.9330];      % '#4DBEEE'	
maroon = [0.6350 0.0780 0.1840];    % '#A2142F'
black = [0 0 0];                    % '#000000'
white = [1 1 1];                    % '#FFFFFF'
grey = [0.6 0.6 0.6];               % '#999999'

% set some default LED* colors...

led1Color = black;
led2Color = blue;
led3Color = red;
led4Color = goldenrod;

% ... and some default kinematic colors...

kxColor = blue;
kyColor = red;
kzColor = goldenrod;

%% Load the trial timetables for everything

load('/Users/dave/Documents/MATLAB/tursiopsTimetables.mat')

CC_OPTICS = CC_OPTICS_DATA;

%% Do Bland-Altman and correlation plots

% linear regression and Bland-Altman comparison of iHr 
% with linWsstHr, rotWsstHr, and ensWsstHr

figSizeBA = [ 100 100 1200 1200 ];

fig1 = figure('Position', figSizeBA, 'Color', white, 'NumberTitle', 'off');
ah1 = subplot(311);
ah2 = subplot(312);
ah3 = subplot(313);

[ rpc1, ah1, stats1 ] = BlandAltman(ah1, ...
    CC_OPTICS.snr_dbFiltLed2(CC_OPTICS.epoch==1), ...
    CC_OPTICS.snr_dbFiltLed3(CC_OPTICS.epoch==1), ...
    {'dbFiltLed2','dbFiltLed3'}, ...
    'Comparison of 1050 nm and 1200 nm optical signals · Epoch 1', ...
    {'filtLed2:filtLed3','LoA'}, ...
    'corrInfo', {'n','eq','r','r2','p'}, ...
    'baInfo', {'RPCnp','LOA','SD','ks'}, ...
    'data1mode', 'truth', ...
    'showFitCI', 'on', ...
    'legend', 'off', ...
    'axesLimits', 'tight', ...
    'baStatsMode','Non-parametric');

[ rpc2, ah2, stats2 ] = BlandAltman(ah2, ...
    CC_OPTICS.dbFiltLed2(CC_OPTICS.epoch==1), ...
    CC_OPTICS.dbFiltLed1(CC_OPTICS.epoch==1), ...
    {'dbFiltLed2','dbFiltLed1'}, ...
    'Comparison of 1050 nm and ambient 1 optical signals · Epoch 1', ...
    {'filtLed2:amb1','LoA'}, ...
    'corrInfo', {'n','eq','r','r2','p'}, ...
    'baInfo', {'RPCnp','LOA','SD','ks'}, ...
    'data1mode', 'truth', ...
    'showFitCI', 'on', ...
    'legend', 'off', ...
    'axesLimits', 'tight', ...
    'baStatsMode','Non-parametric');

[ rpc3, ah3, stats3 ] = BlandAltman(ah3, ...
    CC_OPTICS.dbFiltLed2(CC_OPTICS.epoch==1), ...
    CC_OPTICS.dbFiltLed4(CC_OPTICS.epoch==1), ...
	{'dbFiltLed2','dbFiltLed4'}, ...
    'Comparison of 1050 nm and ambient 2 optical signals · Epoch 1', ...
    {'filtLed2:amb2','LoA'}, ...
    'corrInfo', {'n','eq','r','r2','p'}, ...
    'baInfo', {'RPCnp','LOA','SD','ks'}, ...
    'data1mode', 'truth', ...
    'showFitCI', 'on', ...
    'legend', 'off', ...
    'axesLimits', 'tight', ...
    'baStatsMode','Non-parametric');

fprintf('Displaying correlation and Bland-Altman summary stats:\n');
fprintf('\n-------------\nlinWsstHr : iHr:\n');
disp(stats1);
fprintf('\n-------------\nrotWsstHr : iHr:\n');
disp(stats2);
fprintf('\n-------------\nensWsstHr : iHr:\n');
disp(stats3);

%%

%% Do Bland-Altman and correlation plots

% linear regression and Bland-Altman comparison of iHr 
% with linWsstHr, rotWsstHr, and ensWsstHr

figSizeBA = [ 100 100 1200 1200 ];

fig1 = figure('Position', figSizeBA, 'Color', white, 'NumberTitle', 'off');
ah1 = subplot(311);
ah2 = subplot(312);
ah3 = subplot(313);

[ rpc1, ah1, stats1 ] = BlandAltman(ah1, ...
    CC_OPTICS.dbFiltLed2(CC_OPTICS.epoch==3), ...
    CC_OPTICS.dbFiltLed3(CC_OPTICS.epoch==3), ...
    {'dbFiltLed2','dbFiltLed3'}, ...
    'Comparison of 1050 nm and 1200 nm optical signals · Epoch 3', ...
    {'filtLed2:filtLed3','LoA'}, ...
    'corrInfo', {'n','eq','r','r2','p'}, ...
    'baInfo', {'RPCnp','LOA','SD','ks'}, ...
    'data1mode', 'truth', ...
    'showFitCI', 'on', ...
    'legend', 'off', ...
    'axesLimits', 'tight', ...
    'baStatsMode','Non-parametric');

[ rpc2, ah2, stats2 ] = BlandAltman(ah2, ...
    CC_OPTICS.dbFiltLed2(CC_OPTICS.epoch==3), ...
    CC_OPTICS.dbFiltLed1(CC_OPTICS.epoch==3), ...
    {'dbFiltLed2','dbFiltLed1'}, ...
    'Comparison of 1050 nm and ambient 1 optical signals · Epoch 3', ...
    {'filtLed2:amb1','LoA'}, ...
    'corrInfo', {'n','eq','r','r2','p'}, ...
    'baInfo', {'RPCnp','LOA','SD','ks'}, ...
    'data1mode', 'truth', ...
    'showFitCI', 'on', ...
    'legend', 'off', ...
    'axesLimits', 'tight', ...
    'baStatsMode','Non-parametric');

[ rpc3, ah3, stats3 ] = BlandAltman(ah3, ...
    CC_OPTICS.dbFiltLed2(CC_OPTICS.epoch==3), ...
    CC_OPTICS.dbFiltLed4(CC_OPTICS.epoch==3), ...
	{'dbFiltLed2','dbFiltLed4'}, ...
    'Comparison of 1050 nm and ambient 2 optical signals · Epoch 3', ...
    {'filtLed2:amb2','LoA'}, ...
    'corrInfo', {'n','eq','r','r2','p'}, ...
    'baInfo', {'RPCnp','LOA','SD','ks'}, ...
    'data1mode', 'truth', ...
    'showFitCI', 'on', ...
    'legend', 'off', ...
    'axesLimits', 'tight', ...
    'baStatsMode','Non-parametric');

fprintf('Displaying correlation and Bland-Altman summary stats:\n');
fprintf('\n-------------\nlinWsstHr : iHr:\n');
disp(stats1);
fprintf('\n-------------\nrotWsstHr : iHr:\n');
disp(stats2);
fprintf('\n-------------\nensWsstHr : iHr:\n');
disp(stats3);


%% Do Bland-Altman and correlation plots

% linear regression and Bland-Altman comparison of iHr 
% with linWsstHr, rotWsstHr, and ensWsstHr

figSizeBA = [ 100 100 1200 1200 ];

fig1 = figure('Position', figSizeBA, 'Color', white, 'NumberTitle', 'off');
ah1 = subplot(311);
ah2 = subplot(312);
ah3 = subplot(313);

[ rpc1, ah1, stats1 ] = BlandAltman(ah1, ...
    CC_OPTICS.dbFiltLed2(CC_OPTICS.epoch==5), ...
    CC_OPTICS.dbFiltLed3(CC_OPTICS.epoch==5), ...
    {'dbFiltLed2','dbFiltLed3'}, ...
    'Comparison of 1050 nm and 1200 nm optical signals · Epoch 5', ...
    {'filtLed2:filtLed3','LoA'}, ...
    'corrInfo', {'n','eq','r','r2','p'}, ...
    'baInfo', {'RPCnp','LOA','SD','ks'}, ...
    'data1mode', 'truth', ...
    'showFitCI', 'on', ...
    'legend', 'off', ...
    'axesLimits', 'tight', ...
    'baStatsMode','Non-parametric');

[ rpc2, ah2, stats2 ] = BlandAltman(ah2, ...
    CC_OPTICS.dbFiltLed2(CC_OPTICS.epoch==5), ...
    CC_OPTICS.dbFiltLed1(CC_OPTICS.epoch==5), ...
    {'dbFiltLed2','dbFiltLed1'}, ...
    'Comparison of 1050 nm and ambient 1 optical signals · Epoch 5', ...
    {'filtLed2:amb1','LoA'}, ...
    'corrInfo', {'n','eq','r','r2','p'}, ...
    'baInfo', {'RPCnp','LOA','SD','ks'}, ...
    'data1mode', 'truth', ...
    'showFitCI', 'on', ...
    'legend', 'off', ...
    'axesLimits', 'tight', ...
    'baStatsMode','Non-parametric');

[ rpc3, ah3, stats3 ] = BlandAltman(ah3, ...
    CC_OPTICS.dbFiltLed2(CC_OPTICS.epoch==5), ...
    CC_OPTICS.dbFiltLed4(CC_OPTICS.epoch==5), ...
	{'dbFiltLed2','dbFiltLed4'}, ...
    'Comparison of 1050 nm and ambient 2 optical signals · Epoch 5', ...
    {'filtLed2:amb2','LoA'}, ...
    'corrInfo', {'n','eq','r','r2','p'}, ...
    'baInfo', {'RPCnp','LOA','SD','ks'}, ...
    'data1mode', 'truth', ...
    'showFitCI', 'on', ...
    'legend', 'off', ...
    'axesLimits', 'tight', ...
    'baStatsMode','Non-parametric');

fprintf('Displaying correlation and Bland-Altman summary stats:\n');
fprintf('\n-------------\nlinWsstHr : iHr:\n');
disp(stats1);
fprintf('\n-------------\nrotWsstHr : iHr:\n');
disp(stats2);
fprintf('\n-------------\nensWsstHr : iHr:\n');
disp(stats3);


%% Have a look at SNR plots with all epoch for tt21_134b

for whichId = 1:6

figure('Color', white);

plot(CC_OPTICS.Time(CC_OPTICS.altId == whichId & ...
    string(CC_OPTICS.opticsCue)=='3'), ...
    CC_OPTICS.snr_dbFiltLed2(CC_OPTICS.altId == whichId & ...
    string(CC_OPTICS.opticsCue)=='3'), ...
    'h', 'Color', green );
hold on;

plot(CC_OPTICS.Time(CC_OPTICS.altId == whichId & ...
    string(CC_OPTICS.opticsCue)=='z'), ...
    CC_OPTICS.snr_dbFiltLed2(CC_OPTICS.altId == whichId & ...
    string(CC_OPTICS.opticsCue)=='z'), ...
    'h', 'Color', red );

hold off;
xlabel('Time, local');
ylabel('dB');

switch(whichId)
    case 1
        titleTxt = 'tt21_128f';
    case 2
        titleTxt = 'tt21_134a';
    case 3
        titleTxt = 'tt21_134b';
    case 4
        titleTxt = 'tt21_141e';
    case 5
        titleTxt = 'tt21_142c';
    case 6
        titleTxt = 'tt21_142d';
end

title(titleTxt, 'Interpreter', 'none');
grid;
legend('ifH match', 'no ifH match')

end

%% stats for snr_dbFiltLed2

mean( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 1 & CC_OPTICS.epoch == 1 ) )
std( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 1 & CC_OPTICS.epoch == 1 ) )

mean( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 1 & CC_OPTICS.epoch == 3 ) )
std( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 1 & CC_OPTICS.epoch == 3 ) )

mean( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 2 & CC_OPTICS.epoch == 1 ) )
std( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 2 & CC_OPTICS.epoch == 1 ) )

mean( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 2 & CC_OPTICS.epoch == 3 ) )
std( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 2 & CC_OPTICS.epoch == 3 ) )

mean( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 2 & CC_OPTICS.epoch == 5 ) )
std( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 2 & CC_OPTICS.epoch == 5 ) )

mean( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 3 & CC_OPTICS.epoch == 1 ) )
std( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 3 & CC_OPTICS.epoch == 1 ) )

mean( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 3 & CC_OPTICS.epoch == 3 ) )
std( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 3 & CC_OPTICS.epoch == 3 ) )

mean( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 3 & CC_OPTICS.epoch == 5 ) )
std( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 3 & CC_OPTICS.epoch == 5 ) )

mean( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 4 & CC_OPTICS.epoch == 1 ) )
std( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 4 & CC_OPTICS.epoch == 1 ) )

mean( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 4 & CC_OPTICS.epoch == 3 ) )
std( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 4 & CC_OPTICS.epoch == 3 ) )

mean( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 4 & CC_OPTICS.epoch == 5 ) )
std( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 4 & CC_OPTICS.epoch == 5 ) )

mean( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 5 & CC_OPTICS.epoch == 1 ) )
std( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 5 & CC_OPTICS.epoch == 1 ) )

mean( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 5 & CC_OPTICS.epoch == 3 ) )
std( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 5 & CC_OPTICS.epoch == 3 ) )

mean( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 5 & CC_OPTICS.epoch == 5 ) )
std( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 5 & CC_OPTICS.epoch == 5 ) )

mean( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 6 & CC_OPTICS.epoch == 1 ) )
std( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 6 & CC_OPTICS.epoch == 1 ) )

mean( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 6 & CC_OPTICS.epoch == 3 ) )
std( CC_OPTICS.snr_dbFiltLed2( CC_OPTICS.altId == 6 & CC_OPTICS.epoch == 3 ) )

%% stats for snr_dbFiltLed3


mean( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 1 & CC_OPTICS.epoch == 1 ) )
std( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 1 & CC_OPTICS.epoch == 1 ) )

mean( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 1 & CC_OPTICS.epoch == 3 ) )
std( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 1 & CC_OPTICS.epoch == 3 ) )

mean( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 2 & CC_OPTICS.epoch == 1 ) )
std( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 2 & CC_OPTICS.epoch == 1 ) )

mean( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 2 & CC_OPTICS.epoch == 3 ) )
std( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 2 & CC_OPTICS.epoch == 3 ) )

mean( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 2 & CC_OPTICS.epoch == 5 ) )
std( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 2 & CC_OPTICS.epoch == 5 ) )

mean( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 3 & CC_OPTICS.epoch == 1 ) )
std( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 3 & CC_OPTICS.epoch == 1 ) )

mean( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 3 & CC_OPTICS.epoch == 3 ) )
std( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 3 & CC_OPTICS.epoch == 3 ) )

mean( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 3 & CC_OPTICS.epoch == 5 ) )
std( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 3 & CC_OPTICS.epoch == 5 ) )

mean( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 4 & CC_OPTICS.epoch == 1 ) )
std( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 4 & CC_OPTICS.epoch == 1 ) )

mean( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 4 & CC_OPTICS.epoch == 3 ) )
std( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 4 & CC_OPTICS.epoch == 3 ) )

mean( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 4 & CC_OPTICS.epoch == 5 ) )
std( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 4 & CC_OPTICS.epoch == 5 ) )

mean( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 5 & CC_OPTICS.epoch == 1 ) )
std( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 5 & CC_OPTICS.epoch == 1 ) )

mean( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 5 & CC_OPTICS.epoch == 3 ) )
std( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 5 & CC_OPTICS.epoch == 3 ) )

mean( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 5 & CC_OPTICS.epoch == 5 ) )
std( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 5 & CC_OPTICS.epoch == 5 ) )

mean( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 6 & CC_OPTICS.epoch == 1 ) )
std( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 6 & CC_OPTICS.epoch == 1 ) )

mean( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 6 & CC_OPTICS.epoch == 3 ) )
std( CC_OPTICS.snr_dbFiltLed3( CC_OPTICS.altId == 6 & CC_OPTICS.epoch == 3 ) )



%%

figure;
plot(CC_OPTICS.altId, CC_OPTICS.snr_dbFiltLed2, 'Color', blue);
hold on;
plot(CC_OPTICS.altId, CC_OPTICS.snr_dbFiltLed3, 'Color', red);
plot(CC_OPTICS.altId, CC_OPTICS.snr_dbFiltLed4, 'Color', goldenrod);
hold off;
grid;
legend;



%% Swarmcharts for SNR LED2 and SNR LED3, only epochs 1 and 5, all animals

figure;
swarmchart(CC_OPTICS.altId(string(CC_OPTICS.opticsCue)=='3'& ...
    CC_OPTICS.epoch ~= 3), ...
    CC_OPTICS.snr_dbFiltLed2( string(CC_OPTICS.opticsCue)=='3' & ...
    CC_OPTICS.epoch ~= 3), 'bo');
hold on;
swarmchart(CC_OPTICS.altId(string(CC_OPTICS.opticsCue)=='z' & ...
    CC_OPTICS.epoch ~=3), ...
    CC_OPTICS.snr_dbFiltLed2(string(CC_OPTICS.opticsCue)=='z' & ...
    CC_OPTICS.epoch ~= 3),'ro');
hold off;
xlabel('Animal ID');
ylabel('SNR, dB');
title('Optical channel 2 · 1050 nm');
grid;
legend('ccTrue','ccFalse');



figure;
swarmchart(CC_OPTICS.altId(string(CC_OPTICS.opticsCue)=='3'& ...
    CC_OPTICS.epoch ~= 3), ...
    CC_OPTICS.snr_dbFiltLed3( string(CC_OPTICS.opticsCue)=='3' & ...
    CC_OPTICS.epoch ~= 3), 'bo');
hold on;
swarmchart(CC_OPTICS.altId(string(CC_OPTICS.opticsCue)=='z' & ...
    CC_OPTICS.epoch ~=3), ...
    CC_OPTICS.snr_dbFiltLed3(string(CC_OPTICS.opticsCue)=='z' & ...
    CC_OPTICS.epoch ~= 3),'ro');
hold off;
xlabel('Animal ID');
ylabel('SNR, dB');
title('Optical channel 3 · 1200 nm');
grid;
legend('ccTrue','ccFalse');

%% Swarmchart for all animals, epoch 3 only


figure;
swarmchart(CC_OPTICS.altId(string(CC_OPTICS.opticsCue)=='3'& ...
    CC_OPTICS.epoch == 3), ...
    CC_OPTICS.snr_dbFiltLed2( string(CC_OPTICS.opticsCue)=='3' & ...
    CC_OPTICS.epoch == 3), 'bo');
hold on;
swarmchart(CC_OPTICS.altId(string(CC_OPTICS.opticsCue)=='z' & ...
    CC_OPTICS.epoch ==3), ...
    CC_OPTICS.snr_dbFiltLed2(string(CC_OPTICS.opticsCue)=='z' & ...
    CC_OPTICS.epoch == 3),'ro');
hold off;
xlabel('Animal ID');
ylabel('SNR, dB');
title('Optical channel 2 · 1050 nm');
grid;
legend('ccTrue','ccFalse');



figure;
swarmchart(CC_OPTICS.altId(string(CC_OPTICS.opticsCue)=='3'& ...
    CC_OPTICS.epoch == 3), ...
    CC_OPTICS.snr_dbFiltLed3( string(CC_OPTICS.opticsCue)=='3' & ...
    CC_OPTICS.epoch == 3), 'bo');
hold on;
swarmchart(CC_OPTICS.altId(string(CC_OPTICS.opticsCue)=='z' & ...
    CC_OPTICS.epoch ==3), ...
    CC_OPTICS.snr_dbFiltLed3(string(CC_OPTICS.opticsCue)=='z' & ...
    CC_OPTICS.epoch == 3),'ro');
hold off;
xlabel('Animal ID');
ylabel('SNR, dB');
title('Optical channel 3 · 1200 nm');
grid;
legend('ccTrue','ccFalse');

%% Swarmcharts for SNR LED2 and SNR LED3, all epochs, all animals

figure('Color', white);

s1 = subplot(121);

swarmchart(CC_OPTICS.altId(string(CC_OPTICS.opticsCue)=='z' & ...
    CC_OPTICS.epoch ==3), ...
    CC_OPTICS.snr_dbFiltLed2(string(CC_OPTICS.opticsCue)=='z' & ...
    CC_OPTICS.epoch == 3), 'x', 'MarkerEdgeColor', goldenrod);

hold on;

swarmchart(CC_OPTICS.altId(string(CC_OPTICS.opticsCue)=='3'& ...
    CC_OPTICS.epoch == 3), ...
    CC_OPTICS.snr_dbFiltLed2( string(CC_OPTICS.opticsCue)=='3' & ...
    CC_OPTICS.epoch == 3), 'd', 'MarkerEdgeColor', green);

swarmchart(CC_OPTICS.altId(string(CC_OPTICS.opticsCue)=='z' & ...
    CC_OPTICS.epoch ~=3), ...
    CC_OPTICS.snr_dbFiltLed2(string(CC_OPTICS.opticsCue)=='z' & ...
    CC_OPTICS.epoch ~= 3), 'Marker', '+', 'MarkerEdgeColor', red);

swarmchart(CC_OPTICS.altId(string(CC_OPTICS.opticsCue)=='3'& ...
    CC_OPTICS.epoch ~= 3), ...
    CC_OPTICS.snr_dbFiltLed2( string(CC_OPTICS.opticsCue)=='3' & ...
    CC_OPTICS.epoch ~= 3), 'Marker', 'o', 'MarkerEdgeColor', blue);

hold off;
xlabel('Animal ID', 'FontSize', 15);
xticks([1:6]);
xticklabels({'63HF','01L5','6JK5','83H1','9ON6','9FL3'});
ylabel('SNR, dB', 'FontSize', 15);
ylim([-40 80]);
title('Optical channel 2 · 1050 nm');
grid;
legend('ccMissed E3', 'ccMatched E3', 'ccMissed E1 & E5', 'ccMatched E1 & E5', ...
    'Location', 'southeast');


s2 = subplot(122);


swarmchart(CC_OPTICS.altId(string(CC_OPTICS.opticsCue)=='z' & ...
    CC_OPTICS.epoch ==3), ...
    CC_OPTICS.snr_dbFiltLed3(string(CC_OPTICS.opticsCue)=='z' & ...
    CC_OPTICS.epoch == 3), 'Marker', 'x', 'MarkerEdgeColor', goldenrod);

hold on;

swarmchart(CC_OPTICS.altId(string(CC_OPTICS.opticsCue)=='3'& ...
    CC_OPTICS.epoch == 3), ...
    CC_OPTICS.snr_dbFiltLed3( string(CC_OPTICS.opticsCue)=='3' & ...
    CC_OPTICS.epoch == 3), 'd', 'MarkerEdgeColor', green );

swarmchart(CC_OPTICS.altId(string(CC_OPTICS.opticsCue)=='z' & ...
    CC_OPTICS.epoch ~=3), ...
    CC_OPTICS.snr_dbFiltLed3(string(CC_OPTICS.opticsCue)=='z' & ...
    CC_OPTICS.epoch ~= 3), 'Marker', '+', 'MarkerEdgeColor', red);

swarmchart(CC_OPTICS.altId(string(CC_OPTICS.opticsCue)=='3'& ...
    CC_OPTICS.epoch ~= 3), ...
    CC_OPTICS.snr_dbFiltLed3( string(CC_OPTICS.opticsCue)=='3' & ...
    CC_OPTICS.epoch ~= 3), 'Marker', 'o', 'MarkerEdgeColor', blue);


hold off;
xlabel('Animal ID', 'FontSize', 15);
xticks([1:6]);
xticklabels({'63HF','01L5','6JK5','83H1','9ON6','9FL3'});
ylabel('SNR, dB', 'FontSize', 15);
ylim([-40 80]);
title('Optical channel 3 · 1200 nm');
grid;
legend('ccMissed E3', 'ccMatched E3', 'ccMissed E1 & E5', 'ccMatched E1 & E5', ...
    'Location', 'southeast');

%%

load('/Users/dave/Desktop/DTAG/raw/tt21_134braw.mat')


figure('Color', white);
subplot(131);
thd(E1_DECIMATED.filtLed2, 50);
xlim([0 20]);
ylim([-40 80]);
title('Optical channel 2 · 1050 nm');
subplot(132);
thd(E1_DECIMATED.filtLed3, 50);
xlim([0 20]);
ylim([-40 80]);
title('Optical channel 3 · 1200 nm');
subplot(133);
thd(E1_DECIMATED.filtLed1, 50);
xlim([0 20]);
ylim([-40 80]);
title('Optical channel 1 · ambient');

%% WSST and ridge finding...

% load('/Users/dave/Desktop/DTAG/raw/tt21_134araw.mat')
load('/Users/dave/Desktop/DTAG/raw/tt21_134braw.mat')
% load('/Users/dave/Desktop/DTAG/raw/tt21_141eraw.mat')
% load('/Users/dave/Desktop/DTAG/raw/tt21_142craw.mat')

fprintf('Working up sqrt(LED2^2 + LED3^2) wsst and ridges...\n')

[wsstE1, fE1] = wsst( sqrt( E1_DECIMATED.filtLed2.^2 + ...
    E1_DECIMATED.filtLed3.^2), 50);
[wsstE3, fE3] = wsst( sqrt( E3_DECIMATED.filtLed2.^2 + ...
    E3_DECIMATED.filtLed3.^2), 50);
[wsstE5, fE5] = wsst( sqrt( E5_DECIMATED.filtLed2.^2 + ...
    E5_DECIMATED.filtLed3.^2), 50);

useMeE1 = logical(fE1(:)> 0.4 & fE1(:)<5);
useMeE3 = logical(fE3(:)> 0.4 & fE3(:)<5);
useMeE5 = logical(fE5(:)> 0.4 & fE5(:)<5);

[ridgeE1, iridgeE1] = wsstridge(wsstE1(useMeE1,:), 60, ...
    fE1(useMeE1), 'NumRidges', 3);
[ridgeE3, iridgeE3] = wsstridge(wsstE3(useMeE3,:), 60, ...
    fE3(useMeE3), 'NumRidges', 3);
[ridgeE5, iridgeE5] = wsstridge(wsstE5(useMeE5,:), 60, ...
    fE5(useMeE5), 'NumRidges', 3);

figure('Color', white);
pcolor(E1_DECIMATED.Time, fE1, abs(wsstE1) ); 
shading interp; 
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Optical channel 2+3');
ylim([0 10]);

figure('Color', white);
pcolor(E3_DECIMATED.Time, fE3, abs(wsstE3) ); 
shading interp; 
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Optical channel 2+3');
ylim([0 10]);

figure('Color', white);
pcolor(E5_DECIMATED.Time, fE5, abs(wsstE5) );
shading interp; 
xlabel('Time, local');
ylabel('Frequency, Hz');
ylim([0 10]);
titleTxt = 'Optical channel 2 · 1050 nm · tt21_134b · E5';
title(titleTxt, 'Interpreter', 'none');


figure('Color', white);
plot(E1_DECIMATED.Time, ridgeE1);
hold on;
plot(E1_DECIMATED.Time, E1_DECIMATED.wsstHr, 'Color', maroon);
hold off;
grid;
legend;

figure('Color', white);
plot(E3_DECIMATED.Time, ridgeE3);
hold on;
plot(E3_DECIMATED.Time, E3_DECIMATED.wsstHr, 'Color', maroon);
hold off;
grid;
legend;

figure('Color', white);
plot(E5_DECIMATED.Time, ridgeE5);
hold on;
plot(E5_DECIMATED.Time, E5_DECIMATED.wsstHr, 'Color', maroon);
hold off;
grid;
legend;




%%

aoTimes = CC_DATA.Time( string(CC_DATA.ID)=='6JK5' & CC_DATA.epoch == 5 );

ventTimes = E5_CC.Time(string(E5_CC.cueType)=='v');

figure('Color', white);

s4a = subplot(311);
p4a = plot(E5_DECIMATED.Time, E5_DECIMATED.filtLed2, 'Color', blue);
hold on;
yPos = s4a.YLim(2) - ((s4a.YLim(2) - s4a.YLim(1)) / 10);
p4aAO = xline(aoTimes, '--', 'Color', goldenrod, 'LineWidth', 1.2); 
p4aV = plot(ventTimes, yPos, 'v-', 'Color', green, 'LineWidth', 1.2, ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', green );
p4aVlines = xline(ventTimes, '-', 'Color', green, 'LineWidth', 1.2);
hold off;
xlabel('Time, local');
ylabel('Reflectance, A.U.');
xlim([E1_DECIMATED.Time(1), E5_DECIMATED.Time(end)]);
title('Optical channel 2 · 1050 nm');
grid;
lgndA = legend([ p4a p4aAO(1) p4aV(1)], ...
    '1050 nm','AO', 'vents', 'Location', 'southeast');

s4b = subplot(312);
% plot(E1_DECIMATED.Time, smooth(E1_DECIMATED.filtLed1), ...
%     'Color', black);
p4b = plot(E5_DECIMATED.Time, smooth(E5_DECIMATED.filtLed3), ...
    'Color', red);
hold on;
xline(aoTimes, '--', 'Color', goldenrod, 'LineWidth', 1.2);
yPos = s4b.YLim(2) - ((s4b.YLim(2) - s4b.YLim(1)) / 10);
p4bAO = xline(aoTimes, '--', 'Color', goldenrod, 'LineWidth', 1.2); 
p4bV = plot(ventTimes, yPos, 'v-', 'Color', green, 'LineWidth', 1.2, ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', green );
p4bVlines = xline(ventTimes, '-', 'Color', green, 'LineWidth', 1.2);
hold off;
xlabel('Time, local');
ylabel('Reflectance, A.U.');
title('Optical channel 1 · ambient');
grid;
lgndB = legend([ p4b p4bAO(1) p4bV(1)], ...
    '1200 nm','AO', 'vents', 'Location', 'southeast');


colororder([0.8500 0.3250 0.0980
    0.4660 0.6740 0.1880]);
s4c = subplot(313);
yyaxis left;
plot(E5_DECIMATED.Time, E5_DECIMATED.ax);
ylabel('bpAx, m·s^-^2');
yyaxis right;
plot(E5_DECIMATED.Time, E5_DECIMATED.gy);
hold on;
xline(aoTimes, '--', 'Color', goldenrod, 'LineWidth', 1.2);
xlabel('Time, local');
ylabel('bpGy, degrees·s^-^1');
title('SCG and GCG · Band-passed x-axis accelerometer and y-axis gyroscope');
grid;
legend('bpAx','bpGy', 'AO', 'Location', 'southeast');

linkaxes([s4a s4b s4c],'x');


%%

figure('Color', white);

s4a = subplot(411);
plot(E1_DECIMATED.Time, E1_DECIMATED.filtLed2, 'Color', blue);
hold on;
yPos = s4a.YLim(2) - ((s4a.YLim(2) - s4a.YLim(1)) / 10);
ftPlotCues(E1_CC.Time(string(E1_CC.cueType)=='v'), ...
    E1_CC.cueType(string(E1_CC.cueType)=='v'), ...
    yPos);
hold off;
xlabel('Time, local');
ylabel('Reflectance, A.U.');
xlim([E1_DECIMATED.Time(1), E1_DECIMATED.Time(end)]);
ylim([-4000 4000]);
title('Optical channel 2 · 1050 nm');
grid;
legend('1050 nm','vents', 'Location', 'southeast');

s4b = subplot(412);
plot(E1_DECIMATED.Time, E1_DECIMATED.filtLed3, 'Color', maroon);
hold on;
yPos = s4b.YLim(2) - ((s4b.YLim(2) - s4b.YLim(1)) / 10);
ftPlotCues(E1_CC.Time(string(E1_CC.cueType)=='v'), ...
    E1_CC.cueType(string(E1_CC.cueType)=='v'), ...
    yPos);
hold off;
xlabel('Time, local');
ylabel('Reflectance, A.U.');
ylim([-4000 4000]);
title('Optical channel 3 · 1200 nm');
grid;
legend('1200 nm','vents', 'Location', 'southeast');

s4c = subplot(413);
plot(E1_DECIMATED.Time, E1_DECIMATED.filtLed1, 'Color', black);
xlabel('Time, local');
ylabel('Reflectance, A.U.');
title('Optical channel 1 · ambient');
grid;


colororder([0.8500 0.3250 0.0980
    0.4660 0.6740 0.1880]);
s4d = subplot(414);
yyaxis left;
plot(E1_DECIMATED.Time, E1_DECIMATED.ax);
ylabel('SCG, m·s^-^2');
yyaxis right;
plot(E1_DECIMATED.Time, E1_DECIMATED.gy);
xlabel('Time, local', 'FontSize', 15);
ylabel('GCG, degrees · s^-^1');
title('SCG and GCG');
grid;
legend('SCG','GCG', 'Location', 'southeast');



linkaxes([s4a s4b s4c s4d],'x');

%%

figure('Color', white);

s4a = subplot(411);
plot(E5_DECIMATED.Time, E5_DECIMATED.filtLed2, 'Color', blue);
hold on;
yPos = s4a.YLim(2) - ((s4a.YLim(2) - s4a.YLim(1)) / 10);
ftPlotCues(E5_CC.Time(string(E5_CC.cueType)=='v'), ...
    E5_CC.cueType(string(E5_CC.cueType)=='v'), ...
    yPos);
hold off;
xlabel('Time, local');
ylabel('Reflectance, A.U.');
xlim([E5_DECIMATED.Time(1), E5_DECIMATED.Time(end)]);
ylim([-2000 2000]);
title('Optical channel 2 · 1050 nm');
grid;
legend('1050 nm','vents', 'Location', 'southeast');

s4b = subplot(412);
plot(E5_DECIMATED.Time, E5_DECIMATED.filtLed3, 'Color', maroon);
hold on;
yPos = s4b.YLim(2) - ((s4b.YLim(2) - s4b.YLim(1)) / 10);
ftPlotCues(E5_CC.Time(string(E5_CC.cueType)=='v'), ...
    E5_CC.cueType(string(E5_CC.cueType)=='v'), ...
    yPos);
hold off;
xlabel('Time, local');
ylabel('Reflectance, A.U.');
ylim([-2000 2000]);
title('Optical channel 3 · 1200 nm');
grid;
legend('1200 nm','vents', 'Location', 'southeast');

s4c = subplot(413);
plot(E5_DECIMATED.Time, E5_DECIMATED.filtLed1, 'Color', black);
xlabel('Time, local');
ylabel('Tag depth, meters');
set(gca,'YDir','reverse');
grid;


colororder([0.8500 0.3250 0.0980
    0.4660 0.6740 0.1880]);
s4d = subplot(414);
yyaxis left;
plot(E5_DECIMATED.Time, E5_DECIMATED.ax);
ylabel('bpAx, m·s^-^2');
yyaxis right;
plot(E5_CC_ENERGY.Time, E5_CC_ENERGY.iTotalKE, 'o');
ylabel('iKE, J · s^-^1');
title('SCG and iKE');
grid;
legend('bpAx','iKE', 'Location', 'southeast');



linkaxes([s4a s4b s4c],'x');


%%

aoTimes = E3_CC.Time(string(E3_CC.cueType)=='ao');

figure('Color', white);

s9a = subplot(511);
plot(E3_DECIMATED.Time, E3_DECIMATED.filtLed2, 'Color', blue);
hold on;
xline(aoTimes, '--', 'Color', maroon)
hold off;
xlabel('Time, local');
ylabel('A.U.');
title('Optical channel 2 · 1050 nm');
grid;
legend('1050 nm','ao', 'Location', 'southeast');

s9b = subplot(512);
plot(E3_DECIMATED.Time, E3_DECIMATED.filtLed3, 'Color', red);
hold on;
xline(aoTimes, '--', 'Color', maroon)
hold off;
xlabel('Time, local');
ylabel('A.U.');
title('Optical channel 3 · 1200 nm');
grid;
legend('1200 nm','ao', 'Location', 'southeast');

s9c = subplot(513);
plot(E3_DECIMATED.Time, E3_DECIMATED.filtLed4, 'Color', goldenrod);
hold on;
xline(aoTimes, '--', 'Color', maroon)
hold off;
xlabel('Time, local');
ylabel('A.U.');
title('Optical channel 4 · inoperative 950 nm');
grid;
legend('non-op 950 nm','ao', 'Location', 'southeast');

s9d = subplot(514);
plot(E3_DECIMATED.Time, E3_DECIMATED.filtLed1, 'Color', black);
hold on;
xline(aoTimes, '--', 'Color', maroon)
hold off;
xlabel('Time, local');
ylabel('A.U.');
title('Optical channel 4 · ambient');
grid;
legend('ambient','ao', 'Location', 'southeast');

s9e = subplot(515);
plot(E3_DECIMATED.Time, E3_DECIMATED.ax, 'Color', maroon);
hold on;
xline(aoTimes, '-.', 'Color', maroon);
hold off;
xlabel('Time, local');
ylabel('m · s^-^2');
title('Band-passed x-axis accelerometer (SCG)');
grid;
legend('bpAx','ao', 'Location', 'southeast');

linkaxes([s9a s9b s9c s9d s9e],'x');

%%

%% Mann-Whitney U test for unpaired samples

% this compares iPower (dbFiltLedX)

statsE1E5 = mwwtest( (CC_OPTICS.dbFiltLed2(CC_HR_DATA.epoch == 1 & ...
(CC_HR_DATA.altId >= 2 & CC_HR_DATA.altId <= 5)))' , ...
(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & ...
(CC_HR_DATA.altId >= 2 & CC_HR_DATA.altId <= 5)))')

statsE3E5 = mwwtest( (CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 3 & ...
(CC_HR_DATA.altId >= 2 & CC_HR_DATA.altId <= 5)))' , ...
(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & ...
(CC_HR_DATA.altId >= 2 & CC_HR_DATA.altId <= 5)))')

statsE1E3 = mwwtest( (CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & ...
(CC_HR_DATA.altId >= 2 & CC_HR_DATA.altId <= 5)))' , ...
(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 3 & ...
(CC_HR_DATA.altId >= 2 & CC_HR_DATA.altId <= 5)))')

effectE1E5 = statsE1E5.Z / sqrt(sum(statsE1E5.n))
effectE3E5 = statsE3E5.Z / sqrt(sum(statsE3E5.n))
effectE1E3 = statsE1E3.Z / sqrt(sum(statsE1E3.n))


%%

figure('Color', white); 

s10a = subplot(434);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 1), ...
    CC_OPTICS.dbFiltLed2(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 1), ...
    'Color', blue);
xlabel('Time, local');
ylabel('SNR, dB');
title('Optical channel 2 · 1050 nm · Epoch 1');
s10b = subplot(435);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 3), ...
    CC_OPTICS.dbFiltLed2(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 3), ...
    'Color', blue);
xlabel('Time, local');
ylabel('SNR, dB');
title('Optical channel 2 · 1050 nm · Epoch 3');
s10c = subplot(436);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 5), ...
    CC_OPTICS.dbFiltLed2(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 5), ...
    'Color', blue);
xlabel('Time, local');
ylabel('SNR, dB');
title('Optical channel 2 · 1050 nm · Epoch 5');
grid;

s10d = subplot(437);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 1), ...
    CC_OPTICS.dbFiltLed3(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 1), ...
    'Color', red);
xlabel('Time, local');
ylabel('SNR, dB');
title('Optical channel 3 · 1200 nm · Epoch 1');
grid;
s10e = subplot(438);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 3), ...
    CC_OPTICS.dbFiltLed3(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 3), ...
    'Color', red);
xlabel('Time, local');
ylabel('SNR, dB');
title('Optical channel 3 · 1200 nm · Epoch 3');
grid;
s10f = subplot(439);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 5), ...
    CC_OPTICS.dbFiltLed3(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 5), ...
    'Color', red);
xlabel('Time, local');
ylabel('SNR, dB');
title('Optical channel 3 · 1200 nm · Epoch 5');
grid;

s10f = subplot(4,3,10);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 1), ...
    CC_OPTICS.dbFiltLed1(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 1), ...
    'Color', black);
xlabel('Time, local');
ylabel('SNR, dB');
title('Ambient optical channel · Epoch 1');
grid;
s10g = subplot(4,3,11);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 3), ...
    CC_OPTICS.dbFiltLed1(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 3), ...
    'Color', black);
xlabel('Time, local');
ylabel('SNR, dB');
title('Ambient optical channel · Epoch 3');
grid;
s10h = subplot(4,3,12);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 5), ...
    CC_OPTICS.dbFiltLed1(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 5), ...
    'Color', black);
xlabel('Time, local');
ylabel('SNR, dB');
title('Ambient optical channel · Epoch 5');
grid;

s10i = subplot(4,3,1);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 1), ...
    CC_OPTICS.wsstHr(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 1), ...
    'Color', purple);
xlabel('Time, local');
ylabel('beats · min^-^1');
title('Ensemble KCG HR · Epoch 1');
grid;
s10j = subplot(4,3,2);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 3), ...
    CC_OPTICS.wsstHr(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 3), ...
    'Color', purple);
xlabel('Time, local');
ylabel('beats · min^-^1');
title('Ensemble KCG HR · Epoch 3');
grid;
s10j = subplot(4,3,3);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 5), ...
    CC_OPTICS.wsstHr(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 5), ...
    'Color', purple);
xlabel('Time, local');
ylabel('beats · min^-^1');
title('Ensemble KCG HR · Epoch 5');
grid;

%%

figure('Color', white); 

s10a = subplot(434);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 1), ...
    CC_OPTICS.dbFiltLed2(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 1), ...
    'Color', blue);
xlabel('Time, local');
ylabel('SNR, dB');
title('Optical channel 2 · 1050 nm · Epoch 1');
s10b = subplot(435);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 3), ...
    CC_OPTICS.dbFiltLed2(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 3), ...
    'Color', blue);
xlabel('Time, local');
ylabel('SNR, dB');
title('Optical channel 2 · 1050 nm · Epoch 3');
s10c = subplot(436);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 5), ...
    CC_OPTICS.dbFiltLed2(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 5), ...
    'Color', blue);
xlabel('Time, local');
ylabel('SNR, dB');
title('Optical channel 2 · 1050 nm · Epoch 5');
grid;

s10d = subplot(437);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 1), ...
    CC_OPTICS.dbFiltLed3(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 1), ...
    'Color', red);
xlabel('Time, local');
ylabel('SNR, dB');
title('Optical channel 3 · 1200 nm · Epoch 1');
grid;
s10e = subplot(438);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 3), ...
    CC_OPTICS.dbFiltLed3(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 3), ...
    'Color', red);
xlabel('Time, local');
ylabel('SNR, dB');
title('Optical channel 3 · 1200 nm · Epoch 3');
grid;
s10f = subplot(439);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 5), ...
    CC_OPTICS.dbFiltLed3(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 5), ...
    'Color', red);
xlabel('Time, local');
ylabel('SNR, dB');
title('Optical channel 3 · 1200 nm · Epoch 5');
grid;

s10f = subplot(4,3,10);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 1), ...
    CC_OPTICS.dbFiltLed1(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 1), ...
    'Color', black);
xlabel('Time, local');
ylabel('SNR, dB');
title('Ambient optical channel · Epoch 1');
grid;
s10g = subplot(4,3,11);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 3), ...
    CC_OPTICS.dbFiltLed1(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 3), ...
    'Color', black);
xlabel('Time, local');
ylabel('SNR, dB');
title('Ambient optical channel · Epoch 3');
grid;
s10h = subplot(4,3,12);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 5), ...
    CC_OPTICS.dbFiltLed1(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 5), ...
    'Color', black);
xlabel('Time, local');
ylabel('SNR, dB');
title('Ambient optical channel · Epoch 5');
grid;

s10i = subplot(4,3,1);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 1), ...
    hampel( CC_OPTICS.iTotalKE(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 1), 3, 3), ...
    'Color', goldenrod);
xlabel('Time, local');
ylabel('beats · min^-^1');
title('iKE · Epoch 1');
grid;
s10j = subplot(4,3,2);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 3), ...
    CC_OPTICS.iTotalKE(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 3), ...
    'Color', goldenrod);
xlabel('Time, local');
ylabel('beats · min^-^1');
title('iKE · Epoch 3');
grid;
s10j = subplot(4,3,3);
plot(CC_OPTICS.Time(CC_OPTICS.altId==3 & CC_OPTICS.epoch == 5), ...
    hampel( CC_OPTICS.iTotalKE(CC_OPTICS.altId==3  & CC_OPTICS.epoch == 5), 3, 3), ...
    'Color', goldenrod);
xlabel('Time, local');
ylabel('beats · min^-^1');
title('iKE · Epoch 5');
grid;


%%

figure('Color', white); 

subplot(231);
scatter( (CC_OPTICS.iTotalKE(CC_OPTICS.epoch == 1)), ...
    CC_OPTICS.snr_dbFiltLed2(CC_OPTICS.epoch == 1), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', green);
hold on;
plot( (CC_OPTICS.iTotalKE(CC_OPTICS.epoch == 1)), ...
    (CC_OPTICS.iTotalKE(CC_OPTICS.epoch == 1) .* 7.71), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', green);
hold off;
xlabel('log_1_0(iKE), J · s^-^1');
ylabel('SNR, db');
title('Optical channel 2 · 1050 nm · Epoch 1');
grid;

subplot(232);
scatter( log10(CC_OPTICS.iTotalKE(CC_OPTICS.epoch == 3)), ...
    CC_OPTICS.snr_dbFiltLed2(CC_OPTICS.epoch == 3), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', blue);
xlabel('log_1_0(iKE), J · s^-^1');
ylabel('SNR, db');
title('Optical channel 2 · 1050 nm · Epoch 3');
grid;

subplot(233);
scatter( (CC_OPTICS.iTotalKE(CC_OPTICS.epoch == 5)), ...
    CC_OPTICS.snr_dbFiltLed2(CC_OPTICS.epoch == 5), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', red);
% iKE: \beta1: 28.43, CI95 = 19.04,37.83; 
%      \beta0 = 243.85, CI95 = 180.68,307.02
hold on;
plot((CC_OPTICS.iTotalKE(CC_OPTICS.epoch == 5)), ...
    (CC_OPTICS.iTotalKE(CC_OPTICS.epoch == 5) * 28.43), ...
    'Color', black);
hold off;
xlabel('log_1_0(iKE), J · s^-^1');
ylabel('SNR, db');
title('Optical channel 2 · 1050 nm · Epoch 5');
grid;


subplot(234);
scatter( (CC_OPTICS.iHR(CC_OPTICS.epoch == 1)), ...
    CC_OPTICS.snr_dbFiltLed2(CC_OPTICS.epoch == 1), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', green);
xlabel('ifH, beats · m^-^1');
ylabel('SNR, db');
title('Optical channel 2 · 1050 nm · Epoch 1');

subplot(235);
scatter( (CC_OPTICS.iHR(CC_OPTICS.epoch == 3)), ...
    CC_OPTICS.snr_dbFiltLed2(CC_OPTICS.epoch == 3), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', blue);
xlabel('ifH, beats · m^-^1');
ylabel('SNR, db');
title('Optical channel 2 · 1050 nm · Epoch 3');

subplot(236);
scatter( (CC_OPTICS.iHR(CC_OPTICS.epoch == 5)), ...
    CC_OPTICS.snr_dbFiltLed2(CC_OPTICS.epoch == 5), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', red);
xlabel('ifH, beats · m^-^1');
ylabel('SNR, db');
title('Optical channel 2 · 1050 nm · Epoch 5');
grid;

%%

figure('Color', white); 

subplot(231);
scatter( (CC_OPTICS.iTotalKE(CC_OPTICS.epoch == 1)), ...
    CC_OPTICS.snr_dbFiltLed3(CC_OPTICS.epoch == 1), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', green);
hold on;
plot( (CC_OPTICS.iTotalKE(CC_OPTICS.epoch == 1)), ...
    (CC_OPTICS.iTotalKE(CC_OPTICS.epoch == 1) .* 10.00 + -121.97), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', green);
hold off;
xlabel('log_1_0(iKE), J · s^-^1');
ylabel('SNR, db');
title('Optical channel 3 · 1200 nm · Epoch 1');
grid;

subplot(232);
scatter( log10(CC_OPTICS.iTotalKE(CC_OPTICS.epoch == 3)), ...
    CC_OPTICS.snr_dbFiltLed3(CC_OPTICS.epoch == 3), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', blue);
xlabel('log_1_0(iKE), J · s^-^1');
ylabel('SNR, db');
title('Optical channel 3 · 1200 nm · Epoch 3');
grid;

subplot(233);
scatter( log10(CC_OPTICS.iTotalKE(CC_OPTICS.epoch == 5)), ...
    CC_OPTICS.snr_dbFiltLed3(CC_OPTICS.epoch == 5), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', red);
xlabel('log_1_0(iKE), J · s^-^1');
ylabel('SNR, db');
title('Optical channel 3 · 1200 nm · Epoch 5');
grid;


subplot(234);
scatter( (CC_OPTICS.iHR(CC_OPTICS.epoch == 1)), ...
    CC_OPTICS.snr_dbFiltLed3(CC_OPTICS.epoch == 1), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', green);
xlabel('ifH, beats · m^-^1');
ylabel('SNR, db');
title('Optical channel 3 · 1200 nm · Epoch 1');

subplot(235);
scatter( (CC_OPTICS.iHR(CC_OPTICS.epoch == 3)), ...
    CC_OPTICS.snr_dbFiltLed3(CC_OPTICS.epoch == 3), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', blue);
xlabel('ifH, beats · m^-^1');
ylabel('SNR, db');
title('Optical channel 3 · 1200 nm · Epoch 3');

subplot(236);
scatter( (CC_OPTICS.iHR(CC_OPTICS.epoch == 5)), ...
    CC_OPTICS.snr_dbFiltLed3(CC_OPTICS.epoch == 5), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', red);
xlabel('ifH, beats · m^-^1');
ylabel('SNR, db');
title('Optical channel 3 · 1200 nm · Epoch 5');
grid;

%%

figure('Color', white); 

subplot(231);
scatter( (CC_OPTICS.mass(CC_OPTICS.epoch == 1)), ...
    CC_OPTICS.snr_dbFiltLed2(CC_OPTICS.epoch == 1), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', green);
xlabel('mass, kg');
ylabel('SNR, db');
title('Optical channel 2 · 1050 nm · Epoch 1');
grid;

subplot(232);
scatter( (CC_OPTICS.mass(CC_OPTICS.epoch == 3)), ...
    CC_OPTICS.snr_dbFiltLed2(CC_OPTICS.epoch == 3), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', blue);
xlabel('mass, kg');
ylabel('SNR, db');
title('Optical channel 2 · 1050 nm · Epoch 3');
grid;

subplot(233);
scatter( (CC_OPTICS.mass(CC_OPTICS.epoch == 5)), ...
    CC_OPTICS.snr_dbFiltLed2(CC_OPTICS.epoch == 5), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', red);
xlabel('mass, kg');
ylabel('SNR, db');
title('Optical channel 2 · 1050 nm · Epoch 5');
grid;


subplot(234);
scatter( (CC_OPTICS.temp(CC_OPTICS.epoch == 1)), ...
    CC_OPTICS.snr_dbFiltLed2(CC_OPTICS.epoch == 1), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', green);
xlabel('Temperature, °C');
ylabel('SNR, db');
title('Optical channel 2 · 1050 nm · Epoch 1');

subplot(235);
scatter( (CC_OPTICS.temp(CC_OPTICS.epoch == 3)), ...
    CC_OPTICS.snr_dbFiltLed2(CC_OPTICS.epoch == 3), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', blue);
xlabel('Temperature, °C');
ylabel('SNR, db');
title('Optical channel 2 · 1050 nm · Epoch 3');

subplot(236);
scatter( (CC_OPTICS.temp(CC_OPTICS.epoch == 5)), ...
    CC_OPTICS.snr_dbFiltLed2(CC_OPTICS.epoch == 5), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', red);
xlabel('Temperature, °C');
ylabel('SNR, db');
title('Optical channel 2 · 1050 nm · Epoch 5');
grid;


%%

figure('Color', white); 

subplot(231);
scatter( (CC_OPTICS.mass(CC_OPTICS.epoch == 1)), ...
    CC_OPTICS.snr_dbFiltLed3(CC_OPTICS.epoch == 1), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', green);
xlabel('mass, kg');
ylabel('SNR, db');
title('Optical channel 3 · 1200 nm · Epoch 1');
grid;

subplot(232);
scatter( (CC_OPTICS.mass(CC_OPTICS.epoch == 3)), ...
    CC_OPTICS.snr_dbFiltLed3(CC_OPTICS.epoch == 3), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', blue);
xlabel('mass, kg');
ylabel('SNR, db');
title('Optical channel 3 · 1200 nm · Epoch 3');
grid;

subplot(233);
scatter( (CC_OPTICS.mass(CC_OPTICS.epoch == 5)), ...
    CC_OPTICS.snr_dbFiltLed3(CC_OPTICS.epoch == 5), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', red);
xlabel('mass, kg');
ylabel('SNR, db');
title('Optical channel 3 · 1200 nm · Epoch 5');
grid;


subplot(234);
scatter( (CC_OPTICS.temp(CC_OPTICS.epoch == 1)), ...
    CC_OPTICS.snr_dbFiltLed3(CC_OPTICS.epoch == 1), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', green);
xlabel('Temperature, °C');
ylabel('SNR, db');
title('Optical channel 3 · 1200 nm · Epoch 1');

subplot(235);
scatter( (CC_OPTICS.temp(CC_OPTICS.epoch == 3)), ...
    CC_OPTICS.snr_dbFiltLed3(CC_OPTICS.epoch == 3), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', blue);
xlabel('Temperature, °C');
ylabel('SNR, db');
title('Optical channel 3 · 1200 nm · Epoch 3');

subplot(236);
scatter( (CC_OPTICS.temp(CC_OPTICS.epoch == 5)), ...
    CC_OPTICS.snr_dbFiltLed3(CC_OPTICS.epoch == 5), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', red);
xlabel('Temperature, °C');
ylabel('SNR, db');
title('Optical channel 3 · 1200 nm · Epoch 5');
grid;

%%

figure('Color', [1 1 1] );

subplot(3,2,[1:2]);
h4 = heatmap(PHASE_TWO_TRIALS, 'location', 'score');
h4.XDisplayLabels = {'Post-nuchal [1]','Anterior dorsal [2]', ...
    'Left dorsal ridge [3L]','Chest [4]','Right dorsal ridge [3R]'};
h4.YDisplayLabels = {'0.00-0.19','0.20-0.39','0.40-0.59',...
    '0.60-0.79', '0.80-1.00'};
h4.XLabel = 'Tag Location';
h4.YLabel = 'OCG Detection Rate';
title('OCG detection rate by tag attachment location');

subplot(325);
h1 = heatmap(PHASE_TWO_TRIALS, 'epoch', 'score','ColorVariable', ...
    'snr_dbFiltLed2');
h1.XDisplayLabels = {'Baseline','Apnea','Recovery'};
h1.YDisplayLabels = {'0.00-0.19','0.20-0.39','0.40-0.59',...
    '0.60-0.79', '0.80-1.00'};
h1.XLabel = 'Epoch';
h1.YLabel = 'OCG Detection Rate';
title('Mean SNR of optical channel 2, 1050 nm');

subplot(326);
h2 = heatmap(PHASE_TWO_TRIALS, 'epoch', 'score', 'ColorVariable', ...
    'snr_dbFiltLed3');
h2.XDisplayLabels = {'Baseline','Apnea','Recovery'};
h2.YDisplayLabels = {'0.00-0.19','0.20-0.39','0.40-0.59',...
    '0.60-0.79', '0.80-1.00'};
h2.XLabel = 'Epoch';
h2.YLabel = 'OCG Detection Rate';
title('Mean SNR of optical channel 3, 1200 nm');


subplot(323);
h3 = heatmap(PHASE_TWO_TRIALS, 'epoch', 'score');
h3.XDisplayLabels = {'Baseline','Apnea','Recovery'};
h3.YDisplayLabels = {'0.00-0.19','0.20-0.39','0.40-0.59',...
    '0.60-0.79', '0.80-1.00'};
h3.XLabel = 'Epoch';
h3.YLabel = 'OCG Detection Rate';
title('OCG detection rate by epoch');

subplot(324);
h5 = heatmap(PHASE_TWO_TRIALS, 'altId', 'score');
h5.XDisplayLabels = {'63HF','01L5','6JK5','83H1','9ON6','9FL3'};
h5.YDisplayLabels = {'0.00-0.19','0.20-0.39','0.40-0.59',...
    '0.60-0.79', '0.80-1.00'};
h5.XLabel = 'Animal ID';
h5.YLabel = 'OCG Detection Rate';
title('OCG detection rate by animal');

writetable(PHASE_TWO_TRIALS, ...
    '/Users/dave/Documents/MATLAB/phaseTwoData.csv');


%%

P2 = PHASE_TWO_TRIALS;

goodRows = ~isnan(P2.score);

mwwtest(P2.score(P2.epoch == 1)', P2.score(P2.epoch(goodRows) == 5)' )
