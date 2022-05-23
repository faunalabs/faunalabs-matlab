%% ftChapterTwoStats.m
%
%   Dave Haas
%   24 September 2021
%   - this is a single file of all the statistical tests run for Chapter
%   Two of my dissertation. These were scattered all over the place, so
%   making this one file where everything that happens in Matlab can be
%   re-run in the future.

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

% screenSize = get(0,'screensize');
% 
% figAnalysisRegion = figure;
% figSize = figAnalysisRegion.Position;
% figSize(3) = screenSize(3);
% figSize(4) = screenSize(4) * 0.65;


%% Part One

% load the entire chest trial timetable series
load('/Users/dave/Documents/MATLAB/tursiopsTimetables.mat')


%% load animal timeseries

figSizeTrialPlots = [ 100 100 1200 400 ];

for trialNum = 1:6

    switch(trialNum)
        
        case 1
            load('/Users/dave/Desktop/DTAG/raw/tt21_128fraw.mat');
            
            figure('Color', white, 'Position', figSizeTrialPlots);

            plot( CC_HR_DATA.Time(CC_HR_DATA.altId == trialNum), ...
                CC_HR_DATA.iHr(CC_HR_DATA.altId == trialNum), 'o', 'Color', black);
            hold on;
            plot( CC_HR_DATA.Time(CC_HR_DATA.altId == trialNum), ...
                CC_HR_DATA.iHr(CC_HR_DATA.altId == trialNum), '--', 'Color', black);
            plot( E1_KINEMATICS.Time, E1_KINEMATICS.hrTimeConsensus, ...
                'LineWidth', 2', 'Color', green);
            plot( E3_KINEMATICS.Time, E3_KINEMATICS.hrTimeConsensus, ...
                'LineWidth', 2', 'Color', blue);
            xline(VENT_HR_DATA.Time(VENT_HR_DATA.altId == trialNum),  ...
                '--', 'Color', goldenrod, 'LineWidth', 1.5);
            
            hold off;
            xlabel('Time, local');
            ylabel('Heart rate, beats · min^-^1');
            titleTxt = sprintf('Heart rate profile, trial %s, animal ID: %s', ...
                TRIAL.tag, TRIAL.SUBJECT.id);
            title(titleTxt, 'Interpreter', 'none');
            grid;
            legend('iƒH','iƒH trendline', 'baseline wsstƒH', ...
                'submerged apnea wsstƒH', 'respiration', 'Location', 'best');
            
        case 2
            load('/Users/dave/Desktop/DTAG/raw/tt21_134araw.mat');
            
            figure('Color', white, 'Position', figSizeTrialPlots);

            plot( CC_HR_DATA.Time(CC_HR_DATA.altId == trialNum), ...
                CC_HR_DATA.iHr(CC_HR_DATA.altId == trialNum), 'o', 'Color', black);
            hold on;
            plot( CC_HR_DATA.Time(CC_HR_DATA.altId == trialNum), ...
                CC_HR_DATA.iHr(CC_HR_DATA.altId == trialNum), '--', 'Color', black);
            plot( E1_KINEMATICS.Time, E1_KINEMATICS.hrTimeConsensus, ...
                'LineWidth', 2', 'Color', green);
            plot( E3_KINEMATICS.Time, E3_KINEMATICS.hrTimeConsensus, ...
                'LineWidth', 2', 'Color', blue);
            plot( E5_KINEMATICS.Time, E5_KINEMATICS.hrTimeConsensus, ...
                'LineWidth', 2', 'Color', red);        
            xline(VENT_HR_DATA.Time(VENT_HR_DATA.altId == trialNum),  ...
                '--', 'Color', goldenrod, 'LineWidth', 1.5);
            hold off;
            xlabel('Time, local');
            ylabel('Heart rate, beats · min^-^1');
            titleTxt = sprintf('Heart rate profile, trial %s, animal ID: %s', ...
                TRIAL.tag, TRIAL.SUBJECT.id);
            title(titleTxt, 'Interpreter', 'none');
            grid;
            legend('iƒH','iƒH trendline', 'baseline wsstƒH', ...
                'submerged apnea wsstƒH', 'recovery wsstƒH', ...
                'respiration', 'Location', 'best');         
            
            
        case 3
            load('/Users/dave/Desktop/DTAG/raw/tt21_134braw.mat');
            
            figure('Color', white, 'Position', figSizeTrialPlots);

            plot( CC_HR_DATA.Time(CC_HR_DATA.altId == trialNum), ...
                CC_HR_DATA.iHr(CC_HR_DATA.altId == trialNum), 'o', 'Color', black);
            hold on;
            plot( CC_HR_DATA.Time(CC_HR_DATA.altId == trialNum), ...
                CC_HR_DATA.iHr(CC_HR_DATA.altId == trialNum), '--', 'Color', black);
            plot( E1_KINEMATICS.Time, E1_KINEMATICS.hrTimeConsensus, ...
                'LineWidth', 2', 'Color', green);
            plot( E3_KINEMATICS.Time, E3_KINEMATICS.hrTimeConsensus, ...
                'LineWidth', 2', 'Color', blue);
            plot( E5_KINEMATICS.Time, E5_KINEMATICS.hrTimeConsensus, ...
                'LineWidth', 2', 'Color', red);        
            xline(VENT_HR_DATA.Time(VENT_HR_DATA.altId == trialNum),  ...
                '--', 'Color', goldenrod, 'LineWidth', 2);
            hold off;
            xlabel('Time, local');
            ylabel('Heart rate, beats · min^-^1');
            titleTxt = sprintf('Heart rate profile, trial %s, animal ID: %s', ...
                TRIAL.tag, TRIAL.SUBJECT.id);
            title(titleTxt, 'Interpreter', 'none');
            grid;
            legend('iƒH','iƒH trendline','baseline wsstƒH', ...
                'submerged apnea wsstƒH', 'recovery wsstƒH', ...
                'respiration', 'Location', 'best');                  
            
            
        case 4
            load('/Users/dave/Desktop/DTAG/raw/tt21_141eraw.mat');
            
            figure('Color', white, 'Position', figSizeTrialPlots);

            plot( CC_HR_DATA.Time(CC_HR_DATA.altId == trialNum), ...
                CC_HR_DATA.iHr(CC_HR_DATA.altId == trialNum), 'o', 'Color', black);
            hold on;
            plot( CC_HR_DATA.Time(CC_HR_DATA.altId == trialNum), ...
                CC_HR_DATA.iHr(CC_HR_DATA.altId == trialNum), '--', 'Color', black);
            plot( E1_KINEMATICS.Time, E1_KINEMATICS.hrTimeConsensus, ...
                'LineWidth', 2', 'Color', green);
            plot( E3_KINEMATICS.Time, E3_KINEMATICS.hrTimeConsensus, ...
                'LineWidth', 2', 'Color', blue);
            plot( E5_KINEMATICS.Time, E5_KINEMATICS.hrTimeConsensus, ...
                    'LineWidth', 2', 'Color', red);        
            xline(VENT_HR_DATA.Time(VENT_HR_DATA.altId == trialNum),  ...
                '--', 'Color', goldenrod, 'LineWidth', 1.5);
            hold off;
            xlabel('Time, local');
            ylabel('Heart rate, beats · min^-^1');
            titleTxt = sprintf('Heart rate profile, trial %s, animal ID: %s', ...
                TRIAL.tag, TRIAL.SUBJECT.id);
            title(titleTxt, 'Interpreter', 'none');
            grid;
            legend('iƒH','iƒH trendline','baseline wsstƒH', ...
                'submerged apnea wsstƒH', 'recovery wsstƒH', ...
                'respiration', 'Location', 'best');                     
            
            
        case 5
            load('/Users/dave/Desktop/DTAG/raw/tt21_142craw.mat');
            
            figure('Color', white, 'Position', figSizeTrialPlots);

            plot( CC_HR_DATA.Time(CC_HR_DATA.altId == trialNum), ...
                CC_HR_DATA.iHr(CC_HR_DATA.altId == trialNum), 'o', 'Color', black);
            hold on;
            plot( CC_HR_DATA.Time(CC_HR_DATA.altId == trialNum), ...
                CC_HR_DATA.iHr(CC_HR_DATA.altId == trialNum), '--', 'Color', black);
            plot( E1_KINEMATICS.Time, E1_KINEMATICS.hrTimeConsensus, ...
                'LineWidth', 2', 'Color', green);
            plot( E3_KINEMATICS.Time, E3_KINEMATICS.hrTimeConsensus, ...
                'LineWidth', 2', 'Color', blue);
            plot( E5_KINEMATICS.Time, E5_KINEMATICS.hrTimeConsensus, ...
                'LineWidth', 2', 'Color', red);        
            xline(VENT_HR_DATA.Time(VENT_HR_DATA.altId == trialNum),  ...
                '--', 'Color', goldenrod, 'LineWidth', 1.5);
            hold off;
            xlabel('Time, local');
            ylabel('Heart rate, beats · min^-^1');
            titleTxt = sprintf('Heart rate profile, trial %s, animal ID: %s', ...
                TRIAL.tag, TRIAL.SUBJECT.id);
            title(titleTxt, 'Interpreter', 'none');
            grid;
            legend('iƒH','iƒH trendline','baseline wsstƒH', ...
                'submerged apnea wsstƒH', 'recovery wsstƒH', ...
                'respiration', 'Location', 'best');                    
            
            
        case 6
            load('/Users/dave/Desktop/DTAG/raw/tt21_142draw.mat');

            figure('Color', white, 'Position', figSizeTrialPlots);

            plot( CC_HR_DATA.Time(CC_HR_DATA.altId == trialNum), ...
                CC_HR_DATA.iHr(CC_HR_DATA.altId == trialNum), 'o', 'Color', black);
            hold on;
            plot( CC_HR_DATA.Time(CC_HR_DATA.altId == trialNum), ...
                CC_HR_DATA.iHr(CC_HR_DATA.altId == trialNum), '--', 'Color', black);
            plot( E1_KINEMATICS.Time, E1_KINEMATICS.hrTimeConsensus, ...
                'LineWidth', 2', 'Color', green);
            plot( E3_KINEMATICS.Time, E3_KINEMATICS.hrTimeConsensus, ...
                'LineWidth', 2', 'Color', blue);
            xline(VENT_HR_DATA.Time(VENT_HR_DATA.altId == trialNum),  ...
                '--', 'Color', goldenrod, 'LineWidth', 1.5);
            hold off;
            xlabel('Time, local');
            ylabel('Heart rate, beats · min^-^1');
            titleTxt = sprintf('Heart rate profile, trial %s, animal ID: %s', ...
                TRIAL.tag, TRIAL.SUBJECT.id);
            title(titleTxt, 'Interpreter', 'none');
            grid;
            legend('iƒH','iƒH trendline','baseline wsstƒH', ...
                'submerged apnea wsstƒH',  ...
                'respiration', 'Location', 'best');                  
            
    end
    

% generate a sample WSST + WSSTRIDGE plot for a baseline, apnea & recovery



end


%% Do Bland-Altman and correlation plots

% linear regression and Bland-Altman comparison of iHr 
% with linWsstHr, rotWsstHr, and ensWsstHr

figSizeBA = [ 100 100 1200 1200 ];

fig1 = figure('Position', figSizeBA, 'Color', white, 'NumberTitle', 'off');
ah1 = subplot(311);
ah2 = subplot(312);
ah3 = subplot(313);

[ rpc1, ah1, stats1 ] = BlandAltman(ah1, ...
    CC_HR_DATA.iHr, CC_HR_DATA.linWsstHr, ...
    {'ifH','linWsstH'}, ...
    'Comparison of linear WSST and instantaneous heart rate estimates (beats · min^-^1)', ...
    {'ifH : linWsstH','LoBF'}, ...
    'corrInfo', {'n','eq','r','r2','p'}, ...
    'baInfo', {'RPCnp','LOA','SD','ks'}, ...
    'data1mode', 'truth', ...
    'showFitCI', 'on', ...
    'legend', 'off', ...
    'axesLimits', 'tight', ...
    'baStatsMode','Non-parametric');

[ rpc2, ah2, stats2 ] = BlandAltman(ah2, ...
    CC_HR_DATA.iHr, CC_HR_DATA.rotWsstHr, ...
    {'ifH','rotWsstH'}, ...
    'Comparison of rotational WSST and instantaneous heart rate estimates (beats · min^-^1)', ...
    {'ifH:rotWsstH','LoBF'}, ...
    'corrInfo', {'n','eq','r','r2','p'}, ...
    'baInfo', {'RPCnp','LOA','SD','ks'}, ...
    'data1mode', 'truth', ...
    'showFitCI', 'on', ...
    'legend', 'off', ...
    'axesLimits', 'tight', ...
    'baStatsMode','Non-parametric');

[ rpc3, ah3, stats3 ] = BlandAltman(ah3, ...
    CC_HR_DATA.iHr, CC_HR_DATA.ensWsstHr, ...
	{'ifH','ensWsstH'}, ...
    'Comparison of ensemble WSST and instantaneous heart rate estimates (beats · min^-^1)', ...
    {'ifH:linWsstH','LoBF'}, ...
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

%% make box plots of the distributions of ifh and various wsstH 

figure('Color', white);

b1a = subplot(221);
histogram(CC_HR_DATA.iHr, 'Normalization', 'pdf', 'BinWidth', 5, ...
    'NumBins', 50, 'FaceColor', [0.6 0.6 0.6] ); 
hold on;
mu = mean(CC_HR_DATA.iHr);
sigma = std(CC_HR_DATA.iHr);
YY = 0:0.1:150;
f = exp(-(YY-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(YY, f, '--', 'Color', maroon, 'LineWidth', 2);
hold off;
xlabel('beats · min^-^1');
grid;
legend('ifH','normPdf');

b1b = subplot(222);
histogram(CC_HR_DATA.ensWsstHr, 'Normalization', 'pdf', ...
    'BinWidth', 5, 'NumBins', 50, 'FaceColor', blue);
hold on;
mu = mean(CC_HR_DATA.ensWsstHr);
sigma = std(CC_HR_DATA.ensWsstHr);
YY = 0:0.1:150;
f = exp(-(YY-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(YY, f, '--', 'Color', maroon, 'LineWidth', 2);
hold off;
xlabel('beats · min^-^1');
grid;
legend('ensWsstH','normPdf');

b1c = subplot(223);
histogram(CC_HR_DATA.linWsstHr, 'Normalization', 'pdf', ...
    'BinWidth', 5, 'NumBins', 50, 'FaceColor', red );
hold on;
mu = mean(CC_HR_DATA.linWsstHr);
sigma = std(CC_HR_DATA.linWsstHr);
YY = 0:0.1:150;
f = exp(-(YY-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(YY, f, '--', 'Color', maroon, 'LineWidth', 2);
hold off;
xlabel('beats · min^-^1');
grid;
legend('linWsstH','normPdf');

b1d = subplot(224);
histogram(CC_HR_DATA.rotWsstHr, 'Normalization', 'pdf', ...
    'BinWidth', 5, 'NumBins', 50, 'FaceColor', green );
hold on;
mu = mean(CC_HR_DATA.linWsstHr);
sigma = std(CC_HR_DATA.linWsstHr);
YY = 0:0.1:150;
f = exp(-(YY-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(YY, f, '--', 'Color', maroon, 'LineWidth', 2);
hold off;
xlabel('beats · min^-^1');
grid;
legend('rotWsstHr','normPdf');

linkaxes([b1a b1b b1c b1d],'x');

%% Look at epoch ~= 3 only

figure('Color', white);

b1a = subplot(221);
histogram(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 3), 'Normalization', 'pdf', 'BinWidth', 5, ...
    'NumBins', 50, 'FaceColor', [0.6 0.6 0.6] ); 
hold on;
mu = mean(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 3));
sigma = std(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 3));
YY = 0:0.1:150;
f = exp(-(YY-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(YY, f, '--', 'Color', maroon, 'LineWidth', 2);
hold off;
xlabel('beats · min^-^1');
grid;
legend('ifH','normPdf');

b1b = subplot(222);
histogram(CC_HR_DATA.ensWsstHr(CC_HR_DATA.epoch == 3), 'Normalization', 'pdf', ...
    'BinWidth', 5, 'NumBins', 50, 'FaceColor', blue);
hold on;
mu = mean(CC_HR_DATA.ensWsstHr(CC_HR_DATA.epoch == 3));
sigma = std(CC_HR_DATA.ensWsstHr(CC_HR_DATA.epoch == 3));
YY = 0:0.1:150;
f = exp(-(YY-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(YY, f, '--', 'Color', maroon, 'LineWidth', 2);
hold off;
xlabel('beats · min^-^1');
grid;
legend('ensWsstH','normPdf');

b1c = subplot(223);
histogram(CC_HR_DATA.linWsstHr(CC_HR_DATA.epoch == 3), 'Normalization', 'pdf', ...
    'BinWidth', 5, 'NumBins', 50, 'FaceColor', red );
hold on;
mu = mean(CC_HR_DATA.linWsstHr(CC_HR_DATA.epoch == 3));
sigma = std(CC_HR_DATA.linWsstHr(CC_HR_DATA.epoch == 3));
YY = 0:0.1:150;
f = exp(-(YY-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(YY, f, '--', 'Color', maroon, 'LineWidth', 2);
hold off;
xlabel('beats · min^-^1');
grid;
legend('linWsstH','normPdf');

b1d = subplot(224);
histogram(CC_HR_DATA.rotWsstHr(CC_HR_DATA.epoch == 3), 'Normalization', 'pdf', ...
    'BinWidth', 5, 'NumBins', 50, 'FaceColor', green );
hold on;
mu = mean(CC_HR_DATA.linWsstHr(CC_HR_DATA.epoch == 3));
sigma = std(CC_HR_DATA.linWsstHr(CC_HR_DATA.epoch == 3));
YY = 0:0.1:150;
f = exp(-(YY-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(YY, f, '--', 'Color', maroon, 'LineWidth', 2);
hold off;
xlabel('beats · min^-^1');
grid;
legend('rotWsstHr','normPdf');

linkaxes([b1a b1b b1c b1d],'x');

%% make normalized box plots with 

figure('Color', white);

histogram(CC_HR_DATA.iHr, 'Normalization', 'pdf');
mu = mean(CC_HR_DATA.iHr);
sigma = std(CC_HR_DATA.iHr);

hold on;
mu = mean(CC_HR_DATA.iHr);
sigma = std(CC_HR_DATA.iHr);
YY = 0:0.1:150;
f = exp(-(YY-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(YY,f,'LineWidth',2);
hold off;

%%

t1 = E1_KE.Time(1);
t2 = E1_KE.Time(end);
tRange = timerange(t1,t2);

figure('Color', white);

s1 = subplot(221);
plot(KINEMATICS.Time(tRange), KINEMATICS.ax(tRange), 'Color', blue );
xlabel('Time, local');
ylabel('m · sec^-^2');
title('Unfiltered x-axis accelerometer');
grid;

s2 = subplot(222);
plot(KINEMATICS.Time(tRange), KINEMATICS.gy(tRange) , 'Color', red);
xlabel('Time, local');
ylabel('deg · sec^-^1');
title('Unfiltered y-axis gyroscope');
grid;

s3 = subplot(223);
plot(E1_KE.Time, E1_KE.bpAx, 'Color', blue);
xlabel('Time, local');
ylabel('m · sec^-^2');
title('Band pass-filtered x-axis accelerometer');
grid;

s4 = subplot(224);
plot(E1_KE.Time, E1_KE.bpGy, 'Color', red);
xlabel('Time, local');
ylabel('deg · sec^-^1');
title('Bandpass-filtered y-axis gyroscope');
grid;

linkaxes([s1 s2 s3 s4], 'x');

%% 
 
HOKU = '63HF';
LIHO = '01L5';
KOLOHE = '6JK5';
HUA = '83H1';
NOA = '9ON6';
LONO = '9FL3';

HokuStahlHr = ...
    unique(241 * HR_DATA.mass(string(HR_DATA.id) == HOKU) .^ -0.25);
LihoStahlHr = ...
    unique(241 * HR_DATA.mass(string(HR_DATA.id) == LIHO) .^ -0.25);
KoloheStahlHr = ...
    unique(241 * HR_DATA.mass(string(HR_DATA.id) == KOLOHE) .^ -0.25);
HuaStahlHr = ...
    unique(241 * HR_DATA.mass(string(HR_DATA.id) == HUA) .^ -0.25);
NoaStahlHr = ...
    unique(241 * HR_DATA.mass(string(HR_DATA.id) == NOA) .^ -0.25);
LonoStahlHr = ...
    unique(241 * HR_DATA.mass(string(HR_DATA.id) == LONO) .^ -0.25);

meanStahlHr = mean( [ HokuStahlHr LihoStahlHr KoloheStahlHr ...
    HuaStahlHr NoaStahlHr LonoStahlHr] );


HOKU_E1_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 1 & contains(HR_DATA.id, {HOKU}) );
HOKU_E3_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 3 & contains(HR_DATA.id, {HOKU}) );

LIHO_E1_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 1 & contains(HR_DATA.id, {LIHO}) );
LIHO_E3_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 3 & contains(HR_DATA.id, {LIHO}) );
LIHO_E5_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 5 & contains(HR_DATA.id, {LIHO}) );

KOLOHE_E1_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 1 & contains(HR_DATA.id, {KOLOHE}) );
KOLOHE_E3_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 3 & contains(HR_DATA.id, {KOLOHE}) );
KOLOHE_E5_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 5 & contains(HR_DATA.id, {KOLOHE}) );   

HUA_E1_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 1 & contains(HR_DATA.id, {HUA}) );
HUA_E3_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 3 & contains(HR_DATA.id, {HUA}) );
HUA_E5_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 5 & contains(HR_DATA.id, {HUA}) );       

NOA_E1_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 1 & contains(HR_DATA.id, {NOA}) );
NOA_E3_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 3 & contains(HR_DATA.id, {NOA}) );
NOA_E5_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 5 & contains(HR_DATA.id, {NOA}) );       

LONO_E1_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 1 & contains(HR_DATA.id, {LONO}) );
LONO_E3_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 3 & contains(HR_DATA.id, {LONO}) );

%% Compile stats

% ========= Hoku stats ==========

mean_Hoku_E1_wsstHR = mean(HOKU_E1_HR);
mean_Hoku_E3_wsstHR = mean(HOKU_E3_HR);
mean_Hoku_E1_iHR = mean(HOKU_E1_CC_DATA.iHR);
mean_Hoku_E3_iHR = mean(HOKU_E3_CC_DATA.iHR);

std_Hoku_E1_wsstHR  = std(HOKU_E1_HR);
std_Hoku_E3_wsstHR  = std(HOKU_E3_HR);
std_Hoku_E1_iHR  = std(HOKU_E1_CC_DATA.iHR);
std_Hoku_E3_iHR  = std(HOKU_E3_CC_DATA.iHR);

min_Hoku_E1_wsstHR = min(HOKU_E1_HR);
max_Hoku_E1_wsstHR = max(HOKU_E1_HR);
min_Hoku_E3_wsstHR = min(HOKU_E3_HR);
max_Hoku_E3_wsstHR = max(HOKU_E3_HR);

min_Hoku_E1_iHR = min(HOKU_E1_CC_DATA.iHR);
max_Hoku_E1_iHR = max(HOKU_E1_CC_DATA.iHR);
min_Hoku_E3_iHR = min(HOKU_E3_CC_DATA.iHR);
max_Hoku_E3_iHR = max(HOKU_E3_CC_DATA.iHR);

meanIBI_Hoku_E1 = mean(HOKU_E1_VENT_DATA.ibi);
stdIBI_Hoku_E1  = std(HOKU_E1_VENT_DATA.ibi);
meanIBI_Hoku_Eall= mean(VENT_HR_DATA.ibi(VENT_HR_DATA.altId == 1));
stdIBI_Hoku_Eall  = std(VENT_HR_DATA.ibi(VENT_HR_DATA.altId == 1));

HOKU_E3_TIME = HR_DATA.Time(HR_DATA.epoch == 3 & HR_DATA.altId == 1);
apneaDuration_Hoku_E3 = seconds( HOKU_E3_TIME(end) - HOKU_E3_TIME(1) );


% ========= Liho stats ==========

mean_Liho_E1_wsstHR = mean(LIHO_E1_HR);
mean_Liho_E3_wsstHR = mean(LIHO_E3_HR);
mean_Liho_E5_wsstHR = mean(LIHO_E5_HR);
mean_Liho_E1_iHR = mean(LIHO_E1_CC_DATA.iHR);
mean_Liho_E3_iHR = mean(LIHO_E3_CC_DATA.iHR);
mean_Liho_E5_iHR = mean(LIHO_E5_CC_DATA.iHR);

std_Liho_E1_wsstHR  = std(LIHO_E1_HR);
std_Liho_E3_wsstHR  = std(LIHO_E3_HR);
std_Liho_E5_wsstHR  = std(LIHO_E5_HR);
std_Liho_E1_iHR  = std(LIHO_E1_CC_DATA.iHR);
std_Liho_E3_iHR  = std(LIHO_E3_CC_DATA.iHR);
std_Liho_E5_iHR  = std(LIHO_E5_CC_DATA.iHR);

min_Liho_E1_wsstHR = min(LIHO_E1_HR);
max_Liho_E1_wsstHR = max(LIHO_E1_HR);
min_Liho_E3_wsstHR = min(LIHO_E3_HR);
max_Liho_E3_wsstHR = max(LIHO_E3_HR);
min_Liho_E5_wsstHR = min(LIHO_E5_HR);
max_Liho_E5_wsstHR = max(LIHO_E5_HR);

min_Liho_E1_iHR = min(LIHO_E1_CC_DATA.iHR);
max_Liho_E1_iHR = max(LIHO_E1_CC_DATA.iHR);
min_Liho_E3_iHR = min(LIHO_E3_CC_DATA.iHR);
max_Liho_E3_iHR = max(LIHO_E3_CC_DATA.iHR);
min_Liho_E5_iHR = min(LIHO_E5_CC_DATA.iHR);
max_Liho_E5_iHR = max(LIHO_E5_CC_DATA.iHR);

meanIBI_Liho_E1 = mean(LIHO_E1_VENT_DATA.ibi);
meanIBI_Liho_E5 = mean(LIHO_E5_VENT_DATA.ibi);
stdIBI_Liho_E1 = std(LIHO_E1_VENT_DATA.ibi);
stdIBI_Liho_E5 = std(LIHO_E5_VENT_DATA.ibi);

meanIBI_Liho_Eall= mean(VENT_HR_DATA.ibi(VENT_HR_DATA.altId == 2));
stdIBI_Liho_Eall  = std(VENT_HR_DATA.ibi(VENT_HR_DATA.altId == 2));

LIHO_E3_TIME = HR_DATA.Time(HR_DATA.epoch == 3 & HR_DATA.altId == 2);
apneaDuration_Liho_E3 = seconds( LIHO_E3_TIME(end) - LIHO_E3_TIME(1) );


% ========= Kolohe stats ==========

mean_Kolohe_E1_HR = mean(KOLOHE_E1_HR);
mean_Kolohe_E3_HR = mean(KOLOHE_E3_HR);
mean_Kolohe_E5_HR = mean(KOLOHE_E5_HR);
mean_Kolohe_E1_iHR = mean(KOLOHE_E1_CC_DATA.iHR);
mean_Kolohe_E3_iHR = mean(KOLOHE_E3_CC_DATA.iHR);
mean_Kolohe_E5_iHR = mean(KOLOHE_E5_CC_DATA.iHR);

std_Kolohe_E1_wsstHR  = std(KOLOHE_E1_HR);
std_Kolohe_E3_wsstHR  = std(KOLOHE_E3_HR);
std_Kolohe_E5_wsstHR  = std(KOLOHE_E5_HR);
std_Kolohe_E1_iHR  = std(KOLOHE_E1_CC_DATA.iHR);
std_Kolohe_E3_iHR  = std(KOLOHE_E3_CC_DATA.iHR);
std_Kolohe_E5_iHR  = std(KOLOHE_E5_CC_DATA.iHR);

min_Kolohe_E1_wsstHR = min(KOLOHE_E1_HR);
max_Kolohe_E1_wsstHR = max(KOLOHE_E1_HR);
min_Kolohe_E3_wsstHR = min(KOLOHE_E3_HR);
max_Kolohe_E3_wsstHR = max(KOLOHE_E3_HR);
min_Kolohe_E5_wsstHR = min(KOLOHE_E5_HR);
max_Kolohe_E5_wsstHR = max(KOLOHE_E5_HR);

min_Kolohe_E1_iHR = min(KOLOHE_E1_CC_DATA.iHR);
max_Kolohe_E1_iHR = max(KOLOHE_E1_CC_DATA.iHR);
min_Kolohe_E3_iHR = min(KOLOHE_E3_CC_DATA.iHR);
max_Kolohe_E3_iHR = max(KOLOHE_E3_CC_DATA.iHR);
min_Kolohe_E5_iHR = min(KOLOHE_E5_CC_DATA.iHR);
max_Kolohe_E5_iHR = max(KOLOHE_E5_CC_DATA.iHR);

meanIBI_Kolohe_E1 = mean(KOLOHE_E1_VENT_DATA.ibi);
meanIBI_Kolohe_E5 = mean(KOLOHE_E5_VENT_DATA.ibi);
stdIBI_Kolohe_E1 = std(KOLOHE_E1_VENT_DATA.ibi);
stdIBI_Kolohe_E5 = std(KOLOHE_E5_VENT_DATA.ibi);

meanIBI_Kolohe_Eall= mean(VENT_HR_DATA.ibi(VENT_HR_DATA.altId == 3));
stdIBI_Kolohe_Eall  = std(VENT_HR_DATA.ibi(VENT_HR_DATA.altId == 3));

KOLOHE_E3_TIME = HR_DATA.Time(HR_DATA.epoch == 3 & HR_DATA.altId == 3);
apneaDuration_Kolohe_E3 = seconds( KOLOHE_E3_TIME(end) - KOLOHE_E3_TIME(1) );


% ========= Hua stats ==========

mean_Hua_E1_wsstHR = mean(HUA_E1_HR);
mean_Hua_E3_wsstHR = mean(HUA_E3_HR);
mean_Hua_E5_wsstHR = mean(HUA_E5_HR);
mean_Hua_E1_iHR = mean(HUA_E1_CC_DATA.iHR);
mean_Hua_E3_iHR = mean(HUA_E3_CC_DATA.iHR);
mean_Hua_E5_iHR = mean(HUA_E5_CC_DATA.iHR);

std_Hua_E1_wsstHR  = std(HUA_E1_HR);
std_Hua_E3_wsstHR  = std(HUA_E3_HR);
std_Hua_E5_wsstHR  = std(HUA_E5_HR);
std_Hua_E1_iHR  = std(HUA_E1_CC_DATA.iHR);
std_Hua_E3_iHR  = std(HUA_E3_CC_DATA.iHR);
std_Hua_E5_iHR  = std(HUA_E5_CC_DATA.iHR);

min_Hua_E1_wsstHR = min(HUA_E1_HR);
max_Hua_E1_wsstHR = max(HUA_E1_HR);
min_Hua_E3_wsstHR = min(HUA_E3_HR);
max_Hua_E3_wsstHR = max(HUA_E3_HR);
min_Hua_E5_wsstHR = min(HUA_E5_HR);
max_Hua_E5_wsstHR = max(HUA_E5_HR);

min_Hua_E1_iHR = min(HUA_E1_CC_DATA.iHR);
max_Hua_E1_iHR = max(HUA_E1_CC_DATA.iHR);
min_Hua_E3_iHR = min(HUA_E3_CC_DATA.iHR);
max_Hua_E3_iHR = max(HUA_E3_CC_DATA.iHR);
min_Hua_E5_iHR = min(HUA_E5_CC_DATA.iHR);
max_Hua_E5_iHR = max(HUA_E5_CC_DATA.iHR);

meanIBI_Hua_E1 = mean(HUA_E1_VENT_DATA.ibi);
meanIBI_Hua_E5 = mean(HUA_E5_VENT_DATA.ibi);
stdIBI_Hua_E1 = std(HUA_E1_VENT_DATA.ibi);
stdIBI_Hua_E5 = std(HUA_E5_VENT_DATA.ibi);

meanIBI_Hua_Eall= mean(VENT_HR_DATA.ibi(VENT_HR_DATA.altId == 4));
stdIBI_Hua_Eall  = std(VENT_HR_DATA.ibi(VENT_HR_DATA.altId == 4));

HUA_E3_TIME = HR_DATA.Time(HR_DATA.epoch == 3 & HR_DATA.altId == 4);
apneaDuration_Hua_E3 = seconds( HUA_E3_TIME(end) - HUA_E3_TIME(1) );


% ========= Noa stats ==========

mean_Noa_E1_wsstHR = mean(NOA_E1_HR);
mean_Noa_E3_wsstHR = mean(NOA_E3_HR);
mean_Noa_E5_wsstHR = mean(NOA_E5_HR);
mean_Noa_E1_iHR = mean(NOA_E1_CC_DATA.iHR);
mean_Noa_E3_iHR = mean(NOA_E3_CC_DATA.iHR);
mean_Noa_E5_iHR = mean(NOA_E5_CC_DATA.iHR);

std_Noa_E1_wsstHR  = std(NOA_E1_HR);
std_Noa_E3_wsstHR  = std(NOA_E3_HR);
std_Noa_E5_wsstHR  = std(NOA_E5_HR);
std_Noa_E1_iHR  = std(NOA_E1_CC_DATA.iHR);
std_Noa_E3_iHR  = std(NOA_E3_CC_DATA.iHR);
std_Noa_E5_iHR  = std(NOA_E5_CC_DATA.iHR);

min_Noa_E1_wsstHR = min(NOA_E1_HR);
max_Noa_E1_wsstHR = max(NOA_E1_HR);
min_Noa_E3_wsstHR = min(NOA_E3_HR);
max_Noa_E3_wsstHR = max(NOA_E3_HR);
min_Noa_E5_wsstHR = min(NOA_E5_HR);
max_Noa_E5_wsstHR = max(NOA_E5_HR);

min_Noa_E1_iHR = min(NOA_E1_CC_DATA.iHR);
max_Noa_E1_iHR = max(NOA_E1_CC_DATA.iHR);
min_Noa_E3_iHR = min(NOA_E3_CC_DATA.iHR);
max_Noa_E3_iHR = max(NOA_E3_CC_DATA.iHR);
min_Noa_E5_iHR = min(NOA_E5_CC_DATA.iHR);
max_Noa_E5_iHR = max(NOA_E5_CC_DATA.iHR);

meanIBI_Noa_E1 = mean(NOA_E1_VENT_DATA.ibi);
meanIBI_Noa_E5 = mean(NOA_E5_VENT_DATA.ibi);
stdIBI_Noa_E1 = std(NOA_E1_VENT_DATA.ibi);
stdIBI_Noa_E5 = std(NOA_E5_VENT_DATA.ibi);

meanIBI_Noa_Eall= mean(VENT_HR_DATA.ibi(VENT_HR_DATA.altId == 5));
stdIBI_Noa_Eall  = std(VENT_HR_DATA.ibi(VENT_HR_DATA.altId == 5));

NOA_E3_TIME = HR_DATA.Time(HR_DATA.epoch == 3 & HR_DATA.altId == 5);
apneaDuration_Noa_E3 = seconds( NOA_E3_TIME(end) - NOA_E3_TIME(1) );


% ========= Lono stats ==========

mean_Lono_E1_wsstHR = mean(LONO_E1_HR);
mean_Lono_E3_wsstHR = mean(LONO_E3_HR);
mean_Lono_E1_iHR = mean(LONO_E1_CC_DATA.iHR);
mean_Lono_E3_iHR = mean(LONO_E3_CC_DATA.iHR);

std_Lono_E1_wsstHR  = std(LONO_E1_HR);
std_Lono_E3_wsstHR  = std(LONO_E3_HR);
std_Lono_E1_iHR  = std(LONO_E1_CC_DATA.iHR);
std_Lono_E3_iHR  = std(LONO_E3_CC_DATA.iHR);

min_Lono_E1_wsstHR = min(LONO_E1_HR);
max_Lono_E1_wsstHR = max(LONO_E1_HR);
min_Lono_E3_wsstHR = min(LONO_E3_HR);
max_Lono_E3_wsstHR = max(LONO_E3_HR);
min_Lono_E1_iHR = min(LONO_E1_CC_DATA.iHR);
max_Lono_E1_iHR = max(LONO_E1_CC_DATA.iHR);
min_Lono_E3_iHR = min(LONO_E3_CC_DATA.iHR);
max_Lono_E3_iHR = max(LONO_E3_CC_DATA.iHR);

meanIBI_Lono_E1 = mean(LONO_E1_VENT_DATA.ibi);
stdIBI_Lono_E1 = std(LONO_E1_VENT_DATA.ibi);

meanIBI_Lono_Eall= mean(VENT_HR_DATA.ibi(VENT_HR_DATA.altId == 6));
stdIBI_Lono_Eall  = std(VENT_HR_DATA.ibi(VENT_HR_DATA.altId == 6));

LONO_E3_TIME = HR_DATA.Time(HR_DATA.epoch == 3 & HR_DATA.altId == 6);
apneaDuration_Lono_E3 = seconds( LONO_E3_TIME(end) - LONO_E3_TIME(1) );

% calculate grand means, std, and min/max for iHr

grandMean_iHR = mean(CC_HR_DATA.iHr);
grandStd_iHR  = std(CC_HR_DATA.iHr);

grandMean_iHR_E1 = mean(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 1));
grandStd_iHR_E1  = std(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 1));
grandMin_iHr_E1 = min(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 1));
grandMax_iHr_E1 = max(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 1));

grandMean_iHR_E3 = mean(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 3));
grandStd_iHR_E3  = std(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 3));
grandMin_iHr_E3 = min(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 3));
grandMax_iHr_E3 = max(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 3));

grandMean_iHR_E5 = mean(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 5));
grandStd_iHR_E5  = std(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 5));
grandMin_iHr_E5 = min(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 5));
grandMax_iHr_E5 = max(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 5));

% calculate grand means, std, and min/max for ensWsstHr

grandMean_wsstHr = mean(WSST_HR_DATA.ensWsstHr);
grandStd_wsstHr  = std(WSST_HR_DATA.ensWsstHr);

grandMean_wsstHr_E1 = mean(WSST_HR_DATA.ensWsstHr(WSST_HR_DATA.epoch == 1));
grandStd_wsstHr_E1  = std(WSST_HR_DATA.ensWsstHr(WSST_HR_DATA.epoch == 1));
grandMin_wsstHr_E1 = min(WSST_HR_DATA.ensWsstHr(WSST_HR_DATA.epoch == 1));
grandMax_wsstHr_E1 = max(WSST_HR_DATA.ensWsstHr(WSST_HR_DATA.epoch == 1));

grandMean_wsstHr_E3 = mean(WSST_HR_DATA.ensWsstHr(WSST_HR_DATA.epoch == 3));
grandStd_wsstHr_E3  = std(WSST_HR_DATA.ensWsstHr(WSST_HR_DATA.epoch == 3));
grandMin_wsstHr_E3 = min(WSST_HR_DATA.ensWsstHr(WSST_HR_DATA.epoch == 3));
grandMax_wsstHr_E3 = max(WSST_HR_DATA.ensWsstHr(WSST_HR_DATA.epoch == 3));

grandMean_wsstHr_E5 = mean(WSST_HR_DATA.ensWsstHr(WSST_HR_DATA.epoch == 5));
grandStd_wsstHr_E5  = std(WSST_HR_DATA.ensWsstHr(WSST_HR_DATA.epoch == 5));
grandMin_wsstHr_E5 = min(WSST_HR_DATA.ensWsstHr(WSST_HR_DATA.epoch == 5));
grandMax_wsstHr_E5 = max(WSST_HR_DATA.ensWsstHr(WSST_HR_DATA.epoch == 5));

%

grandMean_iHR_E1_E5 = mean(CC_HR_DATA.iHr(CC_HR_DATA.epoch ~= 3));
grandStd_iHR_E1_E5  = std(CC_HR_DATA.iHr(CC_HR_DATA.epoch ~= 3));
grandMin_iHr_E1_E5 = min(CC_HR_DATA.iHr(CC_HR_DATA.epoch ~= 3));
grandMax_iHr_E1_E5 = max(CC_HR_DATA.iHr(CC_HR_DATA.epoch ~= 3));

grandMean_wsstHr_E1_E5 = mean(WSST_HR_DATA.ensWsstHr(WSST_HR_DATA.epoch ~= 3));
grandStd_wsstHr_E1_E5  = std(WSST_HR_DATA.ensWsstHr(WSST_HR_DATA.epoch ~= 3));
grandMin_wsstHr_E1_E5 = min(WSST_HR_DATA.ensWsstHr(WSST_HR_DATA.epoch ~= 3));
grandMax_wsstHr_E1_E5 = max(WSST_HR_DATA.ensWsstHr(WSST_HR_DATA.epoch ~= 3));



%



%% Make a nice table of all these values...
%                                        HOKU  LIHO  KOLOHE  HUA  NOA  LONO
% mean_E1_wsstHr    w/ std_E1_wsstHr
% mean_E3_wsstHr    w/ std_E3_wsstHr
% mean_E5_wsstHr    w/ std_E5_wsstHr
% mean_E1_iHr       w/ std_E1_iHr
% mean_E3_iHr       w/ std_E3_iHr
% mean_E5_iHr       w/ std_E5_iHr

% min_E1_wsstHr
% min_E3_wsstHr
% min_E5_wsstHr
% min_E1_iHr
% min_E3_iHr
% min_E5_iHr

% max_E1_wsstHr
% max_E3_wsstHr
% max_E5_wsstHr
% max_E1_iHr
% max_E3_iHr
% max_E5_iHr

% mean_E1_ibi       w/ std_E1_ibi  
% mean_E3_ibi       w/ std_E3_ibi
% mean_E5_ibi       w/ std_E5_ibi



%% IBI, mean ifH, min ifH, and max ifH by animal ID

figure('Color', [1 1 1]);

s1 = subplot(141);
boxplot( VENT_DATA.ibi, VENT_DATA.id, ...
    'Labels', {'63HF','01L5','6JK5','83H1','9ON6','9FL3'} );
xlabel('Animal ID');
ylabel('Time, seconds');
title('Inter-breath interval');
grid;
s2 = subplot(142);
boxplot( VENT_DATA.mean_ifH, VENT_DATA.id );
xlabel('Animal ID');
ylabel('beats · min^-^1');
title('Mean instantaneous heart rate');
grid;
s3 = subplot(143);
boxplot( VENT_DATA.min_ifH, VENT_DATA.id );
xlabel('Animal ID');
ylabel('beats · min^-^1');
title('Minimum instantaneous heart rate');
grid;
s4 = subplot(144);
boxplot( VENT_DATA.max_ifH, VENT_DATA.id);
xlabel('Animal ID');
ylabel('beats · min^-^1');
title('Maximum instantaneous heart rate');
grid;

%% IBI, mean ifH, min ifH, and max ifH by epoch

figure('Color', [1 1 1]);

s1 = subplot(141);
boxplot( VENT_DATA.ibi, VENT_DATA.epoch, ...
    'Labels', {'Baseline','Apnea','Recovery'} );
xlabel('Epoch');
ylabel('Time, seconds');
title('Inter-breath interval');
grid;
s2 = subplot(142);
boxplot( VENT_DATA.mean_ifH, VENT_DATA.epoch, ...
    'Labels', {'Baseline','Apnea','Recovery'} );
xlabel('Epoch');
ylabel('beats · min^-^1');
title('Mean instantaneous heart rate');
grid;
s3 = subplot(143);
boxplot( VENT_DATA.min_ifH, VENT_DATA.epoch, ...
    'Labels', {'Baseline','Apnea','Recovery'} );
xlabel('Epoch');
ylabel('beats · min^-^1');
title('Minimum instantaneous heart rate');
grid;
s4 = subplot(144);
boxplot( VENT_DATA.max_ifH, VENT_DATA.epoch, ...
    'Labels', {'Baseline','Apnea','Recovery'} );
xlabel('Epoch');
ylabel('beats · min^-^1');
title('Maximum instantaneous heart rate');
grid;

%%

figure('Color',[1 1 1]);

subplot(231);
boxplot(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 1), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 1), 'Notch','on');
xlabel('Animal ID');
ylabel('beats · min^-^1');
ylim([30 150]);
title('Baseline Epoch (E1): ifH');
grid;
subplot(232);
boxplot(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 3), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 3), 'Notch','on');
xlabel('Animal ID');
ylabel('beats · min^-^1');
ylim([30 150]);
title('Apnea Epoch (E3): ifH');
grid;
subplot(233);
boxplot(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 5), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 5), 'Notch','on');
xlabel('Animal ID');
ylabel('beats · min^-^1');
ylim([30 150]);
title('Recovery Epoch (E5): ifH');
grid;


subplot(234);
boxplot(CC_HR_DATA.ensWsstHr(CC_HR_DATA.epoch == 1), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 1), 'Notch','on');
xlabel('Animal ID');
ylabel('beats · min^-^1');
ylim([30 150]);
title('Baseline Epoch (E1): wsstfH');
grid;

subplot(235);
boxplot(CC_HR_DATA.ensWsstHr(CC_HR_DATA.epoch == 3), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 3), 'Notch','on');
xlabel('Animal ID');
ylabel('beats · min^-^1');
ylim([30 150]);
title('Apnea Epoch (E3): wsstfH');
grid;

subplot(236);
boxplot(CC_HR_DATA.ensWsstHr(CC_HR_DATA.epoch == 5), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 5), 'Notch','on');
xlabel('Animal ID');
ylabel('beats · min^-^1');
ylim([30 150]);
title('Recovery Epoch (E5): wsstfH');
grid;


%%


%% Run mixed-effects ANOVA: ensWsstHr ~ fixed(epoch) + random(id)

% a two-way ANOVA with unbalanced design to compare heart rates during
% baseline, apnea, and recovery epochs

% mixed-effects (two-way) ANOVA: ensWsstHr ~ fixed(epoch) + random(id)

y1a = CC_HR_DATA.iHr;
y2a = CC_HR_DATA.ensWsstHr;
yfH = [y1a, y2a];
g1a = CC_HR_DATA.epoch;
g2a = CC_HR_DATA.altId;

[pAOV1a, tblAOV1a, statsAOV1a, termsAOV1a] = ...
    anovan( y1a, {g1a, g2a}, ...
    'model', 1, 'random', 2, ... %'interaction', ...
    'varnames', {'Epoch','AnimalID'} );

[pAOV1b, tblAOV1b, statsAOV1b, termsAOV1b] = ...
    anovan( y2a, {g1a, g2a}, ...
    'model', 1, 'random', 2, ... %'interaction', ...
    'varnames', {'Epoch','AnimalID'} );


figure('Color', white);
[cAOV1a, mAOV1a, hAOV1a, nmsAOV1a] = ...
    multcompare(statsAOV1a, 'CType', 'hsd'); %,'Dimension',[1 2]);

figure('Color', white);
[cAOV1b, mAOV1b, hAOV1b, nmsAOV1b] = ...
    multcompare(statsAOV1b, 'CType', 'hsd'); %,'Dimension',[1 2]);


% mixed-effects ANOVA (same as 1a) without Hoku (special dolphin, missed
% E5)

y1b = CC_HR_DATA.ensWsstHr(CC_HR_DATA.altId ~= 1);
y2b = CC_HR_DATA.iHr(CC_HR_DATA.altId ~= 1);
g1b = CC_HR_DATA.epoch(CC_HR_DATA.altId ~= 1);
g2b = CC_HR_DATA.altId(CC_HR_DATA.altId ~= 1);

[pAOV1b, tblAOV1b, statsAOV1b, termsAOV1b] = ...
    anovan( y2b, {g1b, g2b}, ...
    'model', 1, 'random', 2, ... %'interaction', ...
    'varnames', {'Epoch','AnimalID'} );

figure('Color', white);
[cAOV1b, mAOV1b, hAOV1b, nmsAOV1b] = ...
    multcompare(statsAOV1b, 'CType', 'hsd'); %,'Dimension',[1 2]);

% mixed-effects ANOVA (same as 1a and 1b) without Hoku & Lono

y1c = CC_HR_DATA.ensWsstHr(CC_HR_DATA.altId ~= 1 & CC_HR_DATA.altId ~= 6) ;
y2c = CC_HR_DATA.iHr(CC_HR_DATA.altId ~= 1 & CC_HR_DATA.altId ~= 6) ;
g1c = CC_HR_DATA.epoch(CC_HR_DATA.altId ~= 1 & CC_HR_DATA.altId ~= 6) ;
g2c = CC_HR_DATA.altId(CC_HR_DATA.altId ~= 1 & CC_HR_DATA.altId ~= 6) ;

[pAOV1c, tblAOV1c, statsAOV1c, termsAOV1c] = ...
    anovan( y1c, {g1c, g2c}, ...
    'model', 1, 'random', 2, ... %'interaction', ...
    'varnames', {'Epoch','AnimalID'} );

[pAOV1d, tblAOV1d, statsAOV1d, termsAOV1d] = ...
    anovan( y2c, {g1c, g2c}, ...
    'model', 1, 'random', 2, ... %'interaction', ...
    'varnames', {'Epoch','AnimalID'} );

figure('Color', white);
[cAOV1c, mAOV1c, hAOV1c, nmsAOV1c] = ...
    multcompare(statsAOV1c, 'CType', 'hsd'); %,'Dimension',[1 2]);

figure('Color', white);
[cAOV1d, mAOV1d, hAOV1d, nmsAOV1d] = ...
    multcompare(statsAOV1d, 'CType', 'hsd'); %,'Dimension',[1 2]);

% two-way ANOVA: ensWsstHr ~ fixed(epoch * id)

y1a = CC_HR_DATA.iHr;
y2a = CC_HR_DATA.ensWsstHr;
g1a = CC_HR_DATA.epoch;
g2a = CC_HR_DATA.altId;

[pAOV1e, tblAOV1e, statsAOV1e, termsAOV1e] = ...
    anovan( y1a, {g1a, g2a}, ...
    'model', 'interaction', ...
    'varnames', {'Epoch','AnimalID'} );

[pAOV1f, tblAOV1f, statsAOV1f, termsAOV1f] = ...
    anovan( y2a, {g1a, g2a}, ...
    'model', 'interaction', ...
    'varnames', {'Epoch','AnimalID'} );

figure('Color', white);
[cAOV1e, mAOV1e, hAOV1e, nmsAOV1e] = ...
    multcompare(statsAOV1e, 'CType', 'hsd','Dimension',[1 2]);

figure('Color', white);
[cAOV1f, mAOV1f, hAOV1f, nmsAOV1f] = ...
    multcompare(statsAOV1f, 'CType', 'hsd','Dimension',[1 2]);



%% 

figure('Color', white);

s1 = subplot(211);
plot(E1_KE.Time, E1_KE.bpAx, 'Color', blue, 'DisplayName', 'bpAx');
hold on;
plot(E1_VENT_DATA.Time, 0.4, 'v', 'MarkerEdgeColor', black, ...
    'MarkerFaceColor', goldenrod ); 
xline(E1_VENT_DATA.Time, ':', 'Color', goldenrod, 'DisplayName', '');
hold off;
xlabel('Time, local');
ylabel('m · s^-^2');
grid;
legend('bpAx','resp');

s2 = subplot(212);
plot(E1_KE.Time, E1_KE.bpGy, 'Color', red, 'DisplayName', 'bpGy');
hold on;
plot(E1_VENT_DATA.Time, 3.5, 'v', 'MarkerEdgeColor', black, ...
    'MarkerFaceColor', goldenrod ); 
xline(E1_VENT_DATA.Time, ':', 'Color', goldenrod, 'DisplayName', '');
hold off;
xlabel('Time, local');
ylabel('degrees · s^-^1');
grid;
legend('bpGy','resp');

linkaxes([s1 s2], 'x');

%%
 
figCCAVG = figure('Color', white, 'Position', [0 500 1200 800] );

for z = 1:numCCRows

    nullTest = ccVxyz(:,:,z);

    if (nullTest == 0)  % ccVxyz is full of zeros, meaning invalid CC
        fprintf('Skipping cycle %d ...\n', z); 
    else
        fprintf('Processing cycle %d...\n', z); 
    end

    s3a = subplot(331);
    hold on;
    pCCAxyz = plot(cycleTime, ccAxyz(:,1,z), ':','Color', grey );
    hold off;        

    s3b = subplot(332);
    hold on;
    pCCVxyz = plot(cycleTime, ccVxyz(:,1,z), ':','Color', grey );
    hold off;

    s3c = subplot(333);
    hold on;
    pCCGxyz = plot(cycleTime, ccGxyz(:,1,z), ':','Color', grey );
    hold off;    

    s3d = subplot(334);
    hold on;
    pCCAxyz = plot(cycleTime, ccAxyz(:,2,z), ':','Color', grey );
    hold off;        

    s3e = subplot(335);
    hold on;
    pCCVxyz = plot(cycleTime, ccVxyz(:,2,z), ':','Color', grey );
    hold off;

    s3f = subplot(336);
    hold on;
    pCCGxyz = plot(cycleTime, ccGxyz(:,2,z), ':','Color', grey );
    hold off;    

    s3g = subplot(337);
    hold on;
    pCCAxyz = plot(cycleTime, ccAxyz(:,3,z), ':','Color', grey );
    hold off;        

    s3h = subplot(338);
    hold on;
    pCCVxyz = plot(cycleTime, ccVxyz(:,3,z), ':','Color', grey );
    hold off;

    s3i = subplot(339);
    hold on;
    pCCGxyz = plot(cycleTime, ccGxyz(:,3,z), ':','Color', grey );
    hold off;            

end

% overplot mean linear and mean rotational velocities

s3a = subplot(331);
hold on;
pMeanAx = plot(cycleTime, meanLinearAxyz(:,1), 'Color', ...
    blue, 'LineWidth', 3);    
hold off;

s3b = subplot(332);
hold on;
pMeanVx = plot(cycleTime, meanLinearVxyz(:,1), 'Color', ...
    blue, 'LineWidth', 3);
hold off;

s3c = subplot(333);
hold on;
pMeanGx = plot(cycleTime, meanRotationalGxyz(:,1), 'Color', ...
    blue, 'LineWidth', 3);    
hold off;

s3d = subplot(334);
hold on;
pMeanAy = plot(cycleTime, meanLinearAxyz(:,2), 'Color', ...
    red, 'LineWidth', 3);
hold off;

s3e = subplot(335);
hold on;
pMeanVy = plot(cycleTime, meanLinearVxyz(:,2), 'Color', ...
    red, 'LineWidth', 3);
hold off;

s3f = subplot(336);
hold on;
pMeanGy = plot(cycleTime, meanRotationalGxyz(:,2), 'Color', ...
    red, 'LineWidth', 3);    
hold off;

s3g = subplot(337);
hold on;
pMeanAz = plot(cycleTime, meanLinearAxyz(:,3), 'Color', ...
    goldenrod, 'LineWidth', 3);    
hold off;

s3h = subplot(338);
hold on;
pMeanVz = plot(cycleTime, meanLinearVxyz(:,3), 'Color', ...
    goldenrod, 'LineWidth', 3);
hold off;

s3i = subplot(339);
hold on;
pMeanGz = plot(cycleTime, meanRotationalGxyz(:,3), 'Color', ...
    goldenrod, 'LineWidth', 3);        
hold off;

% add in labels, titles and legends

s3a = subplot(331);
hold on;
xlabel('Time, seconds');
ylabel('m · s^-^2');
title('Cardiac cycle linear acceleration');
grid;
legend([ pMeanAx pCCAxyz(1,1) ], {'Ax','all Ax cycles'} );
hold off;    

s3b = subplot(332);
hold on;
xlabel('Time, seconds');
ylabel('m · s^-^1');
title('Cardiac cycle linear velocity');
grid;
legend([ pMeanVx pCCVxyz(1,1) ], {'Vx','all Vx cycles'} );
hold off;

s3c = subplot(333);
hold on;
xlabel('Time, seconds');
ylabel('degrees · s^-^1');
title('Cardiac cycle rotational velocity');
grid;
legend([ pMeanGx pCCGxyz(1,1) ], {'Gx','all Gx cycles'} );
hold off;

s3d = subplot(334);
hold on;
xlabel('Time, seconds');
ylabel('m · s^-^2');
title('Cardiac cycle linear acceleration');
grid;
legend([ pMeanAy pCCAxyz(1,1) ], {'Ay','all Ay cycles'} );
hold off;    

s3e = subplot(335);
hold on;
xlabel('Time, seconds');
ylabel('m · s^-^1');
title('Cardiac cycle linear velocity');
grid;
legend([ pMeanVy pCCGxyz(1,1) ], {'Vy','all Vy cycles'} );
hold off;

s3f = subplot(336);
hold on;
xlabel('Time, seconds');
ylabel('degrees · s^-^1');
title('Cardiac cycle rotational velocity');
grid;
legend([ pMeanGy pCCGxyz(1,1) ], {'Gy','all Gy cycles'} );
hold off;    

s3g = subplot(337);
hold on;
xlabel('Time, seconds');
ylabel('m · s^-^2');
title('Cardiac cycle linear acceleration');
grid;
legend([ pMeanAz pCCGxyz(1,1) ], {'Az','all Az cycles'} );
hold off;    

s3h = subplot(338);
hold on;
xlabel('Time, seconds');
ylabel('m · s^-^1');
title('Cardiac cycle linear velocity');
grid;
legend([ pMeanVz pCCGxyz(1,1) ], {'Vz','all Vz cycles'} );
hold off;

s3i = subplot(339);
hold on;
xlabel('Time, seconds');
ylabel('degrees · s^-^1');
title('Cardiac cycle rotational velocity');
grid;
legend([ pMeanGz pCCGxyz(1,1) ], {'Gz','all Gz cycles'} );
hold off;    

linkaxes([s3a s3b s3c s3d s3e s3f s3g s3h s3i],'x');


%%

[wsstBpOdba, fBpOdba] = wsst(E1_KINEMATICS.bpOdba, 100);
[wsstBpOdav, fBpOdav] = wsst(E1_KINEMATICS.bpOdav, 100);

figSpectra = figure('Color', white);

s1 = subplot(211);
pcolor(E1_KINEMATICS.Time, fBpOdba, abs(wsstBpOdba));
shading interp;
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Band-passed overall dynamic body acceleration (ODBA)');
ylim([0 5]);

s2 = subplot(212);
pcolor(E1_KINEMATICS.Time, fBpOdav, abs(wsstBpOdav));
shading interp;
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Band-passed overall dynamic angular velocity (ODAV)');
ylim([0 5]);

linkaxes([s1 s2],'x');

%%

[wsstBpAx, fBpAx] = wsst( sqrt(E1_KE.bpAx.^2), 100 );
[wsstBpGy, fBpGy] = wsst( sqrt(E1_KE.bpGy.^2), 100 );

[ridgeBpAx, iridgeBpAx] = wsstridge(wsstBpAx, 20, fBpAx);
[ridgeBpGy, iridgeBpGy] = wsstridge(wsstBpGy, 20, fBpGy);

figSpectra = figure('Color', white);

s1 = subplot(221);
pcolor(E1_KE.Time, fBpAx, abs(wsstBpAx));
shading interp;
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Band-passed sqrt(bpAx^2)');
ylim([0 5]);

s2 = subplot(222);
pcolor(E1_KE.Time, fBpGy, abs(wsstBpGy));
shading interp;
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Band-passed sqrt(bpGy^2)');
ylim([0 5]);

s3 = subplot(223);
plot(E1_KE.Time, ridgeBpAx, 'Color', blue);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Band-passed sqrt(bpAx^2)');
ylim([0 5]);

s4 = subplot(224);
plot(E1_KE.Time, ridgeBpGy, 'Color', red);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Band-passed sqrt(bpGy^2)');
ylim([0 5]);


%%

figure('Color', white);

s9a = subplot(221);
plot(E1_KINEMATICS.bpOdba, 'Color', red );
grid;

s9b = subplot(222);
plot(E1_KINEMATICS.bpOdav, 'Color', green );
grid;

s9c = subplot(223);
plot(E1_KINEMATICS.sumImfBpOdba, 'Color', maroon);
hold on;
plot(E1_KINEMATICS.bestImfBpOdba,'Color', black, 'LineWidth', 2);
hold off;
grid;

s9d = subplot(224);
plot(E1_KINEMATICS.sumImfBpOdav, 'Color', purple );
hold on;
plot(E1_KINEMATICS.bestImfBpOdav,'Color', black, 'LineWidth', 2);
hold off;
grid;

linkaxes([s9a s9b s9c s9d],'x');

%% 

figure('Color', white);
s10a = subplot(211);
plot(E1_KE.time_s, E1_KE.bpAx, 'Color', blue);
xlabel('Time, seconds');
ylabel('m · s^-^2');
title('Band pass-filtered accelerometer, x-axis');
grid;
s10b = subplot(212);
plot(E1_KE.time_s, E1_KE.bpGy, 'Color', red);
xlabel('Time, seconds');
ylabel('degrees · s^-^1');
title('Band pass-filtered gyroscope, y-axis');
grid;
linkaxes([s10a s10b],'x');


%%

% Hoku t-tests for all epochs

[h,p,ci,stats] = ttest2(...
    CC_HR_DATA.iHr(CC_HR_DATA.altId == 1 & CC_HR_DATA.epoch == 1), ...
    CC_HR_DATA.ensWsstHr(CC_HR_DATA.altId == 1 & CC_HR_DATA.epoch == 1), ...
    'VarType', 'unequal')

[h,p,ci,stats] = ttest2(...
    CC_HR_DATA.iHr(CC_HR_DATA.altId == 1 & CC_HR_DATA.epoch == 3), ...
    CC_HR_DATA.ensWsstHr(CC_HR_DATA.altId == 1 & CC_HR_DATA.epoch == 3), ...
    'VarType', 'unequal' )

% Liho

[h,p,ci,stats] = ttest2(...
    CC_HR_DATA.iHr(CC_HR_DATA.altId == 2 & CC_HR_DATA.epoch == 1), ...
    CC_HR_DATA.ensWsstHr(CC_HR_DATA.altId == 2 & CC_HR_DATA.epoch == 1), ...
    'VarType', 'unequal')

[h,p,ci,stats] = ttest2(...
    CC_HR_DATA.iHr(CC_HR_DATA.altId == 2 & CC_HR_DATA.epoch == 3), ...
    CC_HR_DATA.ensWsstHr(CC_HR_DATA.altId == 2 & CC_HR_DATA.epoch == 3), ...
    'VarType', 'unequal')

[h,p,ci,stats] = ttest2(...
    CC_HR_DATA.iHr(CC_HR_DATA.altId == 2 & CC_HR_DATA.epoch == 5), ...
    CC_HR_DATA.ensWsstHr(CC_HR_DATA.altId == 2 & CC_HR_DATA.epoch == 5), ...
    'VarType', 'unequal')

% Kolohe

[h,p,ci,stats] = ttest2(...
    CC_HR_DATA.iHr(CC_HR_DATA.altId == 3 & CC_HR_DATA.epoch == 1), ...
    CC_HR_DATA.ensWsstHr(CC_HR_DATA.altId == 3 & CC_HR_DATA.epoch == 1), ...
    'VarType', 'unequal')

[h,p,ci,stats] = ttest2(...
    CC_HR_DATA.iHr(CC_HR_DATA.altId == 3 & CC_HR_DATA.epoch == 3), ...
    CC_HR_DATA.ensWsstHr(CC_HR_DATA.altId == 3 & CC_HR_DATA.epoch == 3), ...
    'VarType', 'unequal')

[h,p,ci,stats] = ttest2(...
    CC_HR_DATA.iHr(CC_HR_DATA.altId == 3 & CC_HR_DATA.epoch == 5), ...
    CC_HR_DATA.ensWsstHr(CC_HR_DATA.altId == 3 & CC_HR_DATA.epoch == 5), ...
    'VarType', 'unequal')

% Hua


[h,p,ci,stats] = ttest2(...
    CC_HR_DATA.iHr(CC_HR_DATA.altId == 4 & CC_HR_DATA.epoch == 1), ...
    CC_HR_DATA.ensWsstHr(CC_HR_DATA.altId == 4 & CC_HR_DATA.epoch == 1), ...
    'VarType', 'unequal')

[h,p,ci,stats] = ttest2(...
    CC_HR_DATA.iHr(CC_HR_DATA.altId == 4 & CC_HR_DATA.epoch == 3), ...
    CC_HR_DATA.ensWsstHr(CC_HR_DATA.altId == 4 & CC_HR_DATA.epoch == 3), ...
    'VarType', 'unequal')

[h,p,ci,stats] = ttest2(...
    CC_HR_DATA.iHr(CC_HR_DATA.altId == 4 & CC_HR_DATA.epoch == 5), ...
    CC_HR_DATA.ensWsstHr(CC_HR_DATA.altId == 4 & CC_HR_DATA.epoch == 5), ...
    'VarType', 'unequal')


% Noa


[h,p,ci,stats] = ttest2(...
    CC_HR_DATA.iHr(CC_HR_DATA.altId == 5 & CC_HR_DATA.epoch == 1), ...
    CC_HR_DATA.ensWsstHr(CC_HR_DATA.altId == 5 & CC_HR_DATA.epoch == 1), ...
    'VarType', 'unequal')

[h,p,ci,stats] = ttest2(...
    CC_HR_DATA.iHr(CC_HR_DATA.altId == 5 & CC_HR_DATA.epoch == 3), ...
    CC_HR_DATA.ensWsstHr(CC_HR_DATA.altId == 5 & CC_HR_DATA.epoch == 3), ...
    'VarType', 'unequal')

[h,p,ci,stats] = ttest2(...
    CC_HR_DATA.iHr(CC_HR_DATA.altId == 5 & CC_HR_DATA.epoch == 5), ...
    CC_HR_DATA.ensWsstHr(CC_HR_DATA.altId == 5 & CC_HR_DATA.epoch == 5), ...
    'VarType', 'unequal')

% Lono

[h,p,ci,stats] = ttest2(...
    CC_HR_DATA.iHr(CC_HR_DATA.altId == 6 & CC_HR_DATA.epoch == 1), ...
    CC_HR_DATA.ensWsstHr(CC_HR_DATA.altId == 6 & CC_HR_DATA.epoch == 1), ...
    'VarType', 'unequal')

[h,p,ci,stats] = ttest2(...
    CC_HR_DATA.iHr(CC_HR_DATA.altId == 6 & CC_HR_DATA.epoch == 3), ...
    CC_HR_DATA.ensWsstHr(CC_HR_DATA.altId == 6 & CC_HR_DATA.epoch == 3), ...
    'VarType', 'unequal')

%%

figure; 
plot(E1_KINEMATICS.Time, E1_KINEMATICS.bestImfBpOdba, 'LineWidth', 2); 
hold on; 
xline(E1_VENT_DATA.Time, '--', 'Color', goldenrod, 'LineWidth', 2); 
plot(E3_KINEMATICS.Time, E3_KINEMATICS.bestImfBpOdba, 'LineWidth', 2);
xline(E3_VENT_DATA.Time, '--', 'Color', goldenrod, 'LineWidth', 2);
plot(E5_KINEMATICS.Time, E5_KINEMATICS.bestImfBpOdba, 'LineWidth', 2); 
xline(E5_VENT_DATA.Time, '--', 'Color', goldenrod, 'LineWidth', 2);
hold off;
grid;
