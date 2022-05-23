%%  ftCardiacEnergyPlots.m

%   Written by Dave Haas between 9 July 2 October 2021

clc;
clear;

global TAG_PATHS;

% named colors
blue = [0 0.4470 0.7410];           % '#0072BD';
red = [0.8500 0.3250 0.0980];       % '#D95319'	
goldenrod = [0.9290 0.6940 0.1250];	% '#EDB120'	
purple = [0.4940 0.1840 0.5560];    % '#7E2F8E'	
green = [0.4660 0.6740 0.1880];     % '#77AC30'	
cyan = [0.3010 0.7450 0.9330];      % '#4DBEEE'	
maroon = [0.6350 0.0780 0.1840];    % '#A2142F'
grey = [0.6 0.6 0.6];               % '#999999'
black = [0 0 0];                    % '#000000
white = [1 1 1];                    % '#ffffff'

% set some default LED* colors...

led1Color = black;
led2Color = blue;
led3Color = red;
led4Color = goldenrod;

% ... and some default kinematic colors...

kxColor = blue;
kyColor = red;
kzColor = goldenrod;

% define makeFilterPlots if you want to see them during early steps...

makeFilterPlots = false;

%% select a tag for analysis

tagStr = input('Enter a tag for analysis (e.g.: tt21_141a): ','s');

fullRawFileName = sprintf('%s/%sraw.mat', TAG_PATHS.RAW, ...
    tagStr);
if ( exist(fullRawFileName, 'file') )
    fprintf('Raw file located for this trial... loading it!\n');
    load(fullRawFileName);
    
    % confirm presence of OPTICS variable
    if ( exist('OPTICS','var') )
        
        fprintf('Key variables are present and loaded... proceeding!\n');

    else
        
        fprintf('Key variables are not present. Check RAW file and restart.\n');
        return;
        
    end
    
else
    
    fprintf('Raw file is not present. Run ftPostProcessing.m to create it.\n');
    return;

end



%% ftCardiacEnergyPlots.m
%
%   Written by Dave Haas
%   18 September 2021

figure('Color', white, 'Position', [0 400 800 500] );

sCCa = subplot(211);
plot(E1_KINEMATICS.Time, E1_KINEMATICS.hrTimeConsensus,'Color',green);
hold on;
plot(E3_KINEMATICS.Time, E3_KINEMATICS.hrTimeConsensus, 'Color', blue);
plot(E5_KINEMATICS.Time, E5_KINEMATICS.hrTimeConsensus, 'Color', red);
plot(E1_CC_ENERGY.Time, E1_CC_ENERGY.iHR, 'o', 'Color', green);
plot(E3_CC_ENERGY.Time, E3_CC_ENERGY.iHR, 'o', 'Color', blue);
plot(E5_CC_ENERGY.Time, E5_CC_ENERGY.iHR, 'o', 'Color', red);
xline(E1_VENT_DATA.Time, ':', 'Color', green, 'LineWidth', 2);
xline(E5_VENT_DATA.Time, ':', 'Color', red, 'LineWidth', 2);
hold off;
xlabel('Time, local');
ylabel('Heart rate,  beats · min^-^1');
titleTxt = sprintf('Heart rate profile for %s', TRIAL.SUBJECT.id);
title(titleTxt);
ylim([0 150]);
grid; 
legend('baseline wsstHR','apnea wsstHR','recovery wsstHR', ...
    'baseline iHR','apnea iHR','recovery iHR', ...
    'baseline vents', 'recovery vents');

sCCb = subplot(212);
plot(E1_CC_ENERGY.Time, E1_CC_ENERGY.iTotalKE, 'x', 'Color', green);
hold on;
plot(E3_CC_ENERGY.Time, E3_CC_ENERGY.iTotalKE, 'x', 'Color', blue);
plot(E5_CC_ENERGY.Time, E5_CC_ENERGY.iTotalKE, 'x', 'Color', red);
xline(E1_VENT_DATA.Time, ':', 'Color', green, 'LineWidth', 2);
xline(E5_VENT_DATA.Time, ':', 'Color', red, 'LineWidth', 2);
hold on;
xlabel('Time, local');
ylabel('Kinetic energy,  J · s');
titleTxt = sprintf('Cardiac energy profile for %s', TRIAL.SUBJECT.id);
title(titleTxt);
grid;
legend('baseline','apnea','recovery', ...
    'baseline vents','recovery vents');

linkaxes([sCCa sCCb],'x');

%%

figure('Color', white, 'Position', [0 600 1200 500] );
plot(E1_CC_ENERGY.Time, E1_CC_ENERGY.iTotalKE, 'x', 'Color', green);
hold on;
plot(E3_CC_ENERGY.Time, E3_CC_ENERGY.iTotalKE, 'x', 'Color', blue);
plot(E5_CC_ENERGY.Time, E5_CC_ENERGY.iTotalKE, 'x', 'Color', red);
xline(E1_VENT_DATA.Time, '-', 'Color', green);
xline(E5_VENT_DATA.Time, '-', 'Color', red);
hold off;
xlabel('Time, seconds');
ylabel('J · s');
titleTxt = 'Change in cardiac energy during pre-apnea (baseline), apnea, and post-apnea (recovery)';
title(titleTxt);
grid;
legend('baseline','apnea','recovery','baseline vents','recovery vents');