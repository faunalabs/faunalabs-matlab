%% ftChapterThreeStats.m
%
%   Dave Haas
%   7 October 2021
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

%% Load the trial timetables for everything

load('/Users/dave/Documents/MATLAB/tursiopsTimetables.mat')


%% Make plots of all the iKE stuff by experimental trial for all animals


for whichTrial = 1:6

    switch(whichTrial)
        
        case 1
            trialStr = 'tt21_128f';
            animalStr = string(unique(CC_HR_DATA.id(CC_HR_DATA.altId == 1)));
        case 2
            trialStr = 'tt21_134a';
            animalStr = string(unique(CC_HR_DATA.id(CC_HR_DATA.altId == 2)));
        case 3
            trialStr = 'tt21_134b';
            animalStr = string(unique(CC_HR_DATA.id(CC_HR_DATA.altId == 3)));
        case 4
            trialStr = 'tt21_141e';
            animalStr = string(unique(CC_HR_DATA.id(CC_HR_DATA.altId == 4)));
        case 5
            trialStr = 'tt21_142c';
            animalStr = string(unique(CC_HR_DATA.id(CC_HR_DATA.altId == 5)));
        case 6
            trialStr = 'tt21_142d';
            animalStr = string(unique(CC_HR_DATA.id(CC_HR_DATA.altId == 6)));
    end
    
    
%     colorsLinKE   = [ blue; maroon; goldenrod; ];
%     colorsRotKE   = [ blue; green; goldenrod; ];
    colorsTotalKE = [ blue; red; goldenrod; ];
    
%     figure('Color', white, 'Position', [0 0 1600 1200]);
    figure('Color', white, 'Position', [0 0 1600 400]);
    
%     s1a = subplot(311);
%     ax1 = gca;
%     colororder(ax1, colorsLinKE);
%     yyaxis left;
%     plot(CC_HR_DATA.Time(CC_HR_DATA.altId == whichTrial), ...
%         CC_HR_DATA.iHr(CC_HR_DATA.altId == whichTrial), '+');
%     xlabel('Time, local');
%     ylabel('beats · min^-^1');
%     titleStr = sprintf('Trial %s, Animal %s - linear kinetic energy', ...
%         trialStr, animalStr);
%     title(titleStr, 'Interpreter','none');
%     yyaxis right;
%     bar(CC_HR_DATA.Time(CC_HR_DATA.altId == whichTrial), ...
%         hampel( CC_HR_DATA.iLinKE(CC_HR_DATA.altId == whichTrial), 3, 5) );
%     ylabel('J · s^-^1');
%     grid;
%     hold on;
%     xline(VENT_HR_DATA.Time(VENT_HR_DATA.altId == whichTrial), ...
%         'Color', goldenrod, 'LineWidth', 2);
%     hold off;
%     
%     legend('ifH','linKE','vent');
%     
% 
%     
%     s1b = subplot(312);
%     ax2 = gca;
%     colororder(ax2, colorsRotKE);
%     yyaxis left;
%     plot(CC_HR_DATA.Time(CC_HR_DATA.altId == whichTrial), ...
%         CC_HR_DATA.iHr(CC_HR_DATA.altId == whichTrial), '+');
%     xlabel('Time, local');
%     ylabel('beats · min^-^1');
%     titleStr = sprintf('Trial %s, Animal %s - rotational kinetic energy', ...
%         trialStr, animalStr);
%     title(titleStr, 'Interpreter','none');
%     yyaxis right;
%     bar(CC_HR_DATA.Time(CC_HR_DATA.altId == whichTrial), ...
%         hampel( CC_HR_DATA.iRotKE(CC_HR_DATA.altId == whichTrial), 3, 5) );
%     ylabel('J · s^-^1');
%     grid;
%     
%     hold on;
%     xline(VENT_HR_DATA.Time(VENT_HR_DATA.altId == whichTrial), ...
%         'Color', goldenrod, 'LineWidth', 2);
%     hold off;    
%     legend('ifH','iRotKE','vent');
%     
%     s1c = subplot(313);

    ax3 = gca;
    colororder(ax3,colorsTotalKE);
 
    yyaxis right;
    bar(CC_HR_DATA.Time(CC_HR_DATA.altId == whichTrial), ...
        hampel( CC_HR_DATA.iTotalKE(CC_HR_DATA.altId == whichTrial), 3, 5) );
    ylabel('J · s^-^1', 'FontSize', 15);
    grid;
    
    hold on;
    xline(VENT_HR_DATA.Time(VENT_HR_DATA.altId == whichTrial), ...
        'Color', goldenrod, 'LineWidth', 2);
    yyaxis left;
    plot(CC_HR_DATA.Time(CC_HR_DATA.altId == whichTrial), ...
        CC_HR_DATA.iHr(CC_HR_DATA.altId == whichTrial), '+');
    xlabel('Time, local');
    ylabel('beats · min^-^1', 'FontSize', 15);
    titleStr = sprintf('Trial %s, Animal %s - total kinetic energy', ...
        trialStr, animalStr);
    title(titleStr, 'Interpreter','none');    
    
    hold off;
    legend('iKE','vent', 'ifH');
    
%     linkaxes([s1a s1b s1c],'x');
    
    
end


%%

figure('Color', white, 'Position', [0 0 1200 400] );

s1 = subplot(131);
boxplot( ...
    hampel(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1), 3, 5), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 1), ...
    'Labels', {'63HF','01L5','6JK5','83H1','9ON6','9FL3'}, 'Notch','on' );
xlabel('Animal ID');
ylabel('J · s^-^1');
title('Epoch 1 · cardiac cycle instantaneous kinetic energy');
grid;

s2 = subplot(132);
boxplot( ...
    hampel(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 3), 3, 5), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 3), ...
    'Labels', {'63HF','01L5','6JK5','83H1','9ON6','9FL3'}, 'Notch','on' );
xlabel('Animal ID');
ylabel('J · s^-^1');
title('Epoch 3 · cardiac cycle instantaneous kinetic energy');
grid;


s3 = subplot(133);
boxplot( ...
    hampel(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5), 3, 5), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 5), ...
    'Labels', {'01L5','6JK5','83H1','9ON6'}, 'Notch','on' );
xlabel('Animal ID');
ylabel('J · s^-^1');
title('Epoch 5 · cardiac cycle instantaneous kinetic energy');
grid;


%%


figure('Color',[1 1 1],'Position',[0 0 1200 1200]);

subplot(331);
boxplot( ...
    hampel( CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1), 3, 5), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 1), 'Notch','on');
xlabel('Animal ID');
ylabel('J · s^-^1');
title('Baseline Epoch (E1): iTotalKE');
grid;

subplot(332);
boxplot( ...
    hampel( CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 3), 3, 5), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 3), 'Notch','on');
xlabel('Animal ID');
ylabel('J · s^-^1');
title('Apnea Epoch (E3): iTotalKE');
grid;

subplot(333);
boxplot( ...
    hampel( CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5), 3, 5), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 5), 'Notch','on');
xlabel('Animal ID');
ylabel('J · s^-^1');
title('Recovery Epoch (E5): iTotalKE');
grid;

%

subplot(334);
boxplot( ...
    hampel( CC_HR_DATA.iLinKE(CC_HR_DATA.epoch == 1), 3, 5), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 1), 'Notch','on');
xlabel('Animal ID');
ylabel('J · s^-^1');
title('Baseline Epoch (E1): iLinKE');
grid;

subplot(335);
boxplot(...
    hampel( CC_HR_DATA.iLinKE(CC_HR_DATA.epoch == 3), 3, 5), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 3), 'Notch','on');
xlabel('Animal ID');
ylabel('J · s^-^1');
title('Apnea Epoch (E3): iLinKE');
grid;

subplot(336);
boxplot( ...
    hampel( CC_HR_DATA.iLinKE(CC_HR_DATA.epoch == 5), 3, 5), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 5), 'Notch','on');
xlabel('Animal ID');
ylabel('J · s^-^1');
title('Recovery Epoch (E5): iLinKE');
grid;

%

subplot(337);
boxplot( ...
    hampel( CC_HR_DATA.iRotKE(CC_HR_DATA.epoch == 1), 3, 5), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 1), 'Notch','on');
xlabel('Animal ID');
ylabel('J · s^-^1');
title('Baseline Epoch (E1): iRotKE');
grid;

subplot(338);
boxplot( ...
    hampel( CC_HR_DATA.iRotKE(CC_HR_DATA.epoch == 3), 3, 5), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 3), 'Notch','on');
xlabel('Animal ID');
ylabel('J · s^-^1');
title('Apnea Epoch (E3): iRotKE');
grid;

subplot(339);
boxplot( ...
    hampel( CC_HR_DATA.iRotKE(CC_HR_DATA.epoch == 5), 3, 5), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 5), 'Notch','on');
xlabel('Animal ID');
ylabel('J · s^-^1');
title('Recovery Epoch (E5): iRotKE');
grid;


%%

for whichTrial = 1:6

    fprintf('Working on trial %d...\n', whichTrial);
    
    switch(whichTrial)
        case 1
            load('/Users/dave/Desktop/DTAG/raw/tt21_128fraw.mat');
        case 2
            load('/Users/dave/Desktop/DTAG/raw/tt21_134araw.mat');
        case 3
            load('/Users/dave/Desktop/DTAG/raw/tt21_134braw.mat');
        case 4
            load('/Users/dave/Desktop/DTAG/raw/tt21_141eraw.mat');
        case 5
            load('/Users/dave/Desktop/DTAG/raw/tt21_142craw.mat');
        case 6
            load('/Users/dave/Desktop/DTAG/raw/tt21_142draw.mat');
    end
    
for whichEpoch = 1:2:5
    
    fprintf('\tProcessing epoch %d...\n', whichEpoch);
    
    switch(whichEpoch)
    
        case 1
            thisKE = E1_KE;
            aSize = height(thisKE);
            bpAx = thisKE.bpAx;
            bpAy = thisKE.bpAy;
            bpAz = thisKE.bpAz;
            bpGxSI = deg2rad(thisKE.bpGx);
            bpGySI = deg2rad(thisKE.bpGy);
            bpGzSI = deg2rad(thisKE.bpGz);            
        case 3
            thisKE = E3_KE;
            aSize = height(thisKE);
            bpAx = thisKE.bpAx;
            bpAy = thisKE.bpAy;
            bpAz = thisKE.bpAz;            
            bpGxSI = deg2rad(thisKE.bpGx);
            bpGySI = deg2rad(thisKE.bpGy);
            bpGzSI = deg2rad(thisKE.bpGz);
        case 5
            if (whichTrial ~= 1 & whichTrial ~= 6)
                thisKE = E5_KE;
                aSize = height(thisKE);
                bpAx = thisKE.bpAx;
                bpAy = thisKE.bpAy;
                bpAz = thisKE.bpAz;
                bpGxSI = deg2rad(thisKE.bpGx);
                bpGySI = deg2rad(thisKE.bpGy);
                bpGzSI = deg2rad(thisKE.bpGz);                
            else
                aSize = 1000;
                bpAx = rand([1000 1]);
                bpAy = rand([1000 1]);
                bpAz = rand([1000 1]);
                bpGxSI = rand([1000 1]);
                bpGySI = rand([1000 1]);
                bpGzSI = rand([1000 1]);
            end
    end
    
    
    if (exist('thisA','var'))
        clear thisA;
    end
    
    thisA(1:aSize,7) = zeros;

    % observations are standardized to have a unit variance

    thisA(1:aSize, 1) = bpAx / std(bpAx);
    thisA(1:aSize, 2) = bpAy / std(bpAy);
    thisA(1:aSize, 3) = bpAz / std(bpAz);
    % thisA(1:aSize, 4) = thisKE.bpGx / std(thisKE.bpGx);
    % thisA(1:aSize, 5) = thisKE.bpGy / std(thisKE.bpGy);
    % thisA(1:aSize, 6) = thisKE.bpGz / std(thisKE.bpGz);
    thisA(1:aSize, 4) = bpGxSI / std(bpGxSI);
    thisA(1:aSize, 5) = bpGySI / std(bpGySI);
    thisA(1:aSize, 6) = bpGzSI / std(bpGzSI);
    thisA(1:aSize, 7) = -1;

    % get the V using svd on thisA

    [thisU, thisS, thisV] = svd(thisA);

    % take the first six coefficients of coeffsV for use in KE equation...

    coeffsV = thisV(1:6,7);

    % confused about what to do next regarding "put V into the kinetic energy
    % engine to compute Ix Iy Iz...

    massTag = 0.403;

    % gyroKEx = 0.5 * massTest * ( (coeffsV(4) * thisKE.bpGx ) .^ 2);
    % gyroKEy = 0.5 * massTest * ( (coeffsV(5) * thisKE.bpGy ) .^ 2);
    % gyroKEz = 0.5 * massTest * ( (coeffsV(6) * thisKE.bpGz ) .^ 2);
    gyroKEx = 0.5 * massTag * ( (coeffsV(4) * bpGxSI ) .^ 2);
    gyroKEy = 0.5 * massTag * ( (coeffsV(5) * bpGySI ) .^ 2);
    gyroKEz = 0.5 * massTag * ( (coeffsV(6) * bpGzSI ) .^ 2);

    % convert each x, y, z into Ix, Iy, Iz using the following logic:

    % gyroKEx = 0.5 * Ix * bpGx .^ 2
    % gyroKEx / 0.5 = Ix * bpGx .^ 2
    % gyroKEx / 0.5 * bpGx .^ 2 = Ix

    Ix = gyroKEx ./ (0.5 * bpGxSI .^ 2 ); 
    Iy = gyroKEy ./ (0.5 * bpGySI .^ 2 ); 
    Iz = gyroKEz ./ (0.5 * bpGzSI .^ 2 ); 

    I = [Ix, Iy, Iz];

    if (exist('Axyz','var'))
       clear Axyz; 
    end

    if (exist('Gxyz','var'))
       clear Gxyz; 
    end

    if (exist('Vxyz','var'))
       clear Vxyz; 
    end
    
    Axyz = [bpAx, bpAy, bpAz];
    Gxyz = [bpGxSI, bpGySI, bpGzSI];
    
    for v = 1:height(Axyz)-1
        Vxyz(v,:) = simps(Axyz(v:v+1,:));
    end
    
    thisLinKE = 0.5 * massTag .* (Vxyz .^ 2);
    thisRotKE = 0.5 * (I .* Gxyz .^ 2 );
    
    sizeLinKE = height(thisLinKE);
    sizeRotKE = height(thisRotKE);
    
    sizeDiff = diff([sizeRotKE, sizeLinKE]);
    
    switch(sizeDiff)
        case -1
            fprintf('\tSizes of linKE < rotKE. Adjusting!\n');            
            thisRotKE = thisRotKE(1:end-1,1:3);
        case -2
            fprintf('\tSizes of linKE < rotKE - off by 2. Adjusting!\n');            
            thisRotKE = thisRotKE(1:end-2,1:3);            
        case 1
           fprintf('\tSizes of linKE > rotKE. Adjusting!\n');            
            thisLinKE = thisLinKE(1:end-1,1:3);            
        otherwise
            fprintf('\tSizes of KE match. Doing nothing!\n');
    end
    
    thisTotKE = thisLinKE + thisRotKE;
    ensTotKE = (thisLinKE + thisRotKE) ./ 2;
    
    figName = sprintf('Trial %d | Epoch %d', whichTrial, whichEpoch);
    
    figure('Color', white', 'Name', figName);
    
    s1a = subplot(411);
	plot(thisLinKE);
    ylabel('J · s^-^1');
    title('Linear kinetic energy');
    grid; 
    
    s1b = subplot(412);
    plot(thisRotKE);
    ylabel('J · s^-^1');
    title('Rotational kinetic energy');
    grid; 

    s1c = subplot(413);
    plot(thisTotKE);
    ylabel('J · s^-^1');
    title('Total kinetic energy');
    grid; 

    s1d = subplot(414);
    plot(ensTotKE);
    ylabel('J · s^-^1');
    title('Ensemble total kinetic energy');
    grid;     
    
    linkaxes([s1a s1b s1c s1d],'x');
    
    fprintf('\tIx: %4.5f | Iy: %4.5f | Iz: %4.5f \n', Ix(1), Iy(1), Iz(1) );    
    
    switch(whichEpoch)
    
        case 1
            I_E1 = mean(I);
        case 3
            I_E3 = mean(I);
        case 5
            I_E5 = mean(I);
    end    

end

fprintf('\tE1 --> Ix: %2.6f| Iy: %2.6f | Iz: %2.6f \n', I_E1(1), I_E1(2), I_E1(3) );
fprintf('\tE3 --> Ix: %2.6f| Iy: %2.6f | Iz: %2.6f \n', I_E3(1), I_E3(2), I_E3(3) );
fprintf('\tE5 --> Ix: %2.6f| Iy: %2.6f | Iz: %2.6f \n', I_E5(1), I_E5(2), I_E5(3) );

MOI(whichTrial).I = [I_E1; I_E3; I_E5];

end

(MOI(1:6).I)

%% moment of inertia calculator for each dolphin

% set the tag mass for correction factor

tagMass = 0.403;        % 403 grams

% estimates from Fish et al 2006; units kg · m^2; these are for spinner
% dolphins will be be rounding errors compared to the 

Iflippers = 0.069;
Ifin      = 0.043;
Iflukes   = 0.007;

% Individual animal calculations for moments of inertia
% all units converted to meters and kg

%
% ---=== HOKU ===---
%


HOKU        = struct;
HOKU.mass   = 188;
HOKU.length = 2.54;         % 254 cm
HOKU.girth  = 1.3081;       % 130.81 cm

% formula for Ibody from Fish et al 2006
% 2/5 * (65 * 0.15 ^ 2) = 0.585, confirmed from paper

HOKU.Ibody  = 2/5 * ( HOKU.mass * (HOKU.girth / (2 * pi) ) ^ 2 ); 
HOKU.Itotal = HOKU.Ibody + Iflippers + Ifin + Iflukes

%
% ---=== LIHO ===---
%

LIHO        = struct;
LIHO.mass   = 165;
LIHO.length = 2.4003;         % 254 cm
LIHO.girth  = 1.2426;       % 130.81 cm

% formula for Ibody from Fish et al 2006
% 2/5 * (65 * 0.15 ^ 2) = 0.585, confirmed from paper

LIHO.Ibody  = 2/5 * ( LIHO.mass * (LIHO.girth / (2 * pi) ) ^ 2 ); 
LIHO.Itotal = LIHO.Ibody + Iflippers + Ifin + Iflukes

%
% ---=== KOLOHE ===---
%

KOLOHE        = struct;
KOLOHE.mass   = 211;
KOLOHE.length = 2.5908;         % 254 cm
KOLOHE.girth  = 1.4478;       % 130.81 cm

% formula for Ibody from Fish et al 2006
% 2/5 * (65 * 0.15 ^ 2) = 0.585, confirmed from paper

KOLOHE.Ibody  = 2/5 * ( KOLOHE.mass * (KOLOHE.girth / (2 * pi) ) ^ 2 ); 
KOLOHE.Itotal = KOLOHE.Ibody + Iflippers + Ifin + Iflukes

%
% ---=== HUA ===---
%

HUA        = struct;
HUA.mass   = 154;
HUA.length = 2.4130;         % 254 cm
HUA.girth  = 1.2573;       % 130.81 cm

% formula for Ibody from Fish et al 2006
% 2/5 * (65 * 0.15 ^ 2) = 0.585, confirmed from paper

HUA.Ibody  = 2/5 * ( HUA.mass * (HUA.girth / (2 * pi) ) ^ 2 ); 
HUA.Itotal = HUA.Ibody + Iflippers + Ifin + Iflukes

%
% ---=== NOA ===---
%

NOA        = struct;
NOA.mass   = 187;
NOA.length = 2.4765;         % 254 cm
NOA.girth  = 1.3716;       % 130.81 cm

% formula for Ibody from Fish et al 2006
% 2/5 * (65 * 0.15 ^ 2) = 0.585, confirmed from paper

NOA.Ibody  = 2/5 * ( NOA.mass * (NOA.girth / (2 * pi) ) ^ 2 ); 
NOA.Itotal = NOA.Ibody + Iflippers + Ifin + Iflukes


%
% ---=== LONO ===---
%

LONO        = struct;
LONO.mass   = 228;
LONO.length = 2.7432;         % 254 cm
LONO.girth  = 1.4542;       % 130.81 cm

% formula for Ibody from Fish et al 2006
% 2/5 * (65 * 0.15 ^ 2) = 0.585, confirmed from paper

LONO.Ibody  = 2/5 * ( LONO.mass * (LONO.girth / (2 * pi) ) ^ 2 ); 
LONO.Itotal = LONO.Ibody + Iflippers + Ifin + Iflukes

%
% ---=== convert degrees · sec^-1 to rads · sec^-1=== ---
%
% 
% bpGxSI = deg2rad(bpGx);
% bpGySI = deg2rad(bpGy);
% bpGzSI = deg2rad(bpGz);

%% 

for whichEpoch = 1:2:5
   
    fprintf('Working on Epoch: %d\n', whichEpoch); 
    
    switch(whichEpoch)
        case 1
            thisKE = E1_KE;        
        case 3
            thisKE = E3_KE;
        case 5
            thisKE = E5_KE;
    end
    
    Axyz = [ thisKE.bpAx thisKE.bpAy thisKE.bpAz ];
    Gxyz = [ thisKE.bpGx thisKE.bpGy thisKE.bpGz ];

    Gxyz_SI = Gxyz .* pi/180;

    % use size(Axyz) to make Vxyz, the dt integral of Axyz

    clear Vxyz;
    sizeVxyz = size(Axyz);
    Vxyz(sizeVxyz(1),sizeVxyz(2)) = zeros;    

    % Simpson's rule for numerical estimation needs even numbers

    if ( mod(length(Vxyz),2) == 0 )            % even
        sampleRange = 1:length(Vxyz);
    else                                                % odd
        sampleRange = 1:length(Vxyz) - 1;
    end

    % compute instantaneous velocities using Simpson's rule

    for v = 1:numel(sampleRange)-1

        Vxyz(v,:) = simps(Axyz(v:v+1,:));

    end

    lenVxyz = length(Vxyz);

    thisLinKE = 0.5 * ( KOLOHE.mass ) * Vxyz.^2 ;
    totalLinKE = sum(thisLinKE,2);

    thisRotKE = 0.5 * ( KOLOHE.mass ) * Gxyz_SI.^2 ;
    totalRotKE = sum(thisRotKE,2);
    
    thisTotalKE = thisLinKE + thisRotKE;
    totalKE = sum(thisTotalKE,2);

    switch(whichEpoch)
        case 1
            e1Time = E1_KE.Time;
            e1LinKE = totalLinKE;
            e1RotKE = totalRotKE;
            e1TotKE = totalKE;
        case 3
            e3Time = E3_KE.Time;
            e3LinKE = totalLinKE;
            e3RotKE = totalRotKE;
            e3TotKE = totalKE;            
        case 5
            e5Time = E5_KE.Time;
            e5LinKE = totalLinKE;
            e5RotKE = totalRotKE;
            e5TotKE = totalKE;            
    end
    
end


figure; 
s1 = subplot(311);
plot(e1Time, e1LinKE);
hold on;
plot(e3Time, e3LinKE);
plot(e5Time, e5LinKE);
hold off;
xlabel('Time, local');
ylabel('J · s');
title('Linear kinetic energy');
grid;
s2 = subplot(312);
plot(e1Time, e1RotKE);
hold on;
plot(e3Time, e3RotKE);
plot(e5Time, e5RotKE);
hold off;
xlabel('Time, local');
ylabel('J · s');
title('Rotational kinetic energy');
grid;
s3 = subplot(313);
plot(e1Time, e1TotKE);
hold on;
plot(e3Time, e3TotKE);
plot(e5Time, e5TotKE);
hold off;
xlabel('Time, local');
ylabel('J · s');
title('Total kinetic energy');
grid;
linkaxes([s1 s2 s3],'x');

%% Run some Kolmogorov-Smirnov



%% Run mixed-effects ANOVA: ensWsstHr ~ fixed(epoch) + random(id)

% a two-way ANOVA with unbalanced design to compare heart rates during
% baseline, apnea, and recovery epochs

% mixed-effects (two-way) ANOVA: ensWsstHr ~ fixed(epoch) + random(id)

y1a = CC_HR_DATA.iTotalKE;
yfH = [y1a];
g1a = CC_HR_DATA.epoch;
g2a = CC_HR_DATA.altId;

[pAOV1a, tblAOV1a, statsAOV1a, termsAOV1a] = ...
    anovan( y1a, {g1a, g2a}, ...
    'model', 1, 'random', 2, ... %'interaction', ...
    'varnames', {'Epoch','AnimalID'} );

figure('Color', white);
[cAOV1a, mAOV1a, hAOV1a, nmsAOV1a] = ...
    multcompare(statsAOV1a, 'CType', 'hsd'); %,'Dimension',[1 2]);


% mixed-effects ANOVA (same as 1a) without Hoku & Lono

% y1c = CC_HR_DATA.ensWsstHr(CC_HR_DATA.altId ~= 1 & CC_HR_DATA.altId ~= 6) ;
y1c = CC_HR_DATA.iTotalKE(CC_HR_DATA.altId ~= 1 & CC_HR_DATA.altId ~= 6) ;
g1c = CC_HR_DATA.epoch(CC_HR_DATA.altId ~= 1 & CC_HR_DATA.altId ~= 6) ;
g2c = CC_HR_DATA.altId(CC_HR_DATA.altId ~= 1 & CC_HR_DATA.altId ~= 6) ;

[pAOV1c, tblAOV1c, statsAOV1c, termsAOV1c] = ...
    anovan( y1c, {g1c, g2c}, ...
    'model', 1, 'random', 2, ... %'interaction', ...
    'varnames', {'Epoch','AnimalID'} );

figure('Color', white);
[cAOV1c, mAOV1c, hAOV1c, nmsAOV1c] = ...
    multcompare(statsAOV1c, 'CType', 'hsd'); %,'Dimension',[1 2]);

% two-way ANOVA: ensWsstHr ~ fixed(epoch * id)

y1a = CC_HR_DATA.iTotalKE;
% y2a = CC_HR_DATA.ensWsstHr;
g1a = CC_HR_DATA.epoch;
g2a = CC_HR_DATA.altId;

[pAOV1e, tblAOV1e, statsAOV1e, termsAOV1e] = ...
    anovan( y1a, {g1a, g2a}, ...
    'model', 'interaction', ...
    'varnames', {'Epoch','AnimalID'} );

% [pAOV1f, tblAOV1f, statsAOV1f, termsAOV1f] = ...
%     anovan( y2a, {g1a, g2a}, ...
%     'model', 'interaction', ...
%     'varnames', {'Epoch','AnimalID'} );

figure('Color', white);
[cAOV1e, mAOV1e, hAOV1e, nmsAOV1e] = ...
    multcompare(statsAOV1e, 'CType', 'hsd','Dimension',[1 2]);

% figure('Color', white);
% [cAOV1f, mAOV1f, hAOV1f, nmsAOV1f] = ...
%     multcompare(statsAOV1f, 'CType', 'hsd','Dimension',[1 2]);


%%

mean(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 1))
std(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 1))
max(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 1))
min(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 1))

mean(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 2))
std(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 2))
max(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 2))
min(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 2))

mean(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 3))
std(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 3))
max(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 3))
min(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 3))

mean(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 4))
std(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 4))
max(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 4))
min(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 4))

mean(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 5))
std(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 5))
max(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 5))
min(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 5))

mean(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 6))
std(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 6))
max(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 6))
min(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 6))

%%

mean(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 2))
std(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 2))
max(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 2))
min(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 2))

mean(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 3))
std(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 3))
max(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 3))
min(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 3))

mean(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 4))
std(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 4))
max(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 4))
min(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 4))

mean(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 5))
std(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 5))
max(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 5))
min(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & CC_HR_DATA.altId == 5))


%%


mean(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 2))
std(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 2))
max(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 2))
min(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 2))
size(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 2))

mean(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 3))
std(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 3))
max(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 3))
min(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 3))
size(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 3))

mean(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 4))
std(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 4))
max(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 4))
min(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 4))
size(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 4))

mean(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 5))
std(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 5))
max(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 5))
min(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 5))
size(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & CC_HR_DATA.altId == 5))

%%

girth(2517,1) = zeros;
girth(1:185,1) = HOKU.girth;
girth(186:720,1) = LIHO.girth;
girth(721:1303,1) = KOLOHE.girth;
girth(1304:1668,1) = HUA.girth;
girth(1669:2200,1) = NOA.girth;
girth(2201:2517,1) = LONO.girth;


%%

figure('Color', white);
swarmchart(CC_HR_DATA.altId, log10(CC_HR_DATA.iLinKE), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', red);
hold on; 
swarmchart(CC_HR_DATA.altId, log10(CC_HR_DATA.iRotKE), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', green);
hold off;
xticklabels({'','63HF','01L5','6JK5','83H1','9ON6','9FL3',''});
ylabel('log10 KE, J·s^-^1');
legend('iLinKE','iRotKE');
grid;

%%

figure('Color', white);
swarmchart(CC_HR_DATA.altId(CC_HR_DATA.epoch ~= 3), ...
    log10(CC_HR_DATA.iLinKE(CC_HR_DATA.epoch ~= 3)), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', red);
hold on; 
swarmchart(CC_HR_DATA.altId(CC_HR_DATA.epoch == 3), ...
    log10(CC_HR_DATA.iLinKE(CC_HR_DATA.epoch == 3)), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', goldenrod);

swarmchart(CC_HR_DATA.altId(CC_HR_DATA.epoch ~= 3), ...
    log10(CC_HR_DATA.iRotKE(CC_HR_DATA.epoch ~= 3)), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', green);
swarmchart(CC_HR_DATA.altId(CC_HR_DATA.epoch == 3), ...
    log10(CC_HR_DATA.iRotKE(CC_HR_DATA.epoch == 3)), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', cyan);


hold off;
xticklabels({'','63HF','01L5','6JK5','83H1','9ON6','9FL3',''});
ylabel('log10 KE, J·s^-^1');
grid;
legend('iLinKE: E1 & E5','iLinKE: E3', 'iRotKE: E1 & E5', 'iRotKE: E3');

%%

figure('Color', white);

subplot(121); 

swarmchart(CC_HR_DATA.mass(CC_HR_DATA.epoch ~= 3), ...
    log10(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch ~= 3)), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', red);
hold on; 
swarmchart(CC_HR_DATA.mass(CC_HR_DATA.epoch == 3), ...
    log10(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 3)), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', cyan);
hold off;
% xticklabels({'','63HF','01L5','6JK5','83H1','9ON6','9FL3',''});
xlabel('Mass, kg');
ylabel('log_1_0 iKE, J·s^-^1');
grid;
legend('iKE: E1 & E5','iKE: E3');

subplot(122);

swarmchart(CC_HR_DATA.girth(CC_HR_DATA.epoch ~= 3), ...
    log10(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch ~= 3)), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', red);
hold on; 
swarmchart(CC_HR_DATA.girth(CC_HR_DATA.epoch == 3), ...
    log10(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 3)), ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', cyan);
hold off;
% xticklabels({'','63HF','01L5','6JK5','83H1','9ON6','9FL3',''});
xlabel('Girth, m');
ylabel('log_1_0 iKE, J·s^-^1');
grid;
legend('iKE: E1 & E5','iKE: E3');

%%

figure('Color', white);

% subplot(121);

swarmchart( log10(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 1)), ...
    log10(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1)), ...
    'd', 'MarkerEdgeColor', black, 'MarkerFaceColor', green);
hold on; 
swarmchart( log10(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 3)), ...
    log10(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 3)), ...
    'v', 'MarkerEdgeColor', black, 'MarkerFaceColor', cyan);
swarmchart( log10(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 5)), ...
    log10(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5)), ...
    '^', 'MarkerEdgeColor', black, 'MarkerFaceColor', red);

plot(log10(CC_HR_DATA.iHr), ( log10(CC_HR_DATA.iHr) .* 1.56) - 4.04, ':', ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
plot(log10(CC_HR_DATA.iHr), ( log10(CC_HR_DATA.iHr) .* 1.29) - 4.04, ':', ...
    'Color', [0.9 0.9 0.9], 'LineWidth', 1.1);
plot(log10(CC_HR_DATA.iHr), ( log10(CC_HR_DATA.iHr) .* 1.83) - 4.04, ':', ...
    'Color', [0.9 0.9 0.9], 'LineWidth', 1.1);

hold off;

xlabel('log_1_0(if_H), beats · min^-^1', 'FontSize', 15);
ylabel('log_1_0(iKE), J·s^-^1', 'FontSize', 15);

grid;
legend('iKE, E1','iKE, E3','iKE, E5', 'log10(iKE) ~ log10(if_H)');

%%

figure('Color', white);

% subplot(121);

swarmchart(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 1), ...
    log10(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1)), ...
    'd', 'MarkerEdgeColor', black, 'MarkerFaceColor', green);
hold on; 
swarmchart(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 3), ...
    log10(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 3)), ...
    'v', 'MarkerEdgeColor', black, 'MarkerFaceColor', cyan);
swarmchart(CC_HR_DATA.iHr(CC_HR_DATA.epoch == 5), ...
    log10(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5)), ...
    '^', 'MarkerEdgeColor', black, 'MarkerFaceColor', red);
% plot(CC_HR_DATA.iHr, (CC_HR_DATA.iHr .* 0.0105181)-1.90, '-', ...
%     'Color', [0.6 0.6 0.6], 'LineWidth', 1.1);\
text(50,-2,'(a)');


hold off;
% xticklabels({'','63HF','01L5','6JK5','83H1','9ON6','9FL3',''});
xlabel('ifH, beats · min^-^1');
ylabel('log_1_0 iKE, J·s^-^1');
% title('(a)', 'FontSize', 14, 'FontWeight', 'normal', ...
%     'FontName', 'Palatino');
% ax = gca;
% ax.TitleHorizontalAlignment = 'right';
grid;
legend('iKE, E1','iKE, E3','iKE, E5');


%%

%subplot(122);
figure('Color', white);

swarmchart(CC_HR_DATA.epoch(CC_HR_DATA.epoch == 1), ...
    log10(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1)), ...
    'd', 'MarkerEdgeColor', black, 'MarkerFaceColor', green);
hold on; 
swarmchart(CC_HR_DATA.epoch(CC_HR_DATA.epoch == 3), ...
    log10(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 3)), ...
    'v', 'MarkerEdgeColor', black, 'MarkerFaceColor', cyan);
swarmchart(CC_HR_DATA.epoch(CC_HR_DATA.epoch == 5), ...
    log10(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5)), ...
    '^', 'MarkerEdgeColor', black, 'MarkerFaceColor', red);
hold off;
xlabel('Epoch');
xticks(1:2:5);
xticklabels({'Baseline','Apnea','Recovery'});
ylabel('log_1_0 iKE, J·s^-^1');
% title('(b)', 'FontSize', 14, 'FontWeight', 'normal', ...
%     'FontName', 'Palatino');
% ax = gca;
% ax.TitleHorizontalAlignment = 'right';
grid;
legend('iKE, E1','iKE, E3','iKE, E5');


%%

figure('Color', white);

swarmchart(CC_HR_DATA.epoch(CC_HR_DATA.epoch == 1 & ...
    (CC_HR_DATA.altId >= 2 & CC_HR_DATA.altId <= 5) ), ...
    log10(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & ...
    (CC_HR_DATA.altId >= 2 & CC_HR_DATA.altId <= 5))), ...
    'd', 'MarkerEdgeColor', black, 'MarkerFaceColor', green);
hold on; 
swarmchart(CC_HR_DATA.epoch(CC_HR_DATA.epoch == 3 & ...
    (CC_HR_DATA.altId >= 2 & CC_HR_DATA.altId <= 5)), ...
    log10(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 3 & ...
    (CC_HR_DATA.altId >= 2 & CC_HR_DATA.altId <= 5))), ...
    'v', 'MarkerEdgeColor', black, 'MarkerFaceColor', cyan);
swarmchart(CC_HR_DATA.epoch(CC_HR_DATA.epoch == 5 & ...
    (CC_HR_DATA.altId >= 2 & CC_HR_DATA.altId <= 5)), ...
    log10(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & ...
    (CC_HR_DATA.altId >= 2 & CC_HR_DATA.altId <= 5))), ...
    '^', 'MarkerEdgeColor', black, 'MarkerFaceColor', red);
hold off;
xlabel('Epoch', 'FontName', 'Palatino', 'FontSize', 13);
xticks(1:2:5);
xticklabels({'Baseline','Apnea','Recovery'});
ylabel('log_1_0 iKE, J·s^-^1', 'FontName', 'Palatino', 'FontSize', 15);
% title('(b)', 'FontSize', 14, 'FontWeight', 'normal', ...
%     'FontName', 'Palatino');
% ax = gca;
% ax.TitleHorizontalAlignment = 'right';
grid;
legend('iKE, E1','iKE, E3','iKE, E5');


%%

figure('Color', white);

swarmchart(CC_HR_DATA.epoch(CC_HR_DATA.epoch == 1 & ...
    (CC_HR_DATA.altId >= 2 & CC_HR_DATA.altId <= 5) ), ...
    (CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & ...
    (CC_HR_DATA.altId >= 2 & CC_HR_DATA.altId <= 5))), ...
    'd', 'MarkerEdgeColor', black, 'MarkerFaceColor', green);
hold on; 
swarmchart(CC_HR_DATA.epoch(CC_HR_DATA.epoch == 3 & ...
    (CC_HR_DATA.altId >= 2 & CC_HR_DATA.altId <= 5)), ...
    (CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 3 & ...
    (CC_HR_DATA.altId >= 2 & CC_HR_DATA.altId <= 5))), ...
    'v', 'MarkerEdgeColor', black, 'MarkerFaceColor', cyan);
swarmchart(CC_HR_DATA.epoch(CC_HR_DATA.epoch == 5 & ...
    (CC_HR_DATA.altId >= 2 & CC_HR_DATA.altId <= 5)), ...
    (CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5 & ...
    (CC_HR_DATA.altId >= 2 & CC_HR_DATA.altId <= 5))), ...
    '^', 'MarkerEdgeColor', black, 'MarkerFaceColor', red);
hold off;
xlabel('Epoch', 'FontName', 'Palatino', 'FontSize', 13);
xticks(1:2:5);
xticklabels({'Baseline','Apnea','Recovery'});
ylabel('iKE, J·s^-^1', 'FontName', 'Palatino', 'FontSize', 15);
% title('(b)', 'FontSize', 14, 'FontWeight', 'normal', ...
%     'FontName', 'Palatino');
% ax = gca;
% ax.TitleHorizontalAlignment = 'right';
grid;
legend('iKE, E1','iKE, E3','iKE, E5');


%% Mann-Whitney U test for unpaired samples

% this compares iKEadj (which is iKE with Hoku and Lono removed, due to
% non-participation in Epoch 5 trial

statsE1E5 = mwwtest( (CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1 & ...
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

%% Mann-Whitney U test for unpaired samples

% this compares iKE w/o removal of Hoku and Lono
% non-participation in Epoch 5 trial

% E1:E5 comparison
statsE1E5 = mwwtest( (CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1))' , ...
(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5))' )

% E1:E3 comparison
statsE1E3 = mwwtest( (CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 1))' , ...
(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 3))' )

% E3:E5 comparison
statsE3E5 = mwwtest( (CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 3))' , ...
(CC_HR_DATA.iTotalKE(CC_HR_DATA.epoch == 5))' )

effectE1E5 = statsE1E5.Z / sqrt(sum(statsE1E5.n))
effectE3E5 = statsE3E5.Z / sqrt(sum(statsE3E5.n))
effectE1E3 = statsE1E3.Z / sqrt(sum(statsE1E3.n))

%% Mann-Whitney U test for unpaired samples

% this compares iRotKE w/o removal of Hoku and Lono
% non-participation in Epoch 5 trial

% E1:E5 comparison
statsE1E5 = mwwtest( (CC_HR_DATA.iRotKE(CC_HR_DATA.epoch == 1))' , ...
(CC_HR_DATA.iRotKE(CC_HR_DATA.epoch == 5))' )

% E1:E3 comparison
statsE1E3 = mwwtest( (CC_HR_DATA.iRotKE(CC_HR_DATA.epoch == 1))' , ...
(CC_HR_DATA.iRotKE(CC_HR_DATA.epoch == 3))' )

% E3:E5 comparison
statsE3E5 = mwwtest( (CC_HR_DATA.iRotKE(CC_HR_DATA.epoch == 3))' , ...
(CC_HR_DATA.iRotKE(CC_HR_DATA.epoch == 5))' )

effectE1E5 = statsE1E5.Z / sqrt(sum(statsE1E5.n))
effectE3E5 = statsE3E5.Z / sqrt(sum(statsE3E5.n))
effectE1E3 = statsE1E3.Z / sqrt(sum(statsE1E3.n))


%% Mann-Whitney U test for unpaired samples

% this compares iLinKE w/o removal of Hoku and Lono
% non-participation in Epoch 5 trial

% E1:E5 comparison
statsE1E5 = mwwtest( (CC_HR_DATA.iLinKE(CC_HR_DATA.epoch == 1))' , ...
(CC_HR_DATA.iLinKE(CC_HR_DATA.epoch == 5))' )

% E1:E3 comparison
statsE1E3 = mwwtest( (CC_HR_DATA.iLinKE(CC_HR_DATA.epoch == 1))' , ...
(CC_HR_DATA.iLinKE(CC_HR_DATA.epoch == 3))' )

% E3:E5 comparison
statsE3E5 = mwwtest( (CC_HR_DATA.iLinKE(CC_HR_DATA.epoch == 3))' , ...
(CC_HR_DATA.iLinKE(CC_HR_DATA.epoch == 5))' )

effectE1E5 = statsE1E5.Z / sqrt(sum(statsE1E5.n))
effectE3E5 = statsE3E5.Z / sqrt(sum(statsE3E5.n))
effectE1E3 = statsE1E3.Z / sqrt(sum(statsE1E3.n))


%% Mann-Whitney U test for unpaired samples

% this compares iLinKE w/o removal of Hoku and Lono
% non-participation in Epoch 5 trial

% E1:E5 comparison
statsE1LinRot = mwwtest( (CC_HR_DATA.iLinKE(CC_HR_DATA.epoch == 1))' , ...
(CC_HR_DATA.iRotKE(CC_HR_DATA.epoch == 1))' )

% E1:E3 comparison
statsE1E3 = mwwtest( (CC_HR_DATA.iLinKE(CC_HR_DATA.epoch == 1))' , ...
(CC_HR_DATA.iLinKE(CC_HR_DATA.epoch == 3))' )

% E3:E5 comparison
statsE3E5 = mwwtest( (CC_HR_DATA.iLinKE(CC_HR_DATA.epoch == 3))' , ...
(CC_HR_DATA.iLinKE(CC_HR_DATA.epoch == 5))' )

effectE1E5 = statsE1E5.Z / sqrt(sum(statsE1E5.n))
effectE3E5 = statsE3E5.Z / sqrt(sum(statsE3E5.n))
effectE1E3 = statsE1E3.Z / sqrt(sum(statsE1E3.n))

