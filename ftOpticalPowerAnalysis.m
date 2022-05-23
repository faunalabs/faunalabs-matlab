%% ftOpticalPowerAnalysis.m
%
%   Written by Dave Haas from 1 September through 17 October 2021
%
%

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

% store screenSize now for later use in auto-sized plots

screenSize = get(0,'screensize');

%%

validChoice = false;

while( ~validChoice)
    
    fprintf('Which trial for optical analysis?\n');
    fprintf('\t1 = tt21_128f\n');
    fprintf('\t2 = tt21_134a\n');
    fprintf('\t3 = tt21_134b\n');
    fprintf('\t4 = tt21_141e\n');
    fprintf('\t5 = tt21_142c\n');
    fprintf('\t6 = tt21_142d\n');
    whichTrial = input('Choose [1-6]: ');
    
    
    switch(whichTrial)
        case 1
            fullPathFileName = '/Users/dave/Desktop/DTAG/raw/tt21_128fraw.mat';
            load(fullPathFileName);     
            validChoice = true;
            id = {'63HF'};
            altId = 1;
            mass   = 184.4091;
            length_m = 2.54;           % 254 cm
            girth  = 1.3081;         % 130.81 cm            
            trialId = 'tt21_128f';
        case 2
            fullPathFileName = '/Users/dave/Desktop/DTAG/raw/tt21_134araw.mat';
            load(fullPathFileName);
            validChoice = true;
            id = {'01L5'};
            altId = 2;
            mass = 164.3636;
            length_m = 2.4003;
            girth  = 1.2446;             
            trialId = 'tt21_134a';
        case 3
            fullPathFileName = '/Users/dave/Desktop/DTAG/raw/tt21_134braw.mat';
            load(fullPathFileName);
            validChoice = true;
            id = {'6JK5'};
            altId = 3;
            mass = 210.9091;
            length_m = 2.5908;
            girth  = 1.4478;                
            trialId = 'tt21_134b';
        case 4
            fullPathFileName = '/Users/dave/Desktop/DTAG/raw/tt21_141eraw.mat';
            load(fullPathFileName);
            validChoice = true;
            id = {'83H1'};
            altId = 4;
            mass = 152.2727;
            length_m = 2.413;
            girth  = 1.2573;              
            trialId = 'tt21_141e';
        case 5
            fullPathFileName = '/Users/dave/Desktop/DTAG/raw/tt21_142craw.mat';
            load(fullPathFileName);
            validChoice = true;
            id = {'9ON6'};
            altId = 5;
            mass = 203.5455;
            length_m = 2.4765;
            girth  = 1.3716;            
            trialId = 'tt21_142c';
        case 6
            fullPathFileName  = '/Users/dave/Desktop/DTAG/raw/tt21_142draw.mat'
            load(fullPathFileName);
            validChoice = true;
            id = {'9FL3'};
            altId = 6;
            mass   = 249.5455;
            length_m = 2.7432;
            girth  = 1.4542;            
            trialId = 'tt21_142d';
        otherwise
            fprintf('Invalid choice. Try again!\n');
            validChoice = false;
    end
end

%%

validChoice = false;

while (~validChoice)

    fprintf('Which epoch for optical analysis?\n');
    fprintf('\t1 = pre-apnea\n');
    fprintf('\t3 = apnea\n');
    fprintf('\t5 = recovery\n');
    whichEpoch = input('Choose [1,3,5]: ');    
    
    switch(whichEpoch)
        case 1
            if (exist('E1_CC','var'))
                fprintf('Epoch 1 selected.\n');
                validChoice = true;
                thisCC = E1_CC;
                thisCCE = E1_CC_ENERGY;
                thisDecData = E1_DECIMATED;                
            else
                fprintf('E1_CC was not found. Figure this out.\n');
                validChoice = false;
            end

        case 3
            if (exist('E3_CC','var'))
                fprintf('Epoch 3 selected.\n');
                validChoice = true;
                thisCC = E3_CC;
                thisCCE = E3_CC_ENERGY;
                thisDecData = E3_DECIMATED;                  
            else
                fprintf('E3_CC was not found. Figure this out.\n');
                validChoice = false;
            end
            
        case 5
            if (exist('E5_CC','var'))
                fprintf('Epoch 5 selected.\n');
                validChoice = true;
                thisCC = E5_CC;
                thisCCE = E5_CC_ENERGY;
                thisDecData = E5_DECIMATED;                   
            else
                fprintf('E5_CC was not found. Figure this out.\n');
                validChoice = false;
            end                  
        otherwise
            fprintf('That was not a valid choice. Try again.\n');
            validChoice = false;
            
    end
end

% fix thisCCE, just in case (this wasn't done previously a lot of times)

thisCCE.Time.Second = round(thisCCE.Time.Second, 2);

%%

validCouplets = 0;

newOCC = struct;

numCCRows = height(thisCC);

for index = 1:numCCRows-1
    
    fprintf('Processing row %d...\n', index);
    
    if (string(thisCC.cueType(index)) == 'ao' || ...
            string(thisCC.cueType(index)) == 'x' )
       
        if ( string(thisCC.cueType(index+1)) == '1' || ...
                string(thisCC.cueType(index+1)) == '2' || ...
                string(thisCC.cueType(index+1)) == '3' || ...
                string(thisCC.cueType(index+1)) == 'z' )
            
            fprintf('Valid couplet at row %d!\n', index);
            fprintf('\tccCue = %s\n', string(thisCC.cueType(index)));
            fprintf('\toptCue = %s\n', string(thisCC.cueType(index+1)));
            
            validCouplets = validCouplets + 1;
        
            newOCC.Time(validCouplets) = thisCC.Time(index);
            newOCC.ccCue(validCouplets) = thisCC.cueType(index);
            newOCC.opticsCue(validCouplets) = thisCC.cueType(index+1);
            
        else
            
            fprintf('No valid couplet at row %d...\n', index);
            
        end
       
    else
        
        fprintf('Something is triggering this at row %d\n', index);
        fprintf('\tthisCC.cueType(i) = %s\n', string(thisCC.cueType(index)));
        fprintf('\tthisCC.cueType(i+1) = %s\n', string(thisCC.cueType(index)));
        
    end
    
end

% assemble into a new time table

optTT = timetable(newOCC.Time(:), newOCC.ccCue(:), newOCC.opticsCue(:), ...
    'VariableNames', {'ccCue','opticsCue'} );



numOptTT = height(optTT);

optMeasures = struct;

optIndex = 0;

for optIndex = 1:numOptTT-1
   
    thisOptRange = timerange(optTT.Time(optIndex), optTT.Time(optIndex+1) );
    
    optMeasures.Time(optIndex) = optTT.Time(optIndex);
    optMeasures.TimeEnd(optIndex) = optTT.Time(optIndex+1);
    
    optMeasures.rmsLed1(optIndex) = rms(thisDecData.led1(thisOptRange));
    optMeasures.rmsLed2(optIndex) = rms(thisDecData.led2(thisOptRange));
    optMeasures.rmsLed3(optIndex) = rms(thisDecData.led3(thisOptRange));
    optMeasures.rmsLed4(optIndex) = rms(thisDecData.led4(thisOptRange));
    
    optMeasures.rmsFiltLed1(optIndex) = rms(thisDecData.filtLed1(thisOptRange));
    optMeasures.rmsFiltLed2(optIndex) = rms(thisDecData.filtLed2(thisOptRange));
    optMeasures.rmsFiltLed3(optIndex) = rms(thisDecData.filtLed3(thisOptRange));
    optMeasures.rmsFiltLed4(optIndex) = rms(thisDecData.filtLed4(thisOptRange));  
    
    optMeasures.powerFiltLed1(optIndex) = ...
        rms(thisDecData.filtLed1(thisOptRange).^2);
    optMeasures.powerFiltLed2(optIndex) = ...
        rms(thisDecData.filtLed2(thisOptRange).^2);    
    optMeasures.powerFiltLed3(optIndex) = ...
        rms(thisDecData.filtLed3(thisOptRange).^2);    
    optMeasures.powerFiltLed4(optIndex) = ...
        rms(thisDecData.filtLed4(thisOptRange).^2);
    
    optMeasures.dbFiltLed1 = pow2db(optMeasures.powerFiltLed1);
    optMeasures.dbFiltLed2 = pow2db(optMeasures.powerFiltLed2);
    optMeasures.dbFiltLed3 = pow2db(optMeasures.powerFiltLed3);
    optMeasures.dbFiltLed4 = pow2db(optMeasures.powerFiltLed4);
    
end

% Plot of all four optical channels in dB

figure('Color', white);
db1 = subplot(221);
plot(optMeasures.Time, optMeasures.dbFiltLed2, 'Color', blue);
xlabel('Time, local');
ylabel('dB');
title('LED Channel 2: \lambda = 1050 nm');
grid;
db2 = subplot(222);
plot(optMeasures.Time, optMeasures.dbFiltLed3, 'Color', red);
xlabel('Time, local');
ylabel('dB');
title('LED Channel 3: \lambda = 1200 nm');
grid;
db3 = subplot(223);
plot(optMeasures.Time, optMeasures.dbFiltLed1, 'Color', black);
xlabel('Time, local');
ylabel('dB');
title('Ambient 1');
grid;
db4 = subplot(224);
plot(optMeasures.Time, optMeasures.dbFiltLed4, 'Color', goldenrod);
xlabel('Time, local');
ylabel('dB');
title('Ambient 2');
grid;

figure('Color', white);
yyaxis left;
plot(optMeasures.Time, optMeasures.dbFiltLed2, 'Color', blue);
xlabel('Time, local');
ylabel('dB');
ylim([30 100]);
yyaxis right
plot(optMeasures.Time, optMeasures.dbFiltLed3, 'Color', red);
xlabel('Time, local');
ylabel('dB');
title('LED Channels 2 and 3: \lambda = 1050 nm & 1200 nm');
grid;
hold on;
plot(optMeasures.Time, optMeasures.dbFiltLed1, 'Color', black);
hold off;
ylim([30 100]);
legend('1050 nm','1200 nm','ambient');

%% Plot filtLed2 all + mean optical waveform

optIndex = 0;
legitRow = 0;

filtLed2Container(1:150,1,numOptTT) = zeros;

led2Plot = figure('Color', white, 'Position', [0 0 600 600] );

for optIndex = 1:numOptTT-1
   
    firstCue = string(optTT.ccCue(optIndex));
    nextCue  = string(optTT.opticsCue(optIndex+1));
    
    if ( contains(firstCue,'ao') && ( contains(nextCue,'1') || ...
            contains(nextCue,'3') ) )
       
        legitRow = legitRow + 1;
        
        thisOptRange = timerange(optTT.Time(optIndex), optTT.Time(optIndex+1) );
                       
        figure(led2Plot);
        hold on;
        plot( thisDecData.filtLed2(thisOptRange), 'Color', grey );
        hold off;
        
        lenThisRange = length(thisDecData.filtLed2(thisOptRange));
    
        filtLed2Container(1:lenThisRange,1,legitRow) = ...
            thisDecData.filtLed2(thisOptRange);
       
    end
    
end

meanFiltLed2 = mean(filtLed2Container(:,1,1:legitRow),3);

numFiltLed2 = legitRow;

legitCcRatio = legitRow / numOptTT;

figure(led2Plot);
hold on;
plot(meanFiltLed2, 'Color', blue, 'LineWidth', 1.5);
hold off;
xlabel('Samples, 50 Hz');
xlim([0 50]);
ylabel('A.U.');
titleTxt = sprintf('Optical channel 2: 1050 nm · Animal %s · Trial %s · Epoch %d', ...
    string(id), trialId, whichEpoch);
title(titleTxt, 'Interpreter','none');
grid;

legendTxt = sprintf('FiltLed2\nn=%d\nratio: %0.2f', legitRow, legitCcRatio);
legend(legendTxt);

fprintf('Mean filtLed2 power (db): %3.1f\n', mean(optMeasures.dbFiltLed2) );


%% Plot filtLed3 all + mean optical waveform%%

optIndex = 0;
legitRow = 0;

filtLed3Container(1:150,1,numOptTT) = zeros;

led3Plot = figure('Color', white, 'Position', [0 0 600 600] );

for optIndex = 1:numOptTT-1
   
    firstCue = string(optTT.ccCue(optIndex));
    nextCue  = string(optTT.opticsCue(optIndex+1));
    
    if ( contains(firstCue,'ao') && ( contains(nextCue,'2') || ...
            contains(nextCue,'3') ) )
       
        legitRow = legitRow + 1;
        
        thisOptRange = timerange(optTT.Time(optIndex), optTT.Time(optIndex+1) );
                       
        figure(led3Plot);
        hold on;
        plot( thisDecData.filtLed3(thisOptRange), 'Color', grey );
        hold off;
        
        lenThisRange = length(thisDecData.filtLed3(thisOptRange));
    
        filtLed3Container(1:lenThisRange,1,legitRow) = ...
            thisDecData.filtLed3(thisOptRange);
       
    end
    
end

meanFiltLed3 = mean(filtLed3Container(:,1,1:legitRow),3);

numFiltLed3 = legitRow;

figure(led3Plot);
hold on;
plot(meanFiltLed3, 'Color', red, 'LineWidth', 1.5);
hold off;
xlabel('Samples, 50 Hz');
xlim([0 50]);
ylabel('A.U.');
titleTxt = sprintf('Optical channel 3: 1200 nm · Animal %s · Trial %s · Epoch %d', ...
    string(id), trialId, whichEpoch);
title(titleTxt, 'Interpreter','none');
grid;


legitCcRatio = legitRow / numOptTT;

legendTxt = sprintf('FiltLed3\nn=%d\nratio: %0.2f', legitRow, legitCcRatio);
legend(legendTxt);

fprintf('Mean filtLed3 power (db): %3.1f\n', mean(optMeasures.dbFiltLed3) );

%%

filtLedsContainer(1:150,4,numOptTT-1) = zeros;

xScale = 0:0.02:2.98;

figure('Color', white);
s1 = subplot(221);
s2 = subplot(222);
s3 = subplot(223);
s4 = subplot(224);
hold on;


for optIndex = 1:numOptTT-1
    
    thisOptRange = timerange(optTT.Time(optIndex), optTT.Time(optIndex+1) );
        
    subplot(s1);
    hold on;
    plot( thisDecData.filtLed2(thisOptRange), 'Color', grey );
    hold off;
    
    subplot(s2);
    hold on;
	plot( thisDecData.filtLed3(thisOptRange), 'Color', grey );
    hold off;
    
    subplot(s3);
    hold on;
    plot( thisDecData.filtLed1(thisOptRange), 'Color', grey );
    hold off;
    
    subplot(s4);
    hold on;
    plot( thisDecData.filtLed4(thisOptRange), 'Color', grey );
    hold off;    
    
    
    lenThisRange = length(thisDecData.filtLed2(thisOptRange));
    
    filtLedsContainer(1:lenThisRange,2,optIndex) = ...
        thisDecData.filtLed2(thisOptRange);
    filtLedsContainer(1:lenThisRange,3,optIndex) = ...
        thisDecData.filtLed3(thisOptRange);
    filtLedsContainer(1:lenThisRange,1,optIndex) = ...
        thisDecData.filtLed1(thisOptRange);
    filtLedsContainer(1:lenThisRange,4,optIndex) = ...
        thisDecData.filtLed4(thisOptRange);
    
end


meanFiltLed2 = mean(filtLedsContainer(:,2,:),3);
meanFiltLed3 = mean(filtLedsContainer(:,3,:),3);
meanFiltLed1 = mean(filtLedsContainer(:,1,:),3);
meanFiltLed4 = mean(filtLedsContainer(:,4,:),3);

subplot(s1);
hold on;
plot( mean(filtLedsContainer(:,2,:),3), 'Color', blue, 'LineWidth', 2);
xlabel('Samples, 50 Hz');
ylabel('A.U.');
title('\lambda = 1050 nm');
xlim([0 50]);
grid;
hold off;

subplot(s2);
hold on;
plot( mean(filtLedsContainer(:,3,:),3), 'Color', red, 'LineWidth', 2);
xlabel('Samples, 50 Hz');
ylabel('A.U.');
title('\lambda = 1200 nm');
xlim([0 50]);
grid;
hold off;

subplot(s3);
hold on;
plot( mean(filtLedsContainer(:,1,:),3), 'Color', black, 'LineWidth', 2);
xlabel('Samples, 50 Hz');
ylabel('A.U.');
title('Ambient channel 1');
xlim([0 50]);
grid;
hold off;

subplot(s4);
hold on;
plot( mean(filtLedsContainer(:,4,:),3), 'Color', goldenrod, 'LineWidth', 2);
xlabel('Samples, 50 Hz');
ylabel('A.U.');
title('Ambient channel 2');
xlim([0 50]);
grid;
hold off;

%%


figure('Color', white); 
yyaxis left;
plot(optMeasures.Time, optMeasures.dbFiltLed2);
ylabel('dB');
yyaxis right;
plot(thisCCE.Time, thisCCE.iHR);
ylabel('beats · min^-^1');
title('Optical channel 2: \lambda = 1050 nm');
grid;
legend('1050 nm', 'ifH');


figure('Color', white); 
yyaxis left;
plot(optMeasures.Time, optMeasures.dbFiltLed3);
ylabel('dB');
yyaxis right;
plot(thisCCE.Time, thisCCE.iHR);
ylabel('beats · min^-^1');
title('Optical channel 3: \lambda = 1200 nm');
grid;
legend('1200 nm', 'ifH');


%%

rmsPowDbTimeTable = timetable( optMeasures.Time', optMeasures.TimeEnd', ...
    optMeasures.rmsLed1', optMeasures.rmsLed2', optMeasures.rmsLed3', ...
    optMeasures.rmsLed4', optMeasures.rmsFiltLed1', ...
    optMeasures.rmsFiltLed2', optMeasures.rmsFiltLed3', ...
    optMeasures.rmsFiltLed4', optMeasures.powerFiltLed1', ...
    optMeasures.powerFiltLed2', optMeasures.powerFiltLed3', ...
    optMeasures.powerFiltLed4', optMeasures.dbFiltLed1', ...
    optMeasures.dbFiltLed2', optMeasures.dbFiltLed3', ...
    optMeasures.dbFiltLed4', 'VariableNames', ...
    {'time_end','rmsLed1','rmsLed2','rmsLed3','rmsLed4',...
    'rmsFiltLed1','rmsFiltLed2','rmsFiltLed3','rmsFiltLed4', ...
    'powerFiltLed1','powerFiltLed2','powerFiltLed3','powerFiltLed4', ...
    'dbFiltLed1','dbFiltLed2','dbFiltLed3','dbFiltLed4' } );


%%

optPowerTT = synchronize(optTT, rmsPowDbTimeTable);

optPowerCCE = synchronize(optPowerTT, thisCCE);

tempTT = synchronize(optPowerCCE, thisDecData, 'first', 'nearest');

noNanRows = ( ~isnan(tempTT.iHR));

nextTempTT = tempTT(noNanRows,:);

noMissingRows = ( ~ismissing(nextTempTT.ccCue) | ...
    ~ismissing(nextTempTT.opticsCue));

finalTempTT = nextTempTT(noMissingRows,:);

finalTempTT.wsstHr = finalTempTT.wsstHr * 60;

finalTempTT = removevars(finalTempTT, {'iLinKE_xyz', 'iRotKE_xyz', ...
    'iTotalKE_xyz', 'mass', 'cleanLed1', 'cleanLed2', ...
    'cleanLed3', 'cleanLed4' } );

heightTT = height(finalTempTT);

altIdVar(heightTT,1)  = zeros;
massVar(heightTT,1)   = zeros;
lengthVar(heightTT,1) = zeros;
girthVar(heightTT,1)  = zeros;

altIdVar(:)  = altId;
massVar(:)   = mass;
lengthVar(:) = length_m;
girthVar(:)  = girth;

finalTempTT = addvars(finalTempTT, altIdVar, massVar, lengthVar, girthVar, ...
    'After', 'ID', 'NewVariableNames', ...
    {'altId', 'mass', 'length', 'girth'} );

CC_OPTICS_DATA = finalTempTT;

%% 

disp(CC_OPTICS_DATA)


fprintf('Review the above to make sure there are no NaTs, NaNs or missings...\n');

lookGood = lower(input('Everything look good (y/n) [n]: ','s'));

switch(lookGood)
    case 'y'
        fprintf('Okay... time to save all this... \n');
        okToProceed = true;
    otherwise
        fprintf('Looks like you saw something wrong. Fix and repeat.\n');
        okToProceed = false;
end


%%

if (okToProceed)

    switch(whichEpoch)

        case 1
            okToAppend = lower(input('Okay to append E1_CC_OPTICS (y/n) [n] ','s'));
            switch(okToAppend)
                case 'y'
                    E1_CC_OPTICS = CC_OPTICS_DATA;
                    fprintf('\nAppending E1_CC_OPTICS timetable to RAW file...\n');
                    save(fullPathFileName, 'E1_CC_OPTICS', '-append');  
                otherwise
                    fprintf('Skipping append. Redo this when ready...\n');
            end

        case 3
            okToAppend = lower(input('Okay to append E3_CC_OPTICS (y/n) [n] ','s'));
            switch(okToAppend)
                case 'y'        
                    E3_CC_OPTICS = CC_OPTICS_DATA;
                    fprintf('\nAppending E3_CC_OPTICS timetable to RAW file...\n');
                    save(fullPathFileName, 'E3_CC_OPTICS', '-append');
                otherwise
                    fprintf('Skipping append. Redo this when ready...\n');
            end                

        case 5
            okToAppend = lower(input('Okay to append E5_CC_OPTICS (y/n) [n] ','s'));
            switch(okToAppend)
                case 'y'           
                    E5_CC_OPTICS = CC_OPTICS_DATA;
                    fprintf('\nAppending E5_CC_OPTICS timetable to RAW file...\n');
                    save(fullPathFileName, 'E5_CC_OPTICS', '-append');  
                otherwise
                    fprintf('Skipping append. Redo this when ready...\n');
            end

        otherwise
            fprintf('Something horrible has happened. Bailing!\n');
            return;

    end

end

