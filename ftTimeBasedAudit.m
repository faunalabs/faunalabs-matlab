%%  ftTimeBasedAudit.m

%   written by Dave Haas on 2-3 July 2021
%
%   This is a data-resize tool that uses time (local or seconds) on the 
%   x-axis, allowing the user to annotate time series data. This is used
%   for aligning cues entered in ftDataEntrySheets.m, and for annotating
%   time series data for behavioral, physioloigcal, and/or acoustic cues.

clc;
%clear;

%%  Load global TAG_PATHS, specify tag data set, and load RAW & PRH files

global TAG_PATHS;

if (exist(TAG_PATHS.RAW,'dir'))
    fprintf('RAW data folder specified in TAG_PATHS is present.\n');
else
    fprintf('RAW data folder specified in TAG_PATHS is not present.\n');
    fprintf('Fix missing RAW data folder and retry.\n');
    return;
end

if (exist(TAG_PATHS.PRH,'dir'))
    fprintf('PRH data folder specified in TAG_PATHS is present.\n');
else
    fprintf('PRH data folder specified in TAG_PATHS is not present.\n');
    fprintf('Fix missing PRH data folder and retry.\n');
    return;
end

if (exist(TAG_PATHS.METADATA,'dir'))
    fprintf('Meta data folder specified in TAG_PATHS is present.\n');
else
    fprintf('Meta data folder specified in TAG_PATHS is not present.\n');
    fprintf('Fix missing meta data folder and retry.\n');
    return;
end


%% Get tag data file of interest and read in the RAW and PRH data

tagStr = input('Which tag file do you want to work with? ','s');

%% Allow the users to specify whether to use RAW or PRH (decimated) data

fprintf('Specify the data set you want to use for this analysis.\n');
fprintf('\tRAW for raw (undecimated) data\n');
fprintf('\tPRH for decimated (25 Hz) optical and kinematic data\n');
workingDataSetStr = lower(input('Your choice (raw/prh) [raw]: ','s'));
switch (workingDataSetStr)
    case 'prh'
        useRaw = false;
        fprintf('Using PRH data set.\n');
    case 'raw'
        useRaw = true;
        fprintf('Using RAW data set.\n');
    case ''
        useRaw = true;
        fprintf('Using PRH data set.\n');
    otherwise
        useRaw = true;
        fprintf('Using PRH data set.\n');
end

if (useRaw)
    
    % use RAW file and set time_s and time using optics
    rawTargetFile = sprintf('%s/%sraw.mat', TAG_PATHS.RAW, tagStr);
    if (exist(rawTargetFile,'file'))
        fprintf('RAW file %s is present. Loading here for analysis...\n', ...
            tagStr);
        load(rawTargetFile);
    else
        fprintf('RAW file %s is not present at %s. Check and retry.\n', ...
            tagStr);
        return;
    end
    
else
    
    prhTargetFile = sprintf('%s/%sprh.mat', TAG_PATHS.PRH, tagStr);
    if (exist(prhTargetFile,'file'))
        fprintf('PRH file %s is present. Loading here for analysis...\n', ...
            tagStr);
        load(prhTargetFile);
        
        if (exist('O','var'))
            fprintf('Building time scales using TIME struct\n');
            % Place TIME struct values into time scales commonly used in analysis
            time = TIME.time;
            time_s = TIME.time_s;
        else
            fprintf('There does not appear to be a TIME struct in PRH.\\n');
            fprintf('Review the PRH file and possible re-run ftPostProcessing.m\n');
            return;
        end
    else
        fprintf('PRH file %s is not present at %s. Check and retry.\n', ...
            tagStr);
        return;
    end
end

% load METADATA file too...
metaDataTargetFile = sprintf('%s/%smeta.mat', TAG_PATHS.METADATA, tagStr);
if (exist(metaDataTargetFile,'file'))
    fprintf('METADATA file %s is present. Loading here for analysis...\n', ...
        tagStr);
    load(metaDataTargetFile);
else
    fprintf('METADATA file %s is not present at %s. Check and retry.\n', ...
        tagStr);
    return;
end


%% Plot TRIAL cues over time series

maxCueIndex = numel(TRIAL.CUES.time);

originalFirstCueTime = TRIAL.CUES.time(1);
originalLastCueTime = TRIAL.CUES.time(end);

f2 = figure;
p1 = plot(OPTICS.time, OPTICS.led2, 'b');
grid;
hold on;
for i = 1:maxCueIndex
    switch(string(TRIAL.CUES.type(i)))
        case '1'
            plot(TRIAL.CUES.time(i), max(p1.YData), 'bd');            
        case '2'
            plot(TRIAL.CUES.time(i), max(p1.YData), 'g^');
        case '3'
            plot(TRIAL.CUES.time(i), max(p1.YData), 'rv');
        case '4'
            plot(TRIAL.CUES.time(i), max(p1.YData), 'm^');
        case '*'
            % plot(TRIAL.CUES.time(i), max(p1.YData), 'b*');
            xline(TRIAL.CUES.time(i), 'b--');
    end
end
hold off;

%% adjust TRIAL.CUES.time vs. opticsTime (raw) or time (prh)

durationDifference = time(1) - TRIAL.CUES.time(1);
secondsDifference = seconds(durationDifference);
fprintf('time(1) vs TRIAL.CUES(1) difference: %d seconds\n', ...
    secondsDifference);
coarseAdjustChoice = lower(input('Perform coarse adjustment (y/n) [y]): ','s'));
if (isempty(coarseAdjustChoice))
    coarseAdjustChoice = 'y';
end
switch(coarseAdjustChoice)
    case 'y'
        for i = 1:maxCueIndex
            TRIAL.CUES.time(i) = TRIAL.CUES.time(i) + durationDifference;
        end
        fprintf('Completed adjustments to TRIAL.CUE.time(*).\n');
    otherwise
        fprintf('Making no adjustments to TRIAL.CUE.time(*).\n');
end

if (isempty(f2.findobj))
    f2 = figure;
else
    figure(f2);    
end

p1 = plot(OPTICS.time, OPTICS.led2, 'b');
grid;
hold on;
for i = 1:maxCueIndex
    switch(string(TRIAL.CUES.type(i)))
        case '1'
            plot(TRIAL.CUES.time(i), max(p1.YData), 'bd');            
        case '2'
            plot(TRIAL.CUES.time(i), max(p1.YData), 'g^');
        case '3'
            plot(TRIAL.CUES.time(i), max(p1.YData), 'rv');
        case '4'
            plot(TRIAL.CUES.time(i), max(p1.YData), 'm^');
        case '*'
            %plot(TRIAL.CUES.time(i), max(p1.YData), 'b*');
            xline(TRIAL.CUES.time,'m--');
    end
end
hold off;

cueAlignComplete = false;
updatePlot = false;

while (~cueAlignComplete)
    
    [x,y,button] = ginput(1);
    
    switch(char(button))
        case '>'
            for i = 1:maxCueIndex
            	TRIAL.CUES.time(i) = TRIAL.CUES.time(i) + durationDifference;
            end
            cueAlignComplete = false
            updatePlot = true;
            
        case '<'
            for i = 1:maxCueIndex
            	TRIAL.CUES.time(i) = TRIAL.CUES.time(i) - durationDifference;
            end            
            cueAlignComplete = false
            updatePlot = true;

        case 'm'
            secondsUpdateText = 'Enter +/- seconds adjustment: ';
            secondsUpdate = str2double(input(secondsUpdateText,'s'));
            if (isnumeric(secondsUpdate))
                fprintf('Updating by %d seconds...\n', secondsUpdate);
                newDurationOffset = duration([ 0 0 secondsUpdate]);
                for i = 1:maxCueIndex
                    TRIAL.CUES.time(i) = TRIAL.CUES.time(i) + newDurationOffset;
                end            
                cueAlignComplete = false
                updatePlot = true;
            else
                fprintf('That was an invalid number of seconds. Try again.\n');
            end
            
        case 'g'
            
            alignDurationDiff = originalFirstCueTime - TRIAL.CUES.time(1);
            fprintf('Alignment difference: %d seconds.\n', ...
                seconds(alignDurationDiff) );
            alignSaveText = 'You said aligment looks good. Want to save? (y/n) [n]: ';
            alignSaveStr = lower(input(alignSaveText, 's'));
            switch(alignSaveStr)
                case 'y'
                    TRIAL.tagTimelogTimeDifference = alignDurationDiff;
                    TRIAL.CUES.tagTimeOffset = alignDurationDiff;
                    TRIAL.CUES.timeSyncComplete = 'yes';
                    TRIAL.CUES.timeSyncBy = 'Dave Haas';
                    save(metaDataTargetFile, 'TRIAL');
                otherwise
                    fprintf('Exiting without saving. Re-run to save.\n');
            end
            
            cueAlignComplete = true;
        case 'q'
            fprintf('Quitting cue alignment tool.\n');
            cueAlignComplete = true;
            updatePlot = false;
        otherwise
            cueAlignComplete = false;
            updatePlot = false;
    end
    
    if (updatePlot)

        if (isempty(f2.findobj))
            f2 = figure;
        else
            figure(f2);    
        end

        p1 = plot(OPTICS.time, OPTICS.led2, 'b');
        grid;
        hold on;
        for i = 1:maxCueIndex
            switch(string(TRIAL.CUES.type(i)))
                case '1'
                    plot(TRIAL.CUES.time(i), max(p1.YData), 'bd');            
                case '2'
                    plot(TRIAL.CUES.time(i), max(p1.YData), 'g^');
                case '3'
                    plot(TRIAL.CUES.time(i), max(p1.YData), 'rv');
                case '4'
                    plot(TRIAL.CUES.time(i), max(p1.YData), 'm^');
                case '*'
                    %plot(TRIAL.CUES.time(i), max(p1.YData), 'b*');
                    xline(TRIAL.CUES.time,'m--');
            end
        end
        hold off;
        grid;
        updatePlot = false;
        
    end
end

%% set time_s and time, depending on whether we're running RAW or PRH

if (useRaw)
    time = OPTICS.time;
    time_s = OPTICS.time_s;
else
    time = TIME.time;
    time_s = TIME.time_s;
end


%% audit trim selection test

% give users the choice of which optics to work with
%   1 = ledX (raw optical data)
%   2 = dcoLedX (DC removal optical data)
%   3 = lpLedX (low-pass, high frequency-removed optical data)
%   4 = hrLedX (raw - Savitzy-Golay filtered optical data)
%   5 = bpHrLedX (band-passed 0.33-3.5 Hz optical data)

clc;
fprintf('Choose one of these optical data sets to use for auditing:\n');
opticsDataChoiceText = '[1]ledX [2]dcoLedX [3]lpLedx [4]hrLedX [5]bpHrLedX : ';
opticsDataChoice = str2double(input(opticsDataChoiceText,'s'));
switch(opticsDataChoice)
    case 1
        figureName = 'ftAudit using raw optical data (no filtering)';
        fprintf('Loading raw optical data into LEDx variable.\n');
        LED1 = OPTICS.led1;
        LED2 = OPTICS.led2;
        LED3 = OPTICS.led3;
        LED4 = OPTICS.led4;
    case 2
        figureName = 'ftAudit using dco filtered optical data';
        fprintf('Loading DC-removed optical data into LEDx variable.\n');
        LED1 = dcoLed1;
        LED2 = dcoLed2;
        LED3 = dcoLed3;
        LED4 = dcoLed4;
    case 3
        figureName = 'ftAudit using dco+lp filtered optical data';
        fprintf('Loading DC-removed + low-passed optical data into LEDx variable.\n');
        LED1 = lpLed1;
        LED2 = lpLed2;
        LED3 = lpLed3;
        LED4 = lpLed4;
    case 4
        figureName = 'ftAudit using dco+lp-SG filtered optical data';
        fprintf('Loading -DC+LP-SG optical data into LEDx variable.\n');
        LED1 = hrLed1;
        LED2 = hrLed2;
        LED3 = hrLed3;
        LED4 = hrLed4;
    case 5
        figureName = 'ftAudit using dco+lp-SG+band-pass filtered optical data';
        fprintf('Loading -DC+LP-SG+BP optical data into LEDx variable.\n');
        LED1 = bpHrLed1;
        LED2 = bpHrLed2;
        LED3 = bpHrLed3;
        LED4 = bpHrLed4;
    otherwise
        opticsDataChoice = 1;
        figureName = 'ftAudit using RAW optical data';
        fprintf('Loading raw optical data into LEDx variable.\n');
        LED1 = OPTICS.led1;
        LED2 = OPTICS.led2;
        LED3 = OPTICS.led3;
        LED4 = OPTICS.led4;
end

% define a default windowSize (e.g.: 30 seconds) and windowDuration
windowSizeSeconds = 15;
windowDuration = duration([0 0 windowSizeSeconds]);

% create a time slice of the first windDuration seconds of the record
timeSlice  = [ time(1)  (time(1) + windowDuration) ];

% create a duration (h m s) called recordSize 
recordTime_s = time(end) - time(1);

cueStart = find(time == timeSlice(1));
cueEnd = find(time == timeSlice(2));

xLeft = time_s(1,cueStart);
xRight = time_s(1,cueEnd);
            
screenSize = get(0,'screensize');
markTime_s = zeros(1:2);

% do some work to plot existing cues
timeRange = [time(cueStart) time(cueEnd)];
currentCues = find(TRIAL.CUES.time > timeRange(1) & ...
    TRIAL.CUES.time <= timeRange(2));
if (~isempty(currentCues))
    currentCueTimes = TRIAL.CUES.time(currentCues);
    currentCueIds = TRIAL.CUES.type(currentCues);
    numCues = numel(currentCueTimes);
    cueX = zeros(1,numCues);
    for loopIndex = 1:numCues
        cueX(loopIndex) = find(time == currentCueTimes(loopIndex));
    end
    plotCues = true;
end

% Make the figure

f1 = figure;
f1.Name = figureName;
figSize = f1.Position;
figSize(3) = screenSize(3);
figSize(4) = screenSize(4) * 0.7;
f1.Position = figSize;  

subLed2 = subplot(3,3,[1:2]);
plot(time_s(cueStart:cueEnd), LED2(cueStart:cueEnd), 'b-');
if (plotCues)
    hold on;
    for i = 1:numel(cueX)
        yPos = subLed2.YLim(2) - ( (subLed2.YLim(2)-subLed2.YLim(1)) * 0.05 );
        text( time_s(cueX(i)), yPos, string(currentCueIds(i)), ...
            'FontSize',14,'FontWeight','bold','Color',[1 0 0] );
        xline( time_s(cueX(i)), 'r:' );
    end
    hold off;
end
xlabel('Time, seconds');
ylabel('Intensity (1050) nm');
xlim([xLeft xRight]);
grid;

subLed2Wsst = subplot(3,3,3);
[wsstLed2,f] = wsst(LED2(1,cueStart:cueEnd), opticsFs);
pcolor(time_s(1,cueStart:cueEnd), f, abs(wsstLed2) );
shading interp
xlabel('Seconds');
ylabel('Frequency (Hz)');
title('WSST - HR LED2 1050 nm');
ylim([0 3]);
grid;

subLed3 = subplot(3,3,[4:5]);
plot(time_s(cueStart:cueEnd), LED3(cueStart:cueEnd), 'r-');
xlabel('Time, seconds');
ylabel('Intensity (1200 nm)');
xlim([xLeft xRight]);
grid;

subLed3Wsst = subplot(3,3,6);
[wsstLed3,f] = wsst(LED3(1,cueStart:cueEnd), opticsFs);
pcolor(time_s(1,cueStart:cueEnd), f, abs(wsstLed3) );
shading interp
xlabel('Seconds');
ylabel('Frequency (Hz)');
title('WSST - HR LED3 1200 nm');
ylim([0 3]);
grid;

subLed1 = subplot(3,3,[7:8]);
plot(time_s(cueStart:cueEnd), LED1(cueStart:cueEnd), 'k-');
xlabel('Time, seconds');
ylabel('Intensity (ambient)');
xlim([xLeft xRight]);
grid;

subLed1Wsst = subplot(3,3,9);
[wsstLed1,f] = wsst(LED1(1,cueStart:cueEnd), opticsFs);
pcolor(time_s(1,cueStart:cueEnd), f, abs(wsstLed1) );
shading interp
xlabel('Seconds');
ylabel('Frequency (Hz)');
title('WSST - Ambient');
ylim([0 3]);
grid;



% if old leftX and rightX are sitting around in memory, nuke 'em...
if ( exist('leftX','var') && exist('rightX', 'var') )
    clear leftX rightX;
end

trimComplete = false;
updatePlot = false;

fprintf('[n]ext [p]revious [s]tart [e]nd [w]in(seconds) [r]ecord [q]uit\n');
while (~trimComplete)
    
    [x,y,button] = ginput(1);
    
	if (button == 1)
        button = char('1');     % left click
    elseif (button == 2)
        button = char('2');     % middle click
    elseif (button == 3)
        button = char('3');     % right click
    end
            
    switch(char(button))
        
        case '1'    % left click
            
            clickTimeSeconds = round(x,2);
            markTimeSeconds = time_s(1) + clickTimeSeconds;
            markCueStart = find( round(time_s,2) == markTimeSeconds);
            
            if (pingPong == 0)   % 1st audit markTime is not present
                markTime_s(1) = time_s(markCueStart(1));
                markTime(1) = time(markCueStart(1));
                fprintf('Audit start cue: %.1f seconds | %s\n', ...
                    round(markTime_s(1),1) , markTime(1) );
                pingPong = ~pingPong;
            else
                markTime_s(2) = time_s(markCueStart(1));
                markTime(2) = time(markCueStart(1));
                fprintf('Audit end cue: %.1f seconds | %s\n', ...
                    round(markTime_s(2),1) , markTime(2) );
                pingPong = ~pingPong;
            end
            
            figure(f1);
            hold on;
            xline(markTime_s,'m-','LineWidth',2);
            hold off;

            
        case '2'    % middle click
            
        case '3'    % right click

        case 'r'
            if (markTime_s(2) > markTime_s(1))  % cues are time+
                auditDuration = duration([0 0 (markTime_s(2) - ...
                    markTime_s(1))]);
            else                                % cues are time-
                auditDuration = duration([0 0 (markTime_s(1) - ...
                    markTime_s(2))]);
            end
            fprintf('Audit duration: %s\n', auditDuration);
            
        case 'n'
            timeSlice(1) = timeSlice(1) + windowDuration;
            timeSlice(2) = timeSlice(2) + windowDuration;
            cueStart = find(time == timeSlice(1));
            cueEnd = find(time == timeSlice(2));
            if (timeSlice(1) >= time(end) || timeSlice(2) >= time(end) )
                fprintf('Reached end of record, so showing last %d seconds.\n', ...
                    seconds(windowDuration) );
                timeSlice(1) = time(end) - windowDuration;
                timeSlice(2) = time(end);
                cueStart = find(time == timeSlice(1));
                cueEnd = find(time == timeSlice(2));
            else
                % do nothing, can probably just make this an if-end
            end
            fprintf('New cueStart: %d | cueEnd: %d\n', ...
                    cueStart, cueEnd);
            xLeft = time_s(1,cueStart);
            xRight = time_s(1,cueEnd);
            updatePlot = true;

        case 'p'
            timeSlice(1) = timeSlice(1) - windowDuration;
            timeSlice(2) = timeSlice(2) - windowDuration;
            cueStart = find(time == timeSlice(1));
            cueEnd = find(time == timeSlice(2));            
            if (timeSlice(1) <= time(1) || timeSlice(2) <= time(1) )
                fprintf('Reached start of record. Showing first %d seconds.\n', ...
                    seconds(windowDuration) );
                timeSlice(1) = time(1);
                timeSlice(2) = time(1) + windowDuration;
                cueStart = find(time == timeSlice(1));
                cueEnd = find(time == timeSlice(2));
            else
                % do nothing, can probably make this an if-end
            end
            fprintf('New cueStart: %d | cueEnd: %d\n', ...
                    cueStart, cueEnd);
            xLeft = time_s(1,cueStart);
            xRight = time_s(1,cueEnd);
            updatePlot = true;
        
        case 'j'
            jumpText = 'Jump to time in seconds: ';
            jumpCue = str2double(input(jumpText,'s'));
            if (jumpCue >= time_s(1) && jumpCue <= time_s(end))
                cueStart = find(time_s == jumpCue);
                timeSlice(1) = time(cueStart);
                timeSlice(2) = timeSlice(1) + windowDuration;
                cueEnd = find(time == timeSlice(2));
                xLeft = time_s(1,cueStart);
                xRight = time_s(1,cueEnd);
                updatePlot = true;
            else
                fprintf('Jump time is outside the scope of this trial.\n');
            end
            
        case 'w'
            windowSizeStr = input('Enter new window size in seconds: ','s');
            newWindowSize = str2double(windowSizeStr);
            if (newWindowSize >= 5)
                windowSize = newWindowSize;
                windowSizeSeconds = newWindowSize;
                windowDuration = duration([0 0 windowSizeSeconds]);
                timeSlice  = [ timeSlice(1) (timeSlice(1) + windowDuration) ];
                cueStart = find(time == timeSlice(1));
                cueEnd = find(time == timeSlice(2));
                xLeft = time_s(1,cueStart);
                xRight = time_s(1,cueEnd);
                updatePlot = true;
            else
                fprintf('Make sure time is 5 or more seconds.\n');
                updatePlot = false;
            end

        case 's'
            timeSlice(1) = time(1);
            timeSlice(2) = time(1) + windowDuration;
            cueStart = find(time == timeSlice(1));
            cueEnd = find(time == timeSlice(2));
            xLeft = time_s(1,cueStart);
            xRight = time_s(1,cueEnd);
            updatePlot = true;            
            
        case 'e'
            timeSlice(1) = time(end) - windowDuration;
            timeSlice(2) = time(end);
            cueStart = find(time == timeSlice(1));
            cueEnd = find(time == timeSlice(2));
            xLeft = time_s(1,cueStart);
            xRight = time_s(1,cueEnd);
            updatePlot = true;
            
        case 'q'
            fprintf('*** Quitting ftTimeAudit ***\n');
            close(f1);
            return;   
            
        otherwise
            fprintf('use [n] [p] [w] [e] [l] [r] or [q]\n');
    end

    if ( exist('leftX','var') && exist('rightX', 'var') )
        fprintf('Left and right trim coordinates captured. Continuing...\n');
        trimComplete = true;
    end
    
    if (updatePlot)
        
        % do some work to plot existing cues
        timeRange = [time(cueStart) time(cueEnd)];
        currentCues = find(TRIAL.CUES.time > timeRange(1) & ...
            TRIAL.CUES.time <= timeRange(2));
        if (~isempty(currentCues))
            currentCueTimes = TRIAL.CUES.time(currentCues);
            currentCueIds = TRIAL.CUES.type(currentCues);
            numCues = numel(currentCueTimes);
            fprintf('\tFound %d TRIAL.CUES to add...\n', numCues);
            cueX = zeros(1,numCues);
            for loopIndex = 1:numCues
                cueX(loopIndex) = find(time == currentCueTimes(loopIndex));
            end
            plotCues = true;
        end        
        
        figure(f1);
       
        subLed2 = subplot(3,3,[1:2]);
        plot(time_s(cueStart:cueEnd), LED2(cueStart:cueEnd), 'b-'); 
        if (plotCues)
            hold on;
            for i = 1:numel(cueX)
                yPos = subLed2.YLim(2) - ( (subLed2.YLim(2)-subLed2.YLim(1)) * 0.05 );
                text( time_s(cueX(i)), yPos, string(currentCueIds(i)), ...
                    'FontSize',14,'FontWeight','bold','Color',[1 0 0] );
                xline( time_s(cueX(i)), 'r:' );
            end
            hold off;
        end
        xlabel('Time, seconds');
        ylabel('Intensity (1050) nm');
        xlim([xLeft xRight]);
        grid;

        subLed2Wsst = subplot(3,3,3);
        [wsstLed2,f] = wsst(LED2(1,cueStart:cueEnd), opticsFs);
        pcolor(time_s(1,cueStart:cueEnd), f, abs(wsstLed2) );
        shading interp
        xlabel('Seconds');
        ylabel('Frequency (Hz)');
        title('WSST - HR LED2 1050 nm');
        ylim([0 3]);
        grid;

        subLed3 = subplot(3,3,[4:5]);
        plot(time_s(cueStart:cueEnd), LED3(cueStart:cueEnd), 'r-');
        xlabel('Time, seconds');
        ylabel('Intensity (1200 nm)');
        xlim([xLeft xRight]);
        grid;

        subLed3Wsst = subplot(3,3,6);
        [wsstLed3,f] = wsst(LED3(1,cueStart:cueEnd), opticsFs);
        pcolor(time_s(1,cueStart:cueEnd), f, abs(wsstLed3) );
        shading interp
        xlabel('Seconds');
        ylabel('Frequency (Hz)');
        title('WSST - HR LED3 1200 nm');
        ylim([0 3]);
        grid;

        subLed1 = subplot(3,3,[7:8]);
        plot(time_s(cueStart:cueEnd), LED1(cueStart:cueEnd), 'k-');
        xlabel('Time, seconds');
        ylabel('Intensity (ambient)');
        xlim([xLeft xRight]);
        grid;

        subLed1Wsst = subplot(3,3,9);
        [wsstLed1,f] = wsst(LED1(1,cueStart:cueEnd), opticsFs);
        pcolor(time_s(1,cueStart:cueEnd), f, abs(wsstLed1) );
        shading interp
        xlabel('Seconds');
        ylabel('Frequency (Hz)');
        title('WSST - Ambient');
        ylim([0 3]);
        grid;



        updatePlot = false;
        
    end
end


hold off;

inFmt = 'ss.SSS';
outFmt = 'hh:mm:ss.SSS';
leftD = duration([0 0 leftX], 'InputFormat', inFmt, 'Format', outFmt);
rightD = duration([0 0 rightX],'InputFormat',inFmt,'Format',outFmt);


