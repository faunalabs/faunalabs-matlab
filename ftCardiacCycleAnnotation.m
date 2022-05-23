%%  ftCardiacCycleAnnotation.m

%   Written by Dave Haas between 9 July and 18 September 2021

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

%% Show entire trial with cue audits
% 
% screenSize = get(0,'screensize');
% 
% figAnalysisRegion = figure;
% figSize = figAnalysisRegion.Position;
% figSize(3) = screenSize(3);
% figSize(4) = screenSize(4) * 0.7;
% figAnalysisRegion.Position = figSize;  
% 
% p1a = subplot(211);
% plot(KINEMATICS.Time, KINEMATICS.odba, 'Color', green);
% hold on;
% yPos = p1a.YLim(2) - ((p1a.YLim(2) - p1a.YLim(1)) / 10);
% ftPlotCues(CUE.Time, CUE.type, yPos);
% hold off;
% xlabel('Time, local');
% ylabel('m/s^2');
% title('ODBA');
% grid;
% 
% p1b = subplot(212);
% plot(PRESSURE.Time, PRESSURE.depth, 'Color', blue');
% hold on;
% yPos = p1b.YLim(2) - ((p1b.YLim(2) - p1b.YLim(1)) / 10);
% ftPlotCues(CUE.Time, CUE.type, yPos);
% hold off;
% xlabel('Time, local');
% ylabel('Depth, meters');
% title('Tag Depth');
% set(gca,'YDir','reverse');
% grid;
% 
% linkaxes([p1a p1b],'x');

%% Select which of the five ranges to use for analysis

validEpoch = false;

while (validEpoch == false)

    fprintf('Select one of the five epochs for analysis... \n');
    fprintf('\t1 = preApnea (chill + free-breathe)\n');
    fprintf('\t2 = transition 1 (between preApnea and apnea)\n');
    fprintf('\t3 = apnea (breath-hold)\n');
    fprintf('\t4 = transition 2 (between apnea and postApnea)\n');
    fprintf('\t5 = postApnea (recovery)\n');
    epochChoice = str2double(input('Choose epoch for analysis: ','s'));
    
    switch(epochChoice)
       
        case 1      % PRE-APNEA
            
            if (exist('E1_KINEMATICS','var'))
                confText = 'Epoch analysis results found! Continue? (y/n): ';
                confStr = lower(input(confText,'s'));
                switch(confStr)
                    case 'y'
                        analysisRange = EPOCH.preApnea;
                        thisEpoch = E1_KINEMATICS;
                        validEpoch = true; 
                        if (exist('E1_CC','var'))
                            thisCC = E1_CC;
                            validCC = true;
                        end
                        if (exist('E1_KE','var'))
                           thisKE = E1_KE;
                           validKE = true;
                        else 
                            validKE = false;
                        end
                        
                    otherwise
                        fprintf('Skipping re-analysis. Choose a new epoch!\n');
                        validEpoch = false;
                end                
            else
                fprintf('EPOCH not found. Analyze EPOCH first or choose a new epoch!\n');
                validEpoch = false;       
            end
            
        case 2      % TRANSITION 1
            
            if (exist('E2_KINEMATICS','var'))
                confText = 'Epoch analysis results found! Continue? (y/n): ';
                confStr = lower(input(confText,'s'));
                switch(confStr)
                    case 'y'
                        analysisRange = EPOCH.t1;
                        thisEpoch = E2_KINEMATICS;
                        validEpoch = true;   
                        if (exist('E2_CC','var'))
                            thisCC = E2_CC;
                            validCC = true;  
                        end
                        if (exist('E2_KE','var'))
                           thisKE = E2_KE;
                           validKE = true;
                        else 
                            validKE = false;                           
                        end                        
                    otherwise
                        fprintf('Skipping re-analysis. Choose a new epoch!\n');
                        validEpoch = false;
                end                
            else
                fprintf('EPOCH not found. Analyze EPOCH first or choose a new epoch!\n');
                validEpoch = false;              
            end
            
            
        case 3      % APNEA
            
            if (exist('E3_KINEMATICS','var'))
                confText = 'Epoch analysis results found! Continue? (y/n): ';
                confStr = lower(input(confText,'s'));
                switch(confStr)
                    case 'y'
                        analysisRange = EPOCH.apnea;
                        thisEpoch = E3_KINEMATICS;
                        validEpoch = true;  
                        if (exist('E3_CC','var'))
                            thisCC = E3_CC;
                            validCC = true;                        
                        end
                        if (exist('E3_KE','var'))
                           thisKE = E3_KE;
                           validKE = true;
                        else
                            validKE = false;
                        end                        
                        
                    otherwise
                        fprintf('Choose a new epoch!\n');
                        validEpoch = false;
                end                
            else
                fprintf('EPOCH not found. Analyze EPOCH first or choose a new epoch!\n');
                validEpoch = false;                     
            end
            
        case 4      % TRANSITION 2
            
            if (exist('E4_KINEMATICS','var'))
                confText = 'Epoch analysis results found! Continue? (y/n): ';
                confStr = lower(input(confText,'s'));
                switch(confStr)
                    case 'y'
                        analysisRange = EPOCH.t2;
                        thisEpoch = E4_KINEMATICS;
                        validEpoch = true;
                        if (exist('E4_CC','var'))
                            thisCC = E4_CC;
                            validCC = true;
                        end
                        if (exist('E4_KE','var'))
                           thisKE = E4_KE;
                           validKE = true;
                        else 
                            validKE = false;                           
                        end                          
                        
                    otherwise
                        fprintf('Choose a new epoch!\n');
                        validEpoch = false;
                end                
            else
                fprintf('EPOCH not found. Analyze EPOCH first or choose a new epoch!\n');
                validEpoch = false;               
            end
            
        case 5      % POST-APNEA
            
            if (exist('E5_KINEMATICS','var'))
                confText = 'Epoch analysis results found! Continue? (y/n): ';
                confStr = lower(input(confText,'s'));
                switch(confStr)
                    case 'y'
                        analysisRange = EPOCH.postApnea;
                        thisEpoch = E5_KINEMATICS;
                        validEpoch = true;    
                        if (exist('E5_CC','var'))
                            thisCC = E5_CC;
                            validCC = true;   
                        end
                        if (exist('E5_KE','var'))
                           thisKE = E5_KE;
                           validKE = true;
                        else 
                            validKE = false;
                        end                                
                    otherwise
                        fprintf('Choose a new epoch!\n');
                        validEpoch = false;
                end                
            else
                fprintf('EPOCH not found. Analyze EPOCH first or choose a new epoch!\n');
                validEpoch = false;               
            end
            
        otherwise
            fprintf('That was not a valid epoch choice. Try again.\n');
    end    
end

%% plot some stuff for quickly IDing vent events
% 
% epochTimeRange = timerange(thisEpoch.Time(1), thisEpoch.Time(end));
% 
% figure('Color', white, 'Position', [0 500 500 500] );
% s1 = subplot(211);
% plot(KINEMATICS.Time(epochTimeRange), KINEMATICS.odba(epochTimeRange));
% yPos = s1.YLim(2) - ((s1.YLim(2) - s1.YLim(1)) / 10);
% hold on;
% ftPlotCues(CUE.Time(epochTimeRange), CUE.type(epochTimeRange), yPos);
% hold off;
% grid;
% s2 = subplot(212);
% plot(thisEpoch.Time, thisEpoch.consensusHr);
% grid;
% linkaxes([s1 s2],'x');
% 
% closedTest = lower(input('Closed this window before proceeding (y/n) [y]: ','s'));


%% make band-pass kinematics for this region

kFs = 100;

Time = thisEpoch.Time;
kT_s = KINEMATICS.time_s(analysisRange);

if (validKE)
   fprintf('Skipping band-pass filtering of epoch kinematic signals.\n');
   bpAx = thisKE.bpAx;
   bpAy = thisKE.bpAy;
   bpAz = thisKE.bpAz;
   bpGx = thisKE.bpGx;
   bpGy = thisKE.bpGy;
   bpGz = thisKE.bpGz;
else
    bandpassFiltOrder   = 6;
    kLowCut             = 5.625;    % low end for max SCG & GCG energy
    kHighCut            = 22.5;     % hi end for max SCG & GCG energy
    kNyquistFs          = 50;

    bandPassHrKinematicFilter = designfilt('bandpassiir', ...
        'FilterOrder',bandpassFiltOrder, ...
        'HalfPowerFrequency1', kLowCut,'HalfPowerFrequency2', kHighCut, ...
        'SampleRate', kNyquistFs);

    fprintf('Band-pass butterworth filter for seismocardiography energies:\n');
    fprintf('\tFilter order: %d\n', bandpassFiltOrder);
    fprintf('\tLow cut frequency: %d Hz\n', kLowCut);
    fprintf('\tHigh cut frequency: %d Hz\n', kHighCut);

    bpAx   = filtfilt(bandPassHrKinematicFilter', KINEMATICS.ax(analysisRange));
    bpAy   = filtfilt(bandPassHrKinematicFilter', KINEMATICS.ay(analysisRange));
    bpAz   = filtfilt(bandPassHrKinematicFilter', KINEMATICS.az(analysisRange));
    bpGx   = filtfilt(bandPassHrKinematicFilter', KINEMATICS.gx(analysisRange));
    bpGy   = filtfilt(bandPassHrKinematicFilter', KINEMATICS.gy(analysisRange));
    bpGz   = filtfilt(bandPassHrKinematicFilter', KINEMATICS.gz(analysisRange)); 
    
    thisKE = timetable(Time, kT_s, bpAx, bpAy, bpAz, bpGx, bpGy, bpGz, ...
    'VariableNames',{'time_s','bpAx','bpAy','bpAz','bpGx','bpGy','bpGz'} );

end

kT = Time;
clear Time;

%% check for pre-existing CC timetable and clean NATs before proceeding...

if (exist('thisCC','var'))
    
    cardiacCyclesDefined = true;

    clear natPositive;
    clear natIndex;
    clear numNats;
    
    fprintf('Checking for NATs in thisCC before proceeding...\n');
    
    q = 1;
    for n = 1:height(thisCC)
       if ( isnat(thisCC.Time(n)) )
            fprintf('Found a NAT in thisCC row %d. Adding to natIndex...\n', ...
                n);
            natIndex(q) = n;
            q = q + 1;
       end
    end
    
    if (~exist('natIndex','var'))
        
        fprintf('\tNo NATs detected in thisCC. Proceeding!\n')
        natIndex = 0;
        
    else

        numNats = numel(natIndex);

        if ( numNats > 1 )

            fprintf('NATs are consecutive. Dropping last %d rows.\n', ...
                numel(natIndex) );
            newEnd = natIndex(1) - 1;
            Time = thisCC.Time(1:newEnd);
            cueType = thisCC.cueType(1:newEnd);
            thisCC = timetable(Time,cueType,'VariableNames',{'cueType'});
            thisCC = sortrows(thisCC);

        elseif (numNats == 1)

            fprintf('NAT last entry is NAT. Fixing...\n');
            newEnd = natIndex(1) - 1;
            Time = thisCC.Time(1:newEnd);
            cueType = thisCC.cueType(1:newEnd);
            thisCC = timetable(Time,cueType,'VariableNames',{'cueType'});
            thisCC = sortrows(thisCC);        

        elseif (numNats == 0)

            fprintf('\tNo NATs in thisCC! Proceeding...\n');

        else
            fprintf('NATs are in nonconsecutive order. Handle manually.\n');
            return;
        end
        
    end
    
else
    
    cardiacCyclesDefined = false;

end

%% now go to work...

if (cardiacCyclesDefined)

    fprintf('thisCC timetable is already in place!\n');
    startThereTxt = 'Do you want to start at end of CC Time (y/n) [n]: ';
    startWhere = lower(input(startThereTxt,'s'));
    switch(startWhere)
        case 'n'

            startOverTxt = 'Do you want to clear CC and start over (y/n) [n]: ';
            startOver = lower(input(startOverTxt,'s'));
            switch(startOver)
                case 'y'
                    fprintf('Okay, nuking CC and starting at beginning.\n');
                    clear thisCC;
                    clear Time;
                    clear cueType;                    
                otherwise
                    fprintf('Okay, leaving CC in place.\n');
            end

        otherwise       % deal with pre-started CC
            startHere = thisCC.Time(end);
            secondsJump = seconds(startHere - datetime(thisEpoch.Time(1)));
            sampleJump = kFs * round(secondsJump, 2);
            fprintf('Okay, starting at sample %d | time: %s\n', ...
                sampleJump, startHere);

    end                 % end startWhere switch-case
end                     % end pre-existing case of CC timetable

% Create the epoch sample window for cardiac cycle annotation

timeWindow   = 5;                   % five seconds
sampleWindow = timeWindow * kFs;    % number of samples per window
lastSample = numel(thisKE.bpAx);    
    
if (exist('sampleJump','var'))
    startTime = round(sampleJump,0);
else 
    startTime = 1;
	timeWindow   = 5;                   % five seconds
    sampleWindow = timeWindow * kFs;    % number of samples per window
    lastSample = numel(thisKE.bpAx); 
end

endTime = startTime + sampleWindow;
if ( endTime > numel(kT) )
    endTime = numel(kT); 
    startTime = endTime - sampleWindow;
end

if (exist('thisCC','var'))
    useCC = true;
    ccTimes = thisCC.Time(thisCC.Time >= kT(startTime) & thisCC.Time <= kT(endTime) );
    ccCues = thisCC.cueType(thisCC.Time >= kT(startTime) & thisCC.Time <= kT(endTime) );
else
    useCC = false;
end

screenSize = get(0,'screensize');

figAnalysisRegion = figure('Color', white);
figSize = figAnalysisRegion.Position;
figSize(3) = screenSize(3);
figSize(4) = screenSize(4) * 0.7;
figAnalysisRegion.Position = figSize;  

p1a = subplot(411);
plot(thisEpoch.Time(startTime:endTime), bpAx(startTime:endTime), ...
    'Color', red);
hold on;
yPos = p1a.YLim(2) - ((p1a.YLim(2) - p1a.YLim(1)) / 10);
ftPlotCues(CUE.Time, CUE.type, yPos);
if (useCC)
    ftPlotCues(ccTimes, ccCues, yPos);
end
hold off;
xlabel('Time, local');
ylabel('m · s^-^2');
title('Band-passed Ax');
xlim([ kT(startTime) kT(endTime)]);
grid;

p1b = subplot(412);
plot(thisEpoch.Time(startTime:endTime), bpGy(startTime:endTime), ...
    'Color', green');
hold on;
yPos = p1b.YLim(2) - ((p1b.YLim(2) - p1b.YLim(1)) / 10);
ftPlotCues(CUE.Time, CUE.type, yPos);
if (useCC)
    ftPlotCues(ccTimes, ccCues, yPos);
end
hold off;
xlabel('Time, local');
ylabel('degrees · s^-^1');
xlim([ kT(startTime) kT(endTime)]);
title('Band-passed Gy');
grid;

p1c = subplot(413);
plot(thisEpoch.Time(startTime:endTime), ...    
    50 * cmddenoise(bpAx(startTime:endTime), 'sym5', 5, 's', 2), ...
    'Color', red);
hold on;
yPos = p1b.YLim(2) - ((p1b.YLim(2) - p1b.YLim(1)) / 10);
ftPlotCues(CUE.Time, CUE.type, yPos);
if (useCC)
    ftPlotCues(ccTimes, ccCues, yPos);
end
hold off;
xlabel('Time, local');
ylabel('20 · m · s^-^2');
xlim([ kT(startTime) kT(endTime)]);
title('bpAx w/ Daubechies least-asymmetric wavelet sym5 vanishing moments down to level 3.');
grid;

p1d = subplot(414);
plot(thisEpoch.Time(startTime:endTime), ...
    cmddenoise(bpGy(startTime:endTime), 'sym5', 4, 's', 2), ...
    'Color', green');
hold on;
yPos = p1b.YLim(2) - ((p1b.YLim(2) - p1b.YLim(1)) / 10);
ftPlotCues(CUE.Time, CUE.type, yPos);
if (useCC)
    ftPlotCues(ccTimes, ccCues, yPos);
end
hold off;
xlabel('Time, local');
ylabel('degrees · s^-^1');
xlim([ kT(startTime) kT(endTime)]);
title('bpGy w/ Daubechies least-asymmetric wavelet sym5 vanishing moments down to level 3.');
grid;

linkaxes([p1a p1b p1c p1d],'x');

% cardiac cycle annotation loop starts here

annotationComplete = false;

if (exist('thisCC','var'))
    Time = thisCC.Time';
    cueType = thisCC.cueType';
    i = length(thisCC.Time) + 1;
else
    i = 1;
    clear CC;
    clear ccTime; 
    clear ccCues;
end


while(~annotationComplete)
    
    fprintf('\n---\n');
    fprintf('\t[a] - AV open / MV close  | [c] - AV close / MV open\n');
    fprintf('\t[v] - ventilation event   | [x] - exclusion (auto-ends at next cue)\n');
    fprintf('\t[u] - undo last cue       | [r] - remove cue(s) within L-R range\n');
    fprintf('\t[b] - back five seconds   | [f] - forward five seconds\n');
    fprintf('\t[,] - back by one second  | [.] - forward by one second\n');
    fprintf('\t[<] - go to epoch start   | [>] - go to epoch end\n');
    fprintf('\t[j] - jump to (hh:mm:ss)  | [d] - clean & dump thisCC to console\n');
    fprintf('\t[s] - save & end annotate | [q] - quit w/o save\n');
    
    
    [x,~,button] = ginput(1);
    
    xAxis = gca;
    
    switch(button)
        case 'a'                                    % AO - aortic valve open
            Time(i) = num2ruler(x, xAxis.XAxis);
            cueType(i) = {'ao'};
            
            fprintf('AO -> time: %s | sample: %d\n', Time(i), kT_s(i));
            fprintf('Remaining: %s\n', duration(thisEpoch.Time(end) - Time(end)))
            thisCC = timetable(Time',cueType','VariableNames',{'cueType'});
            thisCC = sortrows(thisCC);
            i = i + 1;
            adjustTime = false;
            updatePlot = true;
            annotationComplete = false;
        case 'c'                                    % AC - aortic valve close
            Time(i) = num2ruler(x, xAxis.XAxis);
            cueType(i) = {'ac'};
            fprintf('AC: %s\n', Time(i));
            fprintf('Remaining: %s\n', duration(thisEpoch.Time(end) - Time(end)))
            thisCC = timetable(Time',cueType','VariableNames',{'cueType'});
            thisCC = sortrows(thisCC);
            i = i + 1;
            adjustTime = false;
            updatePlot = true;
            annotationComplete = false;
        case 'x'                                    % x - exclude marker
            Time(i) = num2ruler(x, xAxis.XAxis);
            cueType(i) = {'x'};
            fprintf('Exclude start: %s\n', Time(i));
            fprintf('Remaining: %s\n', duration(thisEpoch.Time(end) - Time(end)))
            thisCC = timetable(Time',cueType','VariableNames',{'cueType'});
            thisCC = sortrows(thisCC);
            i = i + 1;
            fprintf('Marking start of an exclusion area.\n');
            adjustTime = false;
            updatePlot = true;
            annotationComplete = false;
        case 'v'                                    % v - ventilation
            Time(i) = num2ruler(x, xAxis.XAxis);
            cueType(i) = {'v'};
            fprintf('V -> time: %s | sample: %d\n', Time(i), kT_s(i));
            fprintf('Remaining: %s\n', duration(thisEpoch.Time(end) - Time(end)))
            thisCC = timetable(Time',cueType','VariableNames',{'cueType'});
            thisCC = sortrows(thisCC);
            i = i + 1;
            adjustTime = false;
            updatePlot = true;
            annotationComplete = false;  
        case 'r'                                    % r - remove a cue
            fprintf('Click points to left then right of cue to remove...\n');
            clear rCoords;
            clear rTimes;
            clear rRange;
            clear removeRow;
            [rCoords,~,~] = ginput(2);
            rTimes = num2ruler(rCoords, xAxis.XAxis);
            rRange = timerange(rTimes(1), rTimes(2));
            rTargets = thisCC.Time(rRange);
            thisCC(rTargets,:) = [];
            Time = thisCC.Time';
            cueType = thisCC.cueType';
            adjustTime = false;
            updatePlot = true;
            annotationComplete = false;             
        case 'u'                                    % undo last add to thisCC
            Time = Time(1:end-1);
            cueType = cueType(1:end-1);
            i = i - 1;
            thisCC = timetable(Time',cueType','VariableNames',{'cueType'});
            thisCC = sortrows(thisCC);
            i = i + 1; 
        case 'j'
            seekDate = string(datetime(thisCC.Time(1),'Format','dd-MMM-yyyy'));
            jumpTxt = 'Enter jump time as hh:mm:ss: ';
            seekTime = input(jumpTxt,'s');
            jumpTo = sprintf('%s %s', seekDate, seekTime);
            fprintf('Jumping to %s\n', jumpTo);
            sampleJump = ...
                round(seconds(datetime(jumpTo) - thisCC.Time(1)),2) * 100;
            startTime = sampleJump;
            adjustTime = true;
            updatePlot = true;
            annotationComplete = false;
            
        case 'f'                                    % next 5 second window
            fprintf('[f]oward to the next sample window...\n');
            timeJump = seconds(5);
            startTime = startTime + sampleWindow;
            adjustTime = true;
            updatePlot = true;
            annotationComplete = false;
        case 'b'                                    % back 5 seconds
            fprintf('[b]ack to the previous sample window...\n');
            timeJump = seconds(-5);
            startTime = startTime - sampleWindow;
            adjustTime = true;
            updatePlot = true;
            annotationComplete = false;
        case '.'                                    % move +1 second window
            fprintf('[.] advancing the window by one second...\n');
            timeJump = seconds(1);
            startTime = startTime + 100;
            adjustTime = true;
            updatePlot = true;
            annotationComplete = false;
        case ','                                    % move -1 second window
            fprintf('[,] rewinding the window by one second...\n');
            timeJump = seconds(-1);
            startTime = startTime - 100;
            adjustTime = true;
            updatePlot = true;
            annotationComplete = false; 
        case '<'
            fprintf('[<] jumping to the beginning of the epoch...\n');
            startTime = 1;
            adjustTime = true;
            updatePlot = true;
            annotationComplete = false;
        case '>'
            fprintf('[>] jumping to the end of the epoch...\n');
            startTime = lastSample - sampleWindow + 1;
            adjustTime = true;
            updatePlot = true;
            annotationComplete = false;

        case 'd'
            
            clear natPositive;
            clear natIndex;
            clear numNats;

            fprintf('Checking for NATs in thisCC before proceeding...\n');

            q = 1;
            for n = 1:height(thisCC)
               if ( isnat(thisCC.Time(n)) )
                    fprintf('\tFound a NAT in thisCC row %d. Adding to natIndex...\n', ...
                        n);
                    natIndex(q) = n;
                    q = q + 1;
               end
            end
    
            if (~exist('natIndex','var'))

                fprintf('\tno NATs detected in thisCC. Proceeding!\n')
                natIndex = 0;

            else

                numNats = numel(natIndex);

                if ( numNats > 1 )

                    fprintf('NATs are consecutive. Dropping last %d rows.\n', ...
                        numel(natIndex) );
                    newEnd = natIndex(1) - 1;
                    Time = thisCC.Time(1:newEnd);
                    cueType = thisCC.cueType(1:newEnd);
                    thisCC = timetable(Time,cueType,'VariableNames',{'cueType'});
                    thisCC = sortrows(thisCC);

                elseif (numNats == 1)

                    fprintf('NAT last entry is NAT. Fixing...\n');
                    newEnd = natIndex(1) - 1;
                    Time = thisCC.Time(1:newEnd);
                    cueType = thisCC.cueType(1:newEnd);
                    thisCC = timetable(Time,cueType,'VariableNames',{'cueType'});
                    thisCC = sortrows(thisCC);        

                elseif (numNats == 0)

                    fprintf('\tNo NATs found in thisCC!\n');

                else
                    fprintf('NATs are in nonconsecutive order. Handle manually.\n');
                    return;
                end
        
            end   

            fprintf('Pruning thisCC timetable duplicates...\n');
            [thisCC, ia, ic] = unique(thisCC);
            sizeIA = numel(ia);
            sizeIC = numel(ic);
            if (sizeIA - sizeIC == 0)
                fprintf('\tNo duplicate entries found!\n');
            else
                fprintf('\tFixed one or more duplicate entries in thisCC.\n');
            end
            
            fprintf('Dumping contents of CC to console...\n\n');
            disp(thisCC);
            adjustTime = false;
            updatePlot = false;
            annotationComplete = false;            
            
        case 's'

            fprintf('Checking for sneaky NATs before saving...\n');

            clear natPositive;
            clear natIndex;
            clear numNats;


            fprintf('Checking for NATs in thisCC before proceeding...\n');

            q = 1;
            for n = 1:height(thisCC)
               if ( isnat(thisCC.Time(n)) )
                    fprintf('\tFound a NAT in thisCC row %d. Adding to natIndex...\n', ...
                        n);
                    natIndex(q) = n;
                    q = q + 1;
               end
            end
    
            if (~exist('natIndex','var'))

                fprintf('\tno NATs detected in thisCC. Proceeding!\n')
                natIndex = 0;

            else

                numNats = numel(natIndex);

                if ( numNats > 1 )

                    fprintf('NATs are consecutive. Dropping last %d rows.\n', ...
                        numel(natIndex) );
                    newEnd = natIndex(1) - 1;
                    Time = thisCC.Time(1:newEnd);
                    cueType = thisCC.cueType(1:newEnd);
                    thisCC = timetable(Time,cueType,'VariableNames',{'cueType'});
                    thisCC = sortrows(thisCC);

                elseif (numNats == 1)

                    fprintf('NAT last entry is NAT. Fixing...\n');
                    newEnd = natIndex(1) - 1;
                    Time = thisCC.Time(1:newEnd);
                    cueType = thisCC.cueType(1:newEnd);
                    thisCC = timetable(Time,cueType,'VariableNames',{'cueType'});
                    thisCC = sortrows(thisCC);        

                elseif (numNats == 0)

                    fprintf('\tNo NATs found in thisCC!\n');

                else
                    fprintf('NATs are in nonconsecutive order. Handle manually.\n');
                    return;
                end
        
            end   

            fprintf('Pruning thisCC timetable duplicates...\n');
            [thisCC, ia, ic] = unique(thisCC);
            sizeIA = numel(ia);
            sizeIC = numel(ic);
            if (sizeIA - sizeIC == 0)
                fprintf('\tNo duplicate entries found!\n');
            else
                fprintf('\tFixed one or more duplicate entries in thisCC.\n');
            end
            
    
            fprintf('Saving cardiac cycle cues in E*_CC...\n');
            
            switch(epochChoice)
                
                case 1
                    if (exist('E1_CC','var'))
                        confTxt = 'E1_CC exists. Overwrite (y/n) [n]: ';
                        conf = lower(input(confTxt,'s'));
                        switch(conf)
                            case 'y'
                                E1_CC = thisCC;
                            otherwise
                                fprintf('Skipping E*_CC = thisCC. Do manually.\n');
                        end        
                    else
                        E1_CC = thisCC;
                    end
                    
                case 2
                    if (exist('E2_CC','var'))
                        confTxt = 'E2_CC exists. Overwrite (y/n) [n]: ';
                        conf = lower(input(confTxt,'s'));
                        switch(conf)
                            case 'y'
                                E2_CC = thisCC;
                            otherwise
                                fprintf('Skipping E*_CC = CC. Do manually.\n');
                        end        
                    else
                        E2_CC = thisCC;
                    end                    
                    
                case 3
                    if (exist('E3_CC','var'))
                        confTxt = 'E3_CC exists. Overwrite (y/n) [n]: ';
                        conf = lower(input(confTxt,'s'));
                        switch(conf)
                            case 'y'
                                E3_CC = thisCC;
                            otherwise
                                fprintf('Skipping E*_CC = CC. Do manually.\n');
                        end        
                    else
                        E3_CC = thisCC;
                    end                    
                    
                case 4 
                    if (exist('E4_CC','var'))
                        confTxt = 'E4_CC exists. Overwrite (y/n) [n]: ';
                        conf = lower(input(confTxt,'s'));
                        switch(conf)
                            case 'y'
                                E4_CC = thisCC;
                            otherwise
                                fprintf('Skipping E*_CC = CC. Do manually.\n');
                        end        
                    else
                        E4_CC = thisCC;
                    end
                    
                case 5
                    
                    if (exist('E5_CC','var'))
                        confTxt = 'E5_CC exists. Overwrite (y/n) [n]: ';
                        conf = lower(input(confTxt,'s'));
                        switch(conf)
                            case 'y'
                                E5_CC = thisCC;
                            otherwise
                                fprintf('Skipping E*_CC = CC. Do manually.\n');
                        end        
                    else
                        E5_CC = thisCC;
                    end                    
                    
            end
                        
            annotationComplete = true;
            adjustTime = false;
            updatePlot = false;
            cardiacCyclesDefined = true;
            
            if (cardiacCyclesDefined)

                fprintf('Appending E*_CC timetable to RAW file: %s\n', ...
                    fullRawFileName );

                if (exist(fullRawFileName, 'file'))

                    fprintf('RAW file exists... ');

                    switch(epochChoice)
                        case 1
                            prompt = 'OK to append new E1_CC timetable? y/n [n]: ';
                        case 2
                            prompt = 'OK to append new E2_CC timetable? y/n [n]: ';
                        case 3
                            prompt = 'OK to append new E3_CC timetable? y/n [n]: ';
                        case 4
                            prompt = 'OK to append new E4_CC timetable? y/n [n]: ';
                        case 5
                            prompt = 'OK to append new E5_CC timetable? y/n [n]: ';
                    end

                    str = input(prompt,'s');
                    if isempty(str)
                        str = 'n';
                    end

                    if (strcmp(str, 'y'))

                        switch(epochChoice)
                            case 1
                                E1_KE = thisKE;
                                fprintf('\nAppending E1_CC timetable to RAW file...\n');
                                save(fullRawFileName, 'E1_CC', 'E1_KE', '-append');  
                                fprintf('E1_CC added to RAW.\n');
                            case 2
                                E2_KE = thisKE;
                                fprintf('\nAppending E2_CC timetable to RAW file...\n');
                                save(fullRawFileName, 'E2_CC', 'E2_KE', '-append');               
                                fprintf('E2_CC added to RAW.\n');
                            case 3
                                E3_KE = thisKE;
                                fprintf('\nAppending E3_CC timetable to RAW file...\n');
                                save(fullRawFileName, 'E3_CC', 'E3_KE', '-append');               
                                fprintf('E3_CC added to RAW.\n');
                            case 4
                                E4_KE = thisKE;
                                fprintf('\nAppending E4_CC timetable to RAW file...\n');
                                save(fullRawFileName, 'E4_CC', 'E4_KE', '-append');
                                fprintf('E4_CC added to RAW.\n');
                            case 5
                                E5_KE = thisKE;
                                fprintf('\nAppending E5_CC timetable to RAW file...\n');
                                save(fullRawFileName, 'E5_CC', 'E5_KE', '-append');
                                fprintf('E5_CC added to RAW.\n');
                        end

                        fprintf('Ending annotation. Moving on the energy calculations.\n');

                    end 

                else

                    fprintf('Raw file does not exist... this should never happen... Bail!\n');
                    return;

                end

            end
            
        case 'q'
            fprintf('Caught a [q]uit. Save E*_CC = thisCC manually. Exiting!\n');
            annotationComplete = true;
            updatePlot = false;
            adjustTime = false;
        otherwise
            fprintf('Invalid input. Continuing...\n');
            annotationComplete = false;
            updatePlot = false;
    
    end         % end switch(button)
    
    if (adjustTime)
        
        endTime = startTime + sampleWindow;
        if ( endTime > numel(kT) )
            endTime = numel(kT); 
            startTime = endTime - 500;
        end
        updatePlot = true;
    
    end
    
    if (updatePlot)
        
        ccTimes = ...
            thisCC.Time(thisCC.Time >= kT(startTime) & thisCC.Time <= kT(endTime) );
        ccCues = ...
            thisCC.cueType(thisCC.Time >= kT(startTime) & thisCC.Time <= kT(endTime) );
        
        figure(figAnalysisRegion);    
        p1a = subplot(411);
        plot(thisEpoch.Time(startTime:endTime), bpAx(startTime:endTime), ...
            'Color', red);
        hold on;
        yPos = p1a.YLim(2) - ((p1a.YLim(2) - p1a.YLim(1)) / 10);
        ftPlotCues(CUE.Time, CUE.type, yPos);
        ftPlotCues(ccTimes, ccCues, yPos);
        hold off;
        xlabel('Time, local');
        ylabel('m · s^-^2');
        title('Band-passed Ax');
        xlim([ kT(startTime) kT(endTime)]);
        grid;

        p1b = subplot(412);
        plot(thisEpoch.Time(startTime:endTime), bpGy(startTime:endTime), ...
            'Color', green');
        hold on;
        yPos = p1b.YLim(2) - ((p1b.YLim(2) - p1b.YLim(1)) / 10);
        ftPlotCues(CUE.Time, CUE.type, yPos);
        ftPlotCues(ccTimes, ccCues, yPos);
        hold off;
        xlabel('Time, local');
        ylabel('degrees · s^-^1');
        xlim([ kT(startTime) kT(endTime)]);
        title('Band-passed Gy');
        grid;

        p1c = subplot(413);
        plot(thisEpoch.Time(startTime:endTime), ...
            50 * cmddenoise(bpAx(startTime:endTime), 'sym4', 5, 's', 2), ...
            'Color', red);            
        hold on;
        yPos = p1b.YLim(2) - ((p1b.YLim(2) - p1b.YLim(1)) / 10);
        ftPlotCues(CUE.Time, CUE.type, yPos);
        if (useCC)
            ftPlotCues(ccTimes, ccCues, yPos);
        end
        hold off;
        xlabel('Time, local');
        ylabel('20 · m · s^-^2');
        xlim([ kT(startTime) kT(endTime)]);
        title('bpAx w/ Daubechies least-asymmetric wavelet sym5 vanishing moments down to level 3.');
        grid;

        p1d = subplot(414);
        plot(thisEpoch.Time(startTime:endTime), ...
            cmddenoise(bpGy(startTime:endTime), 'sym5', 5, 's', 2), ...
            'Color', green');
        hold on;
        yPos = p1b.YLim(2) - ((p1b.YLim(2) - p1b.YLim(1)) / 10);
        ftPlotCues(CUE.Time, CUE.type, yPos);
        if (useCC)
            ftPlotCues(ccTimes, ccCues, yPos);
        end
        hold off;
        xlabel('Time, local');
        ylabel('°/s');
        xlim([ kT(startTime) kT(endTime)]);
        title('bpGy w/ Daubechies least-asymmetric wavelet sym5 vanishing moments down to level 3.');
        grid;
        
        linkaxes([p1a p1b p1c p1d],'x');    
    
        adjustTime = false;
    
    end
    
end             % end annotationComplete test




%% append CC timetable to RAW file

% if (cardiacCyclesDefined)
%    
%     fprintf('Appending E*_CC timetable to RAW file: %s\n', fullRawFileName );
% 
%     if (exist(fullRawFileName, 'file'))
% 
%         fprintf('RAW file exists... ');
% 
%         switch(epochChoice)
%             case 1
%                 prompt = 'OK to append new E1_CC timetable? y/n [n]: ';
%             case 2
%                 prompt = 'OK to append new E2_CC timetable? y/n [n]: ';
%             case 3
%                 prompt = 'OK to append new E3_CC timetable? y/n [n]: ';
%             case 4
%                 prompt = 'OK to append new E4_CC timetable? y/n [n]: ';
%             case 5
%                 prompt = 'OK to append new E5_CC timetable? y/n [n]: ';
%         end
%         
%         str = input(prompt,'s');
%         if isempty(str)
%             str = 'n';
%         end
% 
%         if (strcmp(str, 'y'))
%             
%             switch(epochChoice)
%                 case 1
%                     E1_KE = thisKE;
%                     fprintf('\nAppending E1_CC timetable to RAW file...\n');
%                     save(fullRawFileName, 'E1_CC', 'E1_KE', '-append');  
%                     fprintf('E1_CC added to RAW.\n');
%                 case 2
%                     E2_KE = thisKE;
%                     fprintf('\nAppending E2_CC timetable to RAW file...\n');
%                     save(fullRawFileName, 'E2_CC', 'E2_KE', '-append');               
%                     fprintf('E2_CC added to RAW.\n');
%                 case 3
%                     E3_KE = thisKE;
%                     fprintf('\nAppending E3_CC timetable to RAW file...\n');
%                     save(fullRawFileName, 'E3_CC', 'E3_KE', '-append');               
%                     fprintf('E3_CC added to RAW.\n');
%                 case 4
%                     E4_KE = thisKE;
%                     fprintf('\nAppending E4_CC timetable to RAW file...\n');
%                     save(fullRawFileName, 'E4_CC', 'E4_KE', '-append');
%                     fprintf('E4_CC added to RAW.\n');
%                 case 5
%                     E5_KE = thisKE;
%                     fprintf('\nAppending E5_CC timetable to RAW file...\n');
%                     save(fullRawFileName, 'E5_CC', 'E5_KE', '-append');
%                     fprintf('E5_CC added to RAW.\n');
%             end
%                                     
%             fprintf('Ending annotation. Moving on the energy calculations.\n');
%             
%         end 
% 
%     else
% 
%         fprintf('Raw file does not exist... this should never happen... Bail!\n');
%         return;
% 
%     end
%         
% end



%% energy calculation

doCardioEnergy = lower(input('Do cardiac energy analysis (y/n) [n]: ','s'));
switch(doCardioEnergy)
    case 'y'
        doCardioEnergy = true;
    otherwise
        doCardioEnergy = false;
end

if (doCardioEnergy)

    numAllRows = height(thisCC);
    ccRows = (string(thisCC.cueType) == 'ao');
    numCCRows =  numel(ccRows(ccRows == 1));

    figV = figure('Color', white, 'Position', [0 0 500 500] );
    figKE = figure('Color', white, 'Position', [0 600 500 500] );

    prePeak = seconds(0.200);

    clear linearCE;
    clear rotationalCE;
    clear T;

    index = 0;

    clear ccAxyz ccVxyz ccGxyz cc_iLinKE cc_iRotKE T

    ccAxyz(1:150,1:3,1:numCCRows) = zeros;
    ccVxyz(1:150,1:3,1:numCCRows) = zeros;
    ccGxyz(1:150,1:3,1:numCCRows) = zeros;
    cc_iLinKE(1:150,1:3,1:numCCRows) = zeros;
    cc_iRotKE(1:150,1:3,1:numCCRows) = zeros;

    T(1:150,1,1:numCCRows) = datetime;

    for rowNum = 1:numAllRows-1
        
        fprintf('Evaluating row number %d\n', rowNum);
        
        thisCue = string(thisCC.cueType(rowNum));
                        
        switch(thisCue)
            
            case 'ao'       % start cue for a good cardiac cycle

                fprintf('\tthisCue: AO\n');
                
                nextCue = string(thisCC.cueType(rowNum + 1));

                switch(nextCue)
                    
                    case 'ao'
                        
                        fprintf('\tnextCue: AO. Valid couplet!\n');
                        
                        validCouplet = true;
                        t1 = thisCC.Time(rowNum) - prePeak;
                        t2 = thisCC.Time(rowNum + 1) - prePeak;
                        
                    case 'x'
                        
                        fprintf('\tnextCue: X. Valid couplet!\n');
                        
                        validCouplet = true;
                        t1 = thisCC.Time(rowNum) - prePeak;
                        t2 = thisCC.Time(rowNum + 1) - prePeak;
                        
                    case 'v'
                        
                        fprintf('\tnextCue: V. Invalid couplet...\n');
                        fprintf('\t\ttrying rowNum + 2...\n');
                        
                        nextCue = string(thisCC.cueType(rowNum + 2));
                        
                        if ( contains(nextCue,'ao') || ...
                                contains(nextCue, 'x') )
                            
                            fprintf('\tnextCue + 1: AO or X. Valid couplet!\n');
                            
                            validCouplet = true;
                            t1 = thisCC.Time(rowNum) - prePeak;
                            t2 = thisCC.Time(rowNum + 2) - prePeak;                         
                            
                        else 
                            
                            validCouplet = false;
                        
                        end
                        
                end
                
            case 'x'        % start of an exclusion zone
                
                fprintf('\tthisCue: X. Not a couplet.\n');
                validCouplet = false;
                
            case 'v'        % start of a ventilation cue
                
                fprintf('\tthisCue: V. Not a couplet.\n');
                validCouplet = false;
                
            otherwise
                
                fprintf('Found an invalid testCue at row number %d.\n', ...
                    rowNum); 
                fprintf('Halting at / around line 1195.\n');
                return;
                
        end

        
        if (validCouplet)

            index = index + 1;

            fprintf('\tWorking on cardiac cycle %d at thisCC.rowNum %d...\n', ...
                index, rowNum);                
            
            fprintf('\ttestCue: %s | nextCue: %s\n', thisCue, nextCue);
            
            durationSeconds(index) = seconds(diff([t1 t2]));
                
            fprintf('\tCC(%d) duration: %1.3f seconds\n', ...
                index, durationSeconds(index) );
                
            t = timerange( t1, t2 );
                
            thisT = thisKE.Time(t);        

            T(:,1,index) = ...
                ( thisT(1):duration(seconds(0.01)):thisT(1)+seconds(1.49) )';

            Axyz = [ thisKE.bpAx(t) thisKE.bpAy(t) thisKE.bpAz(t) ];
            Gxyz = [ thisKE.bpGx(t) thisKE.bpGy(t) thisKE.bpGz(t) ];
            GxyzSI = deg2rad(Gxyz);
            
            % use size(Axyz) to make Vxyz, the dt integral of Axyz

            clear Vxyz;
            sizeVxyz = size(Axyz);
            Vxyz(sizeVxyz(1),sizeVxyz(2)) = zeros;    

            fprintf('\tsizeVxyz: %d\n', sizeVxyz);

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
            lenGxyz = length(Gxyz);

            % when a cardiac cycle is longer than 1.5 seconds, clip at 1.5 s

            if (lenVxyz >= 150)
                lenVxyz = 150;
            end

            if (lenGxyz >= 150)
                lenGxyz = 150;
            end

            ccAxyz(1:lenVxyz,1:3,index) = Axyz(1:lenVxyz,:);
            ccVxyz(1:lenVxyz,1:3,index) = Vxyz(1:lenVxyz,:);
            ccGxyz(1:lenGxyz,1:3,index) = GxyzSI(1:lenGxyz,:);

            figure(figV);

            pLinAx = subplot(331);
            hold on;
            plot(Axyz(:,1),'Color', [0.6 0.6 0.6] );
            hold off;
            
            pLinVx = subplot(332);
            hold on;
            plot(Vxyz(:,1),'Color', blue);
            hold off;

            pAngVx = subplot(333);
            hold on;
            plot(GxyzSI(:,1), 'Color', blue);
            hold off;

            pLinAy = subplot(334);
            hold on;
            plot(Axyz(:,2),'Color', [0.6 0.6 0.6] );
            hold off;            
            
            pLinVy = subplot(335);
            hold on;
            plot(Vxyz(:,2),'Color', red);
            hold off;

            pAngVy = subplot(336);
            hold on;
            plot(GxyzSI(:,2), 'Color', red);
            hold off;

            pLinAz = subplot(337);
            hold on;
            plot(Axyz(:,3),'Color', [0.6 0.6 0.6] );
            hold off;            
            
            pLinVz = subplot(338);
            hold on;
            plot(Vxyz(:,3),'Color', goldenrod);
            hold off;

            pAngVz = subplot(339);
            hold on;
            plot(GxyzSI(:,3), 'Color', goldenrod);
            hold off;

            meanVxyz(index,1:3) = mean(Vxyz);
            meanGxyzSI(index,1:3) = mean(GxyzSI);

            % mass FaunaTag 111: 0.403 kg on 1 May & 1 Sept 2021
            massTag = 0.403;
            
            radiusX = 0.17/2;
            radiusY = 0.17/2;
            radiusZ = 0.023;
            
            Ix = massTag * radiusX ^ 2;
            Iy = massTag * radiusY ^ 2;
            Iz = massTag * radiusZ ^ 2;

%             % radius from dolphin heart center and tag sensor (m); estimate
%             rTagCenter = 0.05; 
% 
%             % moment of inertia constant
%             iTag = massTag * rTagCenter;

            keLin = 0.5 * massTag * ...
                ( Vxyz(:,1).^2 + Vxyz(:,2).^2 + Vxyz(:,3).^2 );

            
% NOTE: I have big doubts that using the Migeotte et al 2015 "mass of
% subject, i.e.: dolphin, is the right way to go for computing Ix Iy Iz.
% I've gone with an estimate for these moments of inertia using the mass
% of the tag times the square of the radius from the center of the tag to
% the edges of the suction cups, which are the anchors for any moments of
% inertia in the x and y rotational dimension. For Iz I've used the
% distance from the animal's body to the top of the tag as the radius,
% although thiss could be incorrect. Until I can implement and test against
% Wang et al 2016 linear regression method of computing Ix Iy and Iz, I
% think my method is the better way, with the caveat that it be
% characterized as an estimate of iRotKE. The thing that makes me the most
% wary of my method is the absurdity proposition that doubling the mass of
% the tag doubles the kinetic energy of cardiac measurements. This whole
% notion needs to be rethought, not just by me but by Migeotte, Yang, and
% all of us who are trying to use GCG kinetic energy estimates. Maybe using
% Power, i.e.:
% linP = F(t) · v(t)
% rotP = torque(t) · angular velocity(t), where torque(t) is:
%        ( I · moment of inertia for object · angular acceleration )   

%             keRot = 0.5 * ...
%                 ( (iTag * Gxyz(:,1).^2) + (iTag * Gxyz(:,1).^2) + ...
%                 (iTag * Gxyz(:,3).^2) );
            
            keRot = 0.5 * ...
                ( (Ix * GxyzSI(:,1).^2) + (Iy * GxyzSI(:,1).^2) + ...
                (Iz * GxyzSI(:,3).^2) );            

            for k = 1:length(Vxyz) - 1

                % calculate instantaneous linear & rotational kinetic energy

                iLinKE(k,:) = simps(0.5 * massTag * ( Vxyz(k:k+1,:).^2 ) );
                % iRotKE(k,:) = simps(0.5 * ( (iTag * Gxyz(k:k+1,:) ).^2 ));
                
%                 iRotKE(k,:) = simps(0.5 * ( Ix * Gxyz(k:k+1,1).^2 + ...
%                     Iy * Gxyz(k:k+1,2).^2 + Iz * Gxyz(k:k+1,3).^2 ) );

                iRotKE(k,1) = simps( 0.5 * (Ix * GxyzSI(k:k+1,1) .^ 2) );
                iRotKE(k,2) = simps( 0.5 * (Iy * GxyzSI(k:k+1,2) .^ 2) );
                iRotKE(k,3) = simps( 0.5 * (Iz * GxyzSI(k:k+1,3) .^ 2) );
                
            end

            % toss iLinKE and iRotKE into cc_iLinKE and cc_iRotKE

            lenILinKE = length(iLinKE);
            if ( length(iLinKE) >= 150 )
                lenILinKE = 150;
            end

            lenIRotKE = length(iRotKE);
            if (length(iRotKE) >= 150 )
                lenIRotKE = 150;
            end

            cc_iLinKE(1:lenILinKE,1:3,index) = iLinKE(1:lenILinKE,:);
            cc_iRotKE(1:lenIRotKE,1:3,index) = iRotKE(1:lenIRotKE,:);

            
            ccTime(index) = thisT(1);
            mean_iLinKE(index,1:3) = mean(iLinKE);
            mean_iRotKE(index,1:3) = mean(iRotKE);

            figure(figKE);
            pLinKE = subplot(211);
            hold on;
            plot(iLinKE, 'Color', red);
            hold off;
            pRotKE = subplot(212);
            hold on;
            plot(iRotKE, 'Color', green);
            hold off;            
            
        else
            
            fprintf('No valid couplet starting at row number %d.\n', rowNum);
        
        end

    end             % end for rowNum = 1:numCC-1

    figure(figV);
    pLinAx = subplot(331);
    hold on;
    xlabel('Time, milliseconds');
    ylabel('m · s^-2');
    title('Linear acceleration x-axis');
    grid;
    hold off;    
    pLinVx = subplot(332);
    hold on;
    xlabel('Time, milliseconds');
    ylabel('m · s^-1');
    title('Linear velocity x-axis');
    grid;
    hold off;
    pAngVx = subplot(333);
    hold on;
    xlabel('Time, milliseconds');
    ylabel('rad · s^-1');
    title('Angular velocity x-axis');
    grid;
    hold off;
    pLinAy = subplot(334);
    hold on;
    xlabel('Time, milliseconds');
    ylabel('m · s^-2');
    title('Linear acceleration y-axis');
    grid;
    hold off;     
    pLinVy = subplot(335);
    hold on;
    xlabel('Time, milliseconds');
    ylabel('m · s^-1');
    title('Linear velocity y-axis');
    grid;
    hold off;
    pAngVy = subplot(336);
    hold on;
    xlabel('Time, milliseconds');
    ylabel('rad · s^-1');
    title('Angular velocity y-axis');
    grid;
    hold off;
	pLinAz = subplot(337);
    hold on;
    xlabel('Time, milliseconds');
    ylabel('m · s^-2');
    title('Linear acceleration z-axis');
    grid;
    hold off; 
    pLinVz = subplot(338);
    hold on;
    xlabel('Time, milliseconds');
    ylabel('m · s^-1');
    title('Linear velocity z-axis');
    grid;
    hold off;
    pAngVz = subplot(339);
    hold on;
    xlabel('Time, milliseconds');
    ylabel('° · s^-1');
    title('Angular velocity z-axis');
    grid;
    hold off;
    linkaxes([pLinAx pLinAy pLinAz pLinVx pAngVx pLinVy pAngVy pLinVz pAngVz],'x');

    figure(figKE);
    pLinKE = subplot(211);
    hold on;
    grid;
    hold off;
    pRotKE = subplot(212);
    hold on;
    grid;
    hold off;
    linkaxes([pLinKE pRotKE],'x');

    mean_iKEtotal = ...
        ( sum(mean_iLinKE(:,1:3),2) + sum(mean_iRotKE(:,1:3),2) );

    if (epochChoice == 1 || epochChoice == 5)
        ventTimes = thisCC.Time(string(thisCC.cueType) == 'v');
    else
        ventTimes = thisCC.Time(1);     % set last breath at apnea onset
    end

    % plot instantaneous kinetic energy of Lin, Rot & Lin+Rot
    figCCEnergy = figure('Color', white, 'Position', [0 0 800 600] );
    
    s2a = subplot(311);
    plot(ccTime, sum(mean_iLinKE(:,1:3),2), '.', 'Color', red, ...
        'MarkerSize', 8);
    hold on;
    xline(ventTimes, '-', 'Color', red);
    hold off;
    xlabel('Time, local');
    ylabel('J · s');
    title('Cardiac cycle - instantaneous linear kinetic energy');
    grid;
    
    s2b = subplot(312);
    plot(ccTime, sum(mean_iRotKE(:,1:3),2), '.', 'Color', green, ...
        'MarkerSize', 8);
    hold on;
    xline(ventTimes, '-', 'Color', red);
    hold off;
    xlabel('Time, local');
    ylabel('J · s');
    title('Cardiac cycle - instantaneous rotational kinetic energy');
    grid;
    
    s2c = subplot(313);
    plot(ccTime, mean_iKEtotal, '.', 'Color', black, 'MarkerSize', 8);
    hold on;
    xline(ventTimes, '-', 'Color', red);
    hold off;
    xlabel('Time, local');
    ylabel('J · s');
    title('Cardiac cycle - instantaneous total kinetic energy');
    grid;
    
    linkaxes([s2a s2b s2c],'x');

    meanLinearAxyz = mean(ccAxyz(:,1:3,1:end),3);
    meanLinearVxyz = mean(ccVxyz(:,1:3,1:end),3);
    
    meanRotationalGxyz = mean(ccGxyz(:,1:3,1:end),3);

    meanCC_iLinKE = mean(cc_iLinKE(:,1:3,1:end),3);
    meanCC_iRotKE = mean(cc_iRotKE(:,1:3,1:end),3);

    cycleTime = 0.01:0.01:1.5;


    figCCKE = figure('Color', white, 'Position', [0 500 1200 800] );

    for z = 1:numCCRows

        nullTest = ccVxyz(:,:,z);

        if (nullTest == 0)  % ccVxyz is full of zeros, meaning invalid CC
            fprintf('Skipping cycle %d ...\n', z); 
        else
            fprintf('Processing cycle %d...\n', z); 
        end

        s3a = subplot(221);
        hold on;
        pCCVxyz = plot(cycleTime, ccVxyz(:,:,z), ':','Color', grey );
        hold off;

        s3b = subplot(222);
        hold on;
        pCCILinKE = plot(cycleTime, cc_iLinKE(:,:,z), ':', 'Color', grey );
        hold off;

        s3c = subplot(223);
        hold on;
        pCCGxyz = plot(cycleTime, ccGxyz(:,:,z), ':','Color', grey );
        hold off;    

        s3d = subplot(224);
        hold on;
        pCCIRotKE = plot(cycleTime, cc_iRotKE(:,:,z), ':', 'Color', grey );
        hold off;

    end

    % overplot mean linear and mean rotational velocities
    s3a = subplot(221);
    hold on;
    pMeanVx = plot(cycleTime, meanLinearVxyz(:,1), 'Color', ...
        blue, 'LineWidth', 3);
    pMeanVy = plot(cycleTime, meanLinearVxyz(:,2), 'Color', ...
        red, 'LineWidth', 3);
    pMeanVz = plot(cycleTime, meanLinearVxyz(:,3), 'Color', ...
        goldenrod, 'LineWidth', 3);
    hold off;

    s3b = subplot(222);
    hold on;
    pMeanLinKeX = plot(cycleTime, meanCC_iLinKE(:,1), 'Color', ...
        blue, 'LineWidth', 3);
    pMeanLinKeY = plot(cycleTime, meanCC_iLinKE(:,2), 'Color', ...
        red, 'LineWidth', 3);
    pMeanLinKeZ = plot(cycleTime, meanCC_iLinKE(:,3), 'Color', ...
        goldenrod, 'LineWidth', 3);
    hold off;

    s3c = subplot(223);
    hold on;
    pMeanGx = plot(cycleTime, meanRotationalGxyz(:,1), 'Color', ...
        blue, 'LineWidth', 3);
    pMeanGy = plot(cycleTime, meanRotationalGxyz(:,2), 'Color', ...
        red, 'LineWidth', 3);
    pMeanGz = plot(cycleTime, meanRotationalGxyz(:,3), 'Color', ...
        goldenrod, 'LineWidth', 3);
    hold off;

    s3d = subplot(224);
    hold on;
    pMeanRotKeX = plot(cycleTime, meanCC_iRotKE(:,1), 'Color', ...
        blue, 'LineWidth', 3);
    pMeanRotKeY = plot(cycleTime, meanCC_iRotKE(:,2), 'Color', ...
        red, 'LineWidth', 3);
    pMeanRotKeZ = plot(cycleTime, meanCC_iRotKE(:,3), 'Color', ...
        goldenrod, 'LineWidth', 3);
    hold off;

    % add in labels, titles and legends

    s3a = subplot(221);
    hold on;
    xlabel('Time, seconds');
    ylabel('m · s^-^1');
    title('Cardiac cycle linear velocity');
    grid;
    legend([ pMeanVx pMeanVy pMeanVz pCCVxyz(1,1) ], ...
        {'Vx','Vy','Vz','all cycles'});
    hold off;

    s3b = subplot(222);
    hold on;
    xlabel('Time, seconds');
    ylabel('J · s');
    title('Cardiac cycle linear kinetic energy');
    grid;
    legend([pMeanLinKeX pMeanLinKeY pMeanLinKeZ pCCILinKE(1,1) ], ...
        {'xLinKE','yLinKE','zLinKE','all cycles'});
    hold off;

    s3c = subplot(223);
    hold on;
    xlabel('Time, seconds');
    ylabel('degrees · s^-^1');
    title('Cardiac cycle angular velocity');
    grid;
    legend([pMeanGx pMeanGy pMeanGz pCCGxyz(1,1) ], ...
        {'Gx','Gy','Gz','all cycles'});
    hold off;

    s3d = subplot(224);
    hold on;
    xlabel('Time, seconds');
    ylabel('J · s');
    title('Cardiac cycle rotational kinetic energy');
    grid;
    legend([pMeanRotKeX pMeanRotKeY pMeanRotKeZ pCCIRotKE(1,1) ], ...
        {'xRotKE','yRotKE','zRotKE','all cycles'});
    hold off;

    linkaxes([s3a s3b s3c s3d],'x');

    % look at durationSeconds plot
    
    figure('Color', white, 'Position', [0 200 500 500] );
    plot(ccTime, durationSeconds, 'Color', blue);
    xlabel('Time, local');
    ylabel('Seconds');
    title('Cardiac cycle duration');
    grid;
    
    % create a total 
    % assemble some means for CC_VEL_KE

    cc_iTotalKE = cc_iLinKE + cc_iRotKE;
        
    mean_ccVxyz(1,1:3,1:numCCRows) = mean(ccVxyz);
    mean_ccGxyz = mean(ccGxyz);

    % sum_ccILinKE = sum(cc_iLinKE);
    % sum_ccIRotKE = sum(cc_iRotKE);
    % sum_ccITotalKE = sum(sum_ccILinKE) + sum(sum_ccIRotKE);

    for zz = 1:numCCRows 
        temp1(1:3,zz) = sum(cc_iLinKE(:,:,zz));
        temp2(1:3,zz) = sum(cc_iRotKE(:,:,zz));
    end

    % create sums of each cardiac cycles lin and rot KE along xyz axes
    sum_ccILinKE_xyz = temp1';
    sum_ccIRotKE_xyz = temp2';
    sum_ccITotalKE_xyz = sum_ccILinKE_xyz + sum_ccIRotKE_xyz;
    
    % create sums of x+y+z axes KE for each cardiac cycle
    sum_ccILinKE_all = sum(sum_ccILinKE_xyz, 2);
    sum_ccIRotKE_all = sum(sum_ccIRotKE_xyz, 2);
    sum_ccITotalKE_all = sum_ccILinKE_all + sum_ccIRotKE_all;

    mean_ccILinKE_xyz = mean(sum_ccILinKE_xyz);
    mean_ccIRotKE_xyz = mean(sum_ccIRotKE_xyz);
    mean_ccITotalKE_xyz = mean(sum_ccITotalKE_xyz);
    
    mean_ccILinKE_all = mean(sum_ccILinKE_all);
    mean_ccIRotKE_all = mean(sum_ccIRotKE_all);
    mean_ccITotalKE_all = mean(sum_ccITotalKE_all);

    % create CC_VEL_KE

    TempTime = Time;
    Time = T;

    % make a TIMETABLE with (150) x (1 or 1:3) x (number of CC rows)
    
    CC_VEL_KE = table(Time, ccVxyz, ccGxyz, cc_iLinKE, cc_iRotKE);

    Time = TempTime;

    % prep animalID, animalMass, and epochID for inclusion in TIMETABLE
    
    subjectID(1:numCCRows,1) = {TRIAL.SUBJECT.id};
    
    switch(TRIAL.SUBJECT.id)
        case '9FL3'                 % Lono
            subjectMass = 251.7;
        case '6JK5'                 % Kolohe
            subjectMass = 200.9;
        case '9ON6'                 % Noa
            subjectMass = 192.8;
        case '63HF'                 % Hoku
            subjectMass = 179.6;
        case '01L5'                 % Liho
            subjectMass = 163.7;
        case '83H1'                 % Hua
            subjectMass = 147.0;
        otherwise
            mass = input('Animal ID not found. Enter mass manually: ','s');
    end
    
    animalMass(1:numCCRows,1) = subjectMass;
    
    epochID(1:numCCRows,1) = epochChoice;
    
    iHR(1:numCCRows,1) = 60 ./ durationSeconds;
    
    tempTime = Time;
    
    newTime = Time(ccRows == 1);
    
    clear Time;
    
    Time = newTime;
    
    % make a TIMETABLE with (number of CC rows) x (1 or 1:3) sizing
    
    CC_ENERGY = timetable(Time', durationSeconds', iHR, ...
        sum_ccILinKE_xyz, sum_ccIRotKE_xyz, sum_ccITotalKE_xyz, ...
        sum_ccILinKE_all, sum_ccIRotKE_all, sum_ccITotalKE_all, ...
        subjectID, animalMass, epochID, ...
        'VariableNames', ...
        {'duration', 'iHR', 'iLinKE_xyz','iRotKE_xyz','iTotalKE_xyz', ...
        'iLinKE','iRotKE', 'iTotalKE', 'ID', 'mass', 'epoch'} );
    
    clear Time;
    
    Time = tempTime;
    
    % clean up CC_ENERGY

    CC_ENERGY = sortrows(CC_ENERGY);

    longDurations = CC_ENERGY.Time(CC_ENERGY.duration > 1.5);
    
    if (isempty(longDurations))
        fprintf('No long durations detected!\n');
    else
        fprintf('Long durations detected. Here are fix times:\n');
        disp(longDurations);
    end
    
    fprintf('Checking for NATs in CC before proceeding...\n');

    natIndex = [];

    q = 1;
    for n = 1:height(CC_ENERGY)
       if ( isnat(CC_ENERGY.Time(n)) )
            fprintf('Found a NAT in CC row %d. Adding to natIndex...\n', ...
                n);
            natIndex(q) = n;
            q = q + 1;
       end
    end

    if (~exist('natIndex','var'))
        natIndex = 0;
    end

    if (isempty(natIndex))

        fprintf('\tNo NATs detected in thisCC. Proceeding!\n')

    else

        numNats = numel(natIndex);

        if ( numNats > 1 )

            fprintf('NATs are consecutive. Dropping last %d rows.\n', ...
                numel(natIndex) );
            newEnd = natIndex(1) - 1;
            CC_ENERGY(natIndex(1):natIndex(end),:) = [];

        elseif (numNats == 1)

            fprintf('NAT last entry is NAT. Fixing...\n');
            CC_ENERGY(natIndex(1):natIndex(end-1),:) = [];

        elseif (numNats == 0)

            fprintf('\tNo NATs in thisCC! Proceeding...\n');

        else
            fprintf('NATs are in nonconsecutive order. Handle manually.\n');
            return;
        end

    end

    Time = TempTime;
    
    switch(epochChoice)
        case 1
            E1_CC_VEL_KE = CC_VEL_KE;
            E1_CC_ENERGY = CC_ENERGY;
            fprintf('\nAppending E1_CC_VEL_KE & E1_CC_ENERGY timetable to RAW file...\n');
            save(fullRawFileName, 'E1_CC_VEL_KE', 'E1_CC_ENERGY', '-append');  
            fprintf('E1_CC_VEL_KE & E1_CC_ENERGY added to RAW.\n');            
        case 2
            E2_CC_VEL_KE = CC_VEL_KE;
            E2_CC_ENERGY = CC_ENERGY;
            fprintf('\nAppending E2_CC_VEL_KE & E2_CC_ENERGY timetable to RAW file...\n');
            save(fullRawFileName, 'E2_CC_VEL_KE', 'E2_CC_ENERGY', '-append');  
            fprintf('E2_CC_VEL_KE & E2_CC_ENERGY added to RAW.\n');            
        case 3
            E3_CC_VEL_KE = CC_VEL_KE;
            E3_CC_ENERGY = CC_ENERGY;
            fprintf('\nAppending E3_CC_VEL_KE & E3_CC_ENERGY timetable to RAW file...\n');
            save(fullRawFileName, 'E3_CC_VEL_KE', 'E3_CC_ENERGY', '-append');  
            fprintf('E3_CC_VEL_KE & E3_CC_ENERGY added to RAW.\n');            
        case 4
            E4_CC_VEL_KE = CC_VEL_KE;
            E4_CC_ENERGY = CC_ENERGY;
            fprintf('\nAppending E4_CC_VEL_KE & E4_CC_ENERGY timetable to RAW file...\n');
            save(fullRawFileName, 'E4_CC_VEL_KE', 'E4_CC_ENERGY', '-append');  
            fprintf('E4_CC_VEL_KE & E4_CC_ENERGY added to RAW.\n');            
        case 5
            E5_CC_VEL_KE = CC_VEL_KE;
            E5_CC_ENERGY = CC_ENERGY;
            fprintf('\nAppending E5_CC_VEL_KE & E5_CC_ENERGY timetable to RAW file...\n');
            save(fullRawFileName, 'E5_CC_VEL_KE', 'E5_CC_ENERGY', '-append');  
            fprintf('E5_CC_VEL_KE & E5_CC_ENERGY added to RAW.\n');
    end
    

    % make a 3x3 plot of Axyz Vxyz and Gxyz mean SCG and GCG waveforms

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
    ylabel('rad · s^-^1');
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
    ylabel('rad · s^-^1');
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
    ylabel('rad · s^-^1');
    title('Cardiac cycle rotational velocity');
    grid;
    legend([ pMeanGz pCCGxyz(1,1) ], {'Gz','all Gz cycles'} );
    hold off;    

    linkaxes([s3a s3b s3c s3d s3e s3f s3g s3h s3i],'x');

end



%% Do some work on RSA intervals if present

doRsa = lower(input('Do RSA analysis (y/n) [n]: ','s'));
switch(doRsa)
    case 'y'
        doRsa = true;
    otherwise
        doRsa = false;
end

% if (doRsa && epochChoice ~= 3)
if (doRsa)

    figTimePlot = figure('Color', white);

    kFs = 100;

    clear ventIndex;
    clear ventEvents;

    ventEvents = (string(thisCC.cueType) == 'v');

    q = 1;
    for n = 1:height(ventEvents)
       if ( ventEvents(n) == 1 )
            fprintf('Found a ventilation event in row %d. Adding to ventIndex...\n', ...
                n);
            ventIndex(q) = n;
            q = q + 1;
       end
    end

    if (epochChoice == 3)
        numSeconds = seconds( thisCC.Time(2) - thisCC.Time(1) );
        time_s = 1:0.01:round(numSeconds);
    else
        time_s = 1:(1/kFs):max(diff(ventIndex));
    end

    fprintf('\n');

    mean_fH = mean(thisEpoch.hrTimeConsensus);

    switch(TRIAL.SUBJECT.id)
        case '9FL3'                 % Lono
            subjectMass   = 249.5455;
            subjectLength = 2.7432;
            subjectGirth  = 1.4542;
            
        case '6JK5'                 % Kolohe
            subjectMass = 210.9091;
            subjectLength = 2.5908;
            subjectGirth  = 1.4478;            
            
        case '9ON6'                 % Noa
            subjectMass = 203.5455;
            subjectLength = 2.4765;
            subjectGirth  = 1.3716;
            
        case '63HF'                         % Hoku
            subjectMass   = 184.4091;
            subjectLength = 2.54;           % 254 cm
            subjectGirth  = 1.3081;         % 130.81 cm
            
        case '01L5'                 % Liho
            subjectMass = 164.3636;
            subjectLength = 2.4003;
            subjectGirth  = 1.2446;            
            
        case '83H1'                 % Hua
            subjectMass = 152.2727;
            subjectLength = 2.413;
            subjectGirth  = 1.2573;            
            
            
        otherwise
            mass = input('Animal ID not found. Enter mass manually: ','s');
    end


    numVents = numel(ventIndex) - 1;
    animalID = TRIAL.SUBJECT.id;

    clear thisEpochData; 
    thisEpochData = struct;

    for vi = 1:numel(ventIndex)-1

        startAt = thisCC.Time(ventIndex(vi));
        endAt   = thisCC.Time(ventIndex(vi+1));

        thisDuration = duration( endAt - startAt);
        ibiDuration = seconds(thisDuration);

        v = timerange(startAt, endAt);

        ifH = thisEpoch.consensusHr(v) * 60;

        len_ifH = length(ifH);
        mean_ifH = mean(ifH);
        max_ifH = max(ifH);
        min_ifH = min(ifH);
        rsaper_thisRsa = ((max_ifH - min_ifH) / mean_ifH) * 100;
        fR_thisRsa = 60 / ibiDuration;

        % build the data table for this IBI cycle

        thisEpochData.start(vi) = startAt;
        thisEpochData.end(vi) = endAt;
        thisEpochData.samples(vi) = len_ifH;
        thisEpochData.mean_ifH(vi) = mean_ifH;
        thisEpochData.min_ifH(vi) = min_ifH;
        thisEpochData.max_ifH(vi) = max_ifH;
        thisEpochData.ibi(vi) = ibiDuration;
        thisEpochData.rsaper(vi) = rsaper_thisRsa;
        thisEpochData.fR(vi) = fR_thisRsa;
        thisEpochData.ID(vi) = {animalID};
        thisEpochData.mass(vi) = subjectMass;
        thisEpochData.length(vi) = subjectLength;
        thisEpochData.girth(vi) = subjectGirth;
        thisEpochData.epoch(vi) = epochChoice;

    %     figure(figOverPlot);
    %     hold on;
    %     plot(time_s(1:lenThisRsa), thisRsa);
    %     hold off;

        figure(figTimePlot);
        hold on;
        plot(thisEpoch.Time(v), ifH, 'LineWidth', 2);
        xline(endAt, '--', 'Color', red, 'LineWidth', 1.5);
        hold off;

    end         % end for vi = ...

    hold on;
    xlabel('Time, seconds');
    ylabel('iHR');
    grid;
    hold off;

    % move Time to placeholder
    TempTime = Time;

    % move thisEpochData.start to Time
    Time = thisEpochData.start;

    switch(epochChoice)
        case 1

            E1_VENT_DATA = timetable(Time', thisEpochData.end', ...
                thisEpochData.samples', thisEpochData.ID', ...
                thisEpochData.mass', thisEpochData.epoch', ...
                thisEpochData.ibi', thisEpochData.fR', ...
                thisEpochData.mean_ifH', thisEpochData.min_ifH', ...
                thisEpochData.max_ifH', thisEpochData.rsaper', ...
                'VariableNames',...
                {'end','numSamples','id','mass','epoch','ibi','fR', ...
                'mean_ifH','min_ifH','max_ifH','rsaper'  });

                fprintf('\nAppending E1_VENT_DATA timetable to RAW file...\n');
                save(fullRawFileName, 'E1_VENT_DATA', '-append');  
                fprintf('E1_VENT_DATA added to RAW.\n');


        case 3
            E3_VENT_DATA = timetable(Time', thisEpochData.end', ...
                thisEpochData.samples', thisEpochData.ID', ...
                thisEpochData.mass', thisEpochData.epoch', ...
                thisEpochData.ibi', thisEpochData.fR', ...
                thisEpochData.mean_ifH', thisEpochData.min_ifH', ...
                thisEpochData.max_ifH', thisEpochData.rsaper', ...
                'VariableNames',...
                {'end','numSamples','id','mass','epoch','ibi','fR', ...
                'mean_ifH','min_ifH','max_ifH','rsaper'  });        

                fprintf('\nAppending E3_VENT_DATA timetable to RAW file...\n');
                save(fullRawFileName, 'E3_VENT_DATA', '-append');  
                fprintf('E3_VENT_DATA added to RAW.\n');        

        case 5
            E5_VENT_DATA = timetable(Time', thisEpochData.end', ...
                thisEpochData.samples', thisEpochData.ID', ...
                thisEpochData.mass', thisEpochData.epoch', ...
                thisEpochData.ibi', thisEpochData.fR', ...
                thisEpochData.mean_ifH', thisEpochData.min_ifH', ...
                thisEpochData.max_ifH', thisEpochData.rsaper', ...
                'VariableNames',...
                {'end','numSamples','id','mass','epoch','ibi','fR', ...
                'mean_ifH','min_ifH','max_ifH','rsaper'  });        

                fprintf('\nAppending E5_VENT_DATA timetable to RAW file...\n');
                save(fullRawFileName, 'E5_VENT_DATA', '-append');  
                fprintf('E5_VENT_DATA added to RAW.\n');        

    end

    % move TempTime back into Time now that timetable is created
    Time = TempTime;

else
    
    fprintf('This is an apneic epoch, so skipping E3 vent work!\n');
    
end

