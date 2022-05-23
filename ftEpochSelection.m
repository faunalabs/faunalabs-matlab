%%  ftEpochSelection.m

%   Written by Dave Haas between 9 July and 4 August 2021

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
black = [0 0 0];                    % '#000000

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
fullPrhFileName = sprintf('%s/%sprh.mat', TAG_PATHS.PRH, ...
    tagStr);
if ( exist(fullRawFileName, 'file') )
    fprintf('Raw file located for this trial... loading it!\n');
    load(fullRawFileName);
    fprintf('PRH file located for this trial... loading it!\n');
    load(fullPrhFileName);
    
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

%% Do some regional selections with wsst and wsstridge finding

screenSize = get(0,'screensize');

figAnalysisRegion = figure;
figSize = figAnalysisRegion.Position;
figSize(3) = screenSize(3);
figSize(4) = screenSize(4) * 0.7;
figAnalysisRegion.Position = figSize;  

p1a = subplot(311);
plot(OPTICS.Time, OPTICS.led2, 'Color', blue);
hold on;
plot(OPTICS.Time, OPTICS.led3, 'Color', red );
yPos = p1a.YLim(2) - ((p1a.YLim(2) - p1a.YLim(1)) / 10);
ftPlotCues(CUE.Time, CUE.type, yPos);
hold off;
xlabel('Time, local');
ylabel('Intensity (A.U.)');
grid;
legend('1050 nm','1200 nm');

p1b = subplot(312);
plot(KINEMATICS.Time, KINEMATICS.odba, 'Color', green);
hold on;
yPos = p1b.YLim(2) - ((p1b.YLim(2) - p1b.YLim(1)) / 10);
ftPlotCues(CUE.Time, CUE.type, yPos);
hold off;
xlabel('Time, local');
ylabel('m/s^2');
title('ODBA');
grid;

p1c = subplot(313);
plot(PRESSURE.Time, PRESSURE.depth, 'Color', blue');
hold on;
yPos = p1c.YLim(2) - ((p1c.YLim(2) - p1c.YLim(1)) / 10);
ftPlotCues(CUE.Time, CUE.type, yPos);
hold off;
xlabel('Time, local');
ylabel('Depth, meters');
title('Tag Depth');
set(gca,'YDir','reverse');
grid;

linkaxes([p1a p1b p1c],'x');

if (exist('EPOCH','var'))
    
    ohTxt = 'It looks like EPOCH already exists. Redo EPOCH selection? (y/n): ';
    decideStr = lower(input(ohTxt,'s'));
    switch(decideStr)
        case 'y'
            fprintf('Clearing EPOCH and starting over...\n');
            clear EPOCH;
            epochsDefined = false;
        otherwise
            fprintf('OK, leaving EPOCH in place and bailing!\n');
            fprintf('Closing on figures on the way out...\n');
            ftCloseAllFigs;
            return;
    end
    
else
    
    epochsDefined = false;

end

while(~epochsDefined)
    
    fprintf('Click to set pre-apnea start (CFB).\n');
    [x1, ~, k1] = ginput(1);

    fprintf('Click to set pre-apnea end (CFB) / transition 1 (T1) start.\n');
    [x2, ~, k2] = ginput(1);

    fprintf('Click to set transition 1 (T1) end / apnea start (BH).\n');
    [x3, ~, k3] = ginput(1);
    
    fprintf('Click to set apnea end / transition 2 (T2) start.\n');
    [x4, ~, k4] = ginput(1);
    
	fprintf('Click to set transition 2 (T2) end / post-apnea recovery (REC).\n');
    [x5, ~, k5] = ginput(1);
    
	fprintf('Click to set post-apnea recovery (REC) end.\n');
    [x6, ~, k6] = ginput(1);    
    
    xAxis = gca;
    
    if (k1 == 1)
        t1 = num2ruler(x1, xAxis.XAxis);
    end
    
    if (k2 == 1)
        t2 = num2ruler(x2, xAxis.XAxis);
    end
    
    if (k3 == 1)
        t3 = num2ruler(x3, xAxis.XAxis);
    end
    
    if (k4 == 1)
        t4 = num2ruler(x4, xAxis.XAxis);
    end
    
    if (k5 == 1)
        t5 = num2ruler(x5, xAxis.XAxis);
    end
    
    if (k6 == 1)
        t6 = num2ruler(x6, xAxis.XAxis);
    end

    preApneaRange       = timerange(t1, t2);
    tOneRange           = timerange(t2, t3);
    apneaRange          = timerange(t3, t4);
    tTwoRange           = timerange(t4, t5);
    postApneaRange      = timerange(t5, t6);
    
    EPOCH = struct;
    
    EPOCH.preApnea  = preApneaRange;
    EPOCH.t1        = tOneRange;
    EPOCH.apnea     = apneaRange;
    EPOCH.t2        = tTwoRange;
    EPOCH.postApnea = postApneaRange;
    
    oT_preApnea  = OPTICS.Time(preApneaRange);
    oT_tOne      = OPTICS.Time(tOneRange);
    oT_apnea     = OPTICS.Time(apneaRange);
    oT_tTwo      = OPTICS.Time(tTwoRange);
    oT_postApnea = OPTICS.Time(postApneaRange);
    
    kT_preApnea  = KINEMATICS.Time(preApneaRange);
    kT_tOne      = KINEMATICS.Time(tOneRange);
    kT_apnea     = KINEMATICS.Time(apneaRange);
    kT_tTwo      = KINEMATICS.Time(tTwoRange);
    kT_postApnea = KINEMATICS.Time(postApneaRange);
    
      
    figure(figAnalysisRegion);

    p1a = subplot(311);
    plot(OPTICS.Time, OPTICS.led2, 'Color', blue);
    hold on;
    plot(OPTICS.Time, OPTICS.led3, 'Color', red );
    xline(oT_preApnea(1), '--', 'Color', green);
    xline(oT_preApnea(end), '-.', 'Color', green);
    xline(oT_tOne(1), 'm--');
    xline(oT_tOne(end), 'm-.');
    xline(oT_apnea(1), 'r--');
    xline(oT_apnea(end), 'r-.');
    xline(oT_tTwo(1), 'm--');
    xline(oT_tTwo(end), 'm-.');
    xline(oT_postApnea(1), 'b--');
    xline(oT_postApnea(end), 'b-.');
    yPos = p1a.YLim(2) - ((p1a.YLim(2) - p1a.YLim(1)) / 10);
    ftPlotCues(CUE.Time, CUE.type, yPos);
    hold off;
    xlabel('Time, local');
    ylabel('Intensity (A.U.)');
    grid;
    legend('1050 nm','1200 nm');

    p1b = subplot(312);
    plot(KINEMATICS.Time, KINEMATICS.odba, 'Color', green);
    hold on;
    xline(kT_preApnea(1), '--', 'Color', green);
    xline(kT_preApnea(end), '-.', 'Color', green);
    xline(kT_tOne(1), 'm--');
    xline(kT_tOne(end), 'm-.');
    xline(kT_apnea(1), 'r--');
    xline(kT_apnea(end), 'r-.');
    xline(kT_tTwo(1), 'm--');
    xline(kT_tTwo(end), 'm-.');
    xline(kT_postApnea(1), 'b--');
    xline(kT_postApnea(end), 'b-.');
    yPos = p1b.YLim(2) - ((p1b.YLim(2) - p1b.YLim(1)) / 10);
    ftPlotCues(CUE.Time, CUE.type, yPos);
    hold off;
    xlabel('Time, local');
    ylabel('m/s^2');
    title('ODBA');
    grid;

    p1c = subplot(313);
    plot(PRESSURE.Time, PRESSURE.depth, 'Color', blue');
    hold on;
    xline(kT_preApnea(1), '--', 'Color', green);
    xline(kT_preApnea(end), '-.', 'Color', green);
    xline(kT_tOne(1), 'm--');
    xline(kT_tOne(end), 'm-.');
    xline(kT_apnea(1), 'r--');
    xline(kT_apnea(end), 'r-.');
    xline(kT_tTwo(1), 'm--');
    xline(kT_tTwo(end), 'm-.');
    xline(kT_postApnea(1), 'b--');
    xline(kT_postApnea(end), 'b-.');
    yPos = p1c.YLim(2) - ((p1c.YLim(2) - p1c.YLim(1)) / 10);
    ftPlotCues(CUE.Time, CUE.type, yPos);
    hold off;    
    xlabel('Time, local');
    ylabel('Depth, meters');
    title('Tag Depth');
    set(gca,'YDir','reverse');
    grid;

    linkaxes([p1a p1b p1c],'x');

    % ask user if this is good enough to select an end or do more zooming
    goodRegionTxt = 'Is this region good for your analysis? (y/n): ';
    thatsGood = lower(input(goodRegionTxt,'s'));        

    switch(thatsGood)
        case 'y'
            epochsDefined = true;
        otherwise
            epochsDefined = false;
    end
    
end


%% append EPOCH structure to RAW file

if (epochsDefined)
   
    fprintf('Appending EPOCH data structure to RAW file: %s\n', fullRawFileName );

    if (exist(fullRawFileName, 'file'))

        fprintf('RAW file exists... ');

        prompt = 'OK to append new EPOCH data? y/n [n]: ';
        str = input(prompt,'s');
        if isempty(str)
            str = 'n';
        end

        if (strcmp(str, 'y'))
            
            fprintf('\nAppending new data to RAW file...\n');
            save(fullRawFileName, 'EPOCH', '-append');   
            save(fullPrhFileName, 'EPOCH', '-append');
            
            fprintf('EPOCHs appear to have been created and added to RAW.\n');

            fprintf('*** END ftEpochSelection ***\n');
            
        end 

    else

        fprintf('Raw file does not exist... this should never happen... Bail!\n');
        return;

    end
        
end

%% 
