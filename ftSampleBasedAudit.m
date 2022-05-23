%% ftSampleBasedAudit.m

%   written by Dave Haas on 2-3 July 2021
%
%   This is a data-resize tool that uses samples on the x-axis (instead of
%   time), allowing the user to trim off pre-tagOn and post-tagOff samples
%   that would intefere with analysis. This tools is limited in scope to
%   thsi purpose. See ftTimeBasedAudit.m for annotating time series data
%   using a time (local) or time (in seconds), useful for aligning cues
%   entered in ftDataEntrySheets.m.
%   
%   

clc;
%clear;

%% Specify a tag for read

tag = 'tt21_134a';



%% Create a time slice for displaying the first 30 seconds of data

% create a time slice of the first 30 seconds of the record
timeSlice  = [ time(1)  (time(1) + duration([0 0 30])) ];

% create a recordSize in seconds
recordTime_s = time(end) - time(1);

cueStart = find(time == timeSlice(1));
cueEnd = find(time == timeSlice(2));

f1 = figure; 
plot(time_s(cueStart:cueEnd), hrBw_rLed2_25Hz(cueStart:cueEnd), 'b-'); 
hold on; 
plot(time_s(cueStart:cueEnd), hrBw_rLed3_25Hz(cueStart:cueEnd), 'r-');

grid;

% if old leftX and rightX are sitting around in memory, nuke 'em...
if ( exist('leftX','var') && exist('rightX', 'var') )
    clear leftX rightX;
end

trimComplete = false;
updatePlot = false;

while (~trimComplete)
    
    [x,y,button] = ginput(1);
    
    switch(char(button))
        case 'l'
            leftX = x;
            leftY = y;
            xline(leftX,'g','LineWidth',4);
        case 'r'
            rightX = x;
            rightY = y;
            xline(rightX,'r','LineWidth',4);
        case 'q'
            fprintf('Caught a `q` - quitting!\n');
            return;
        case 'w'
            windowSizeStr = input('Enter new window size in seconds: ','s');
            newWindowSize = str2double(windowSizeStr);
        case 'e'
            newCues = [ time(end) - duration([ 0 0 60]) time(end) ];
            cueStart = find(time == newCues(1));
            cueEnd = find(time == newCues(2));
            updatePlot = true;
        otherwise
            fprintf('use [l] [r] or [q]\n');
    end

    if ( exist('leftX','var') && exist('rightX', 'var') )
        fprintf('Left and right trim coordinates captured. Continuing...\n');
        trimComplete = true;
    end
    
    if (updatePlot)
        hold off;
        figure(f1);
        plot(time_s(cueStart:cueEnd), hrBw_rLed2_25Hz(cueStart:cueEnd), 'b-'); 
        hold on; 
        plot(time_s(cueStart:cueEnd), hrBw_rLed3_25Hz(cueStart:cueEnd), 'r-');
        grid;
        updatePlot = false;
    end
end


hold off;

inFmt = 'ss.SSS';
outFmt = 'hh:mm:ss.SSS';
leftD = duration([0 0 leftX], 'InputFormat', inFmt, 'Format', outFmt);
rightD = duration([0 0 rightX],'InputFormat',inFmt,'Format',outFmt);

