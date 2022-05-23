%%  ftHeartRate_Optics.m

%   Written by Dave Haas between 9 July and 14 August 2021
%
%   Last updated by Dave Haas on 7 August 2021, to make wsstridges from
%   findpeaks data.

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
black = [0 0 0];                    % '#000000'
white = [1 1 1];                    % '#FFFFFF'

% set some default LED* colors...

led1Color = black;
led2Color = blue;
led3Color = red;
led4Color = goldenrod;

% ... and some default kinematic colors...

kxColor = blue;
kyColor = red;
kzColor = goldenrod;

%% Make introductions...

fprintf('Welcome to the FaunaTag optical heart rate analysis app.\n');


%% select a tag for analysis

tagStr = input('Enter a tag for analysis (e.g.: tt21_141a): ','s');

fullRawFileName = sprintf('%s/%sraw.mat', TAG_PATHS.RAW, ...
    tagStr);
if ( exist(fullRawFileName, 'file') )
    fprintf('Raw file located for this trial... loading it!\n');
    load(fullRawFileName);
    
    % confirm presence of OPTICS variable
    if ( exist('EPOCH','var') )
        
%         oT_s = OPTICS.time_s;
%         oT = OPTICS.Time;
        
        if (exist('OPTICS.fs','var'))
            oFs = OPTICS.Properties.SampleRate;
        else
            oFs = 250;   % assume 250 Hz sampling on FaunaTags
        end
        
        % define an optical Nyquist frequency
        oNyquistFs          = oFs / 2;

%         kT_s = KINEMATICS.time_s;
%         kT = KINEMATICS.Time;

        if (exist('KINEMATICS.fs','var'))
            kFs = KINEMATICS.Properties.SampleRate;
        else
            kFs = 100;   % assume 100 Hz sampling on FaunaTags
        end
        
        % define a kinematic Nyquist frequency
        kNyquistFs          = kFs / 2;
        
        % define a tag trial name
        tag = CUE.Properties.CustomProperties.tag;
        
        fprintf('Key variables are present and loaded... proceeding!\n');

    else
        
        fprintf('Key variables like EPOCH are not present. Check RAW file or run ftEpochSelection.\n');
        return;
        
    end
    
else
    
    fprintf('Raw file is not present. Run ftPostProcessing.m to create it.\n');
    return;

end


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
                fprintf('E1_KINEMATICS found and copied as E for comparisons.\n');
                E = E1_KINEMATICS;
            end
            
            if (exist('E1_OPTICS','var'))
                confText = 'Optical analysis results already exist. Redo? (y/n): ';
                confStr = lower(input(confText,'s'));
                switch(confStr)
                    case 'y'
                        analysisRange = EPOCH.preApnea;
                        validEpoch = true;                
                    otherwise
                        fprintf('Skipping re-analysis. Choose a new epoch!\n');
                        validEpoch = false;
                end                
            else
                analysisRange = EPOCH.preApnea;
                validEpoch = true;                
            end
            
        case 2      % TRANSITION 1
            
            if (exist('E2_KINEMATICS','var'))
                fprintf('E2_KINEMATICS found and copied as E for comparisons.\n');
                E = E2_KINEMATICS;
            end
            if (exist('E2_OPTICS','var'))
                confText = 'Optical analysis results already exist. Redo? (y/n): ';
                confStr = lower(input(confText,'s'));
                switch(confStr)
                    case 'y'
                        analysisRange = EPOCH.t1;
                        validEpoch = true;                
                    otherwise
                        fprintf('Skipping re-analysis. Choose a new epoch!\n');
                        validEpoch = false;
                end                
            else
                analysisRange = EPOCH.t1;
                validEpoch = true;                
            end
            
            
        case 3      % APNEA
            
            if (exist('E3_KINEMATICS','var'))
                fprintf('E3_KINEMATICS found and copied as E for comparisons.\n');
                E = E3_KINEMATICS;
            end
            if (exist('E3_OPTICS','var'))
                confText = 'Optical analysis results already exist. Redo? (y/n): ';
                confStr = lower(input(confText,'s'));
                switch(confStr)
                    case 'y'
                        analysisRange = EPOCH.apnea;
                        validEpoch = true;                
                    otherwise
                        fprintf('Skipping re-analysis. Choose a new epoch!\n');
                        validEpoch = false;
                end                
            else
                analysisRange = EPOCH.apnea;
                validEpoch = true;                
            end
            
        case 4      % TRANSITION 2
            
            if (exist('E4_KINEMATICS','var'))
                fprintf('E4_KINEMATICS found and copied as E for comparisons.\n');
                E = E4_KINEMATICS;
            end
            if (exist('E4_OPTICS','var'))
                confText = 'Analysis results already exist. Redo? (y/n): ';
                confStr = lower(input(confText,'s'));
                switch(confStr)
                    case 'y'
                        analysisRange = EPOCH.t2;
                        validEpoch = true;                
                    otherwise
                        fprintf('Skipping re-analysis. Choose a new epoch!\n');
                        validEpoch = false;
                end                
            else
                analysisRange = EPOCH.t2;
                validEpoch = true;                
            end
            
        case 5      % POST-APNEA
            
            if (exist('E5_KINEMATICS','var'))
                fprintf('E5_KINEMATICS found and copied as E for comparisons.\n');
                E = E5_KINEMATICS;
            end
            if (exist('E5_OPTICS','var'))
                confText = 'Optical analysis results already exist. Redo? (y/n): ';
                confStr = lower(input(confText,'s'));
                switch(confStr)
                    case 'y'
                        analysisRange = EPOCH.postApnea;
                        validEpoch = true;                
                    otherwise
                        fprintf('Skipping re-analysis. Choose a new epoch!\n');
                        validEpoch = false;
                end                
            else
                analysisRange = EPOCH.postApnea;
                validEpoch = true;                
            end
            
        otherwise
            fprintf('That was not a valid epoch choice. Try again.\n');
    end    
end


%% Show data from the selected epoch

screenSize = get(0,'screensize');

figAnalysisRegion = figure('Color', white);
figSize = figAnalysisRegion.Position;
figSize(3) = screenSize(3) * 0.9;
figSize(4) = screenSize(4) * 0.62;
figAnalysisRegion.Position = figSize;  

analysisRegionDefined = false;

while(~analysisRegionDefined)

    analysisTime = OPTICS.Time(analysisRange);
      
    figure(figAnalysisRegion);

    p1a = subplot(311);
    plot(OPTICS.Time(analysisRange), OPTICS.led2(analysisRange), 'Color', blue);
    hold on;
    plot(OPTICS.Time(analysisRange), OPTICS.led3(analysisRange), 'Color', red );
    hold off;
    xlabel('Time, local');
    ylabel('Intensity (A.U.)');
    grid;
    legend('1050 nm','1200 nm');

    p1b = subplot(312);
    plot(KINEMATICS.Time(analysisRange), KINEMATICS.odba(analysisRange), 'Color', green);
    xlabel('Time, local');
    ylabel('m/s^2');
    title('ODBA');
    grid;

    p1c = subplot(313);
    plot(PRESSURE.Time(analysisRange), PRESSURE.depth(analysisRange), 'Color', blue');
    xlabel('Time, local');
    ylabel('Depth, meters');
    title('Tag Depth');
    set(gca,'YDir','reverse');
    grid;

    linkaxes([p1a p1b p1c],'x');

    % ask user if this is good enough to select an end or do more zooming
    goodRegionTxt = 'Is this region good for your analysis? (y/n) [y]: ';
    thatsGood = lower(input(goodRegionTxt,'s'));        

    switch(thatsGood)
        case 'n'
            analysisRegionDefined = false;
            
        otherwise
            analysisRegionDefined = true;
            
            oT = OPTICS.Time(analysisRange);
            oT_s = OPTICS.time_s(analysisRange);

            oLed1 = OPTICS.led1(analysisRange);
            oLed2 = OPTICS.led2(analysisRange);
            oLed3 = OPTICS.led3(analysisRange);
            oLed4 = OPTICS.led4(analysisRange);

            kT = KINEMATICS.Time(analysisRange);
            kT_s = KINEMATICS.time_s(analysisRange);

            kAx = KINEMATICS.ax(analysisRange);
            kAy = KINEMATICS.ay(analysisRange);
            kAz = KINEMATICS.az(analysisRange);

            kGx = KINEMATICS.gx(analysisRange);
            kGy = KINEMATICS.gy(analysisRange);
            kGz = KINEMATICS.gz(analysisRange);

            kPitch = KINEMATICS.pitch(analysisRange);
            kRoll = KINEMATICS.pitch(analysisRange);
            kHeading = KINEMATICS.heading(analysisRange);

            kOdba = KINEMATICS.odba(analysisRange);
            kOdav = KINEMATICS.odav(analysisRange);

            pDepth = PRESSURE.depth(analysisRange);
            pTemperature = PRESSURE.temperature(analysisRange);

            timeSubselect = true;
            
        
    end
    

end


%% Reflectance version of RAW optical data --> only use rLed* after this!

% invert for reflectance typical waveforms

led1 = oLed1;
led2 = oLed2;
led3 = oLed3;
led4 = oLed4;

iLed1 = oLed1 * -1;
iLed2 = oLed2 * -1;
iLed3 = oLed3 * -1;
iLed4 = oLed4 * -1;

% get each LED full scale range
rangeLed1 = range(iLed1);
rangeLed2 = range(iLed2);
rangeLed3 = range(iLed3);
rangeLed4 = range(iLed4);

% rescale each LED from 1 to range(led*)
rLed1 = rescale(iLed1, 1, rangeLed1);
rLed2 = rescale(iLed2, 1, rangeLed2);
rLed3 = rescale(iLed3, 1, rangeLed3);
rLed4 = rescale(iLed4, 1, rangeLed4);

% Apply some denoise tests...

oLength = numel(oLed2);
kLength = numel(kAx);
okRatio = round( (oLength / kLength), 2);

if ( okRatio ~= 2.5000 )
   fprintf('Time series mismatch. oLength / kLength = %1.2f.\n', ...
       okRatio);
   fprintf('Check figure for off-by-one-sample and approve if all is okay.\n');
   
   figure('Position',[100 100 500 500], 'Color', white);
   plot(oT_s, normalize(oLed2), kT_s, kAx);
   xlabel('Time, seconds');
   ylabel('oLed2 & kAx');
   grid;
   
   okToProceed = lower(input('Okay to proceed (y/n) [y]: ','s'));
   switch(okToProceed)
       case 'n'
           fprintf('Halting. Fix whatever is wrong and try again later.\n');
           return;
       otherwise
           fprintf('Okay! Proceeding with analysis...\n');
   end
   
   
else
    fprintf('Time series matches. O:K ratio = %1.2f.\n', okRatio);
end


%% New bandpass worked great for tt21_142d E3

switch(epochChoice)
    case 3                  % apnea
        bandPassFiltOrder  = 3;      % 3rd order butterworth filter
        oLowCut            = 1.0;    % ~60 beats per minute
        oHighCut           = 1.5;    % ~90 beats per minute
        opticalTestFs      = 8;      % 8 produces good pulses on optics apnea        
    otherwise               % pre- or post-apnea (E1, E2, E4, or E5)
        bandPassFiltOrder  = 3;      % 6th order butterworth filter
        oLowCut            = 0.5;    % ~30 beats per minute
        oHighCut           = 3.5;   % ~210 beats per minute <- 3.5 test 12.5
        %oHighCut           = 12.5;   % ~210 beats per minute <- 3.5 test 12.5
        opticalTestFs      = oHighCut * 2;     % 8 produces good pulses on optics apnea        
        %opticalTestFs      = 50;               % 50 is a bit noisier
        % opticalTestFs       = 125;
end

[Aobp, Bobp, Cobp, Dobp] = butter(bandPassFiltOrder, ...
    [oLowCut oHighCut] / (opticalTestFs) );

if (mod(bandPassFiltOrder,2) == 0)
    designFiltOrder = bandPassFiltOrder;
else
    designFiltOrder = bandPassFiltOrder + 1;
end

bandPassOpticalFilter = designfilt('bandpassiir','FilterOrder', designFiltOrder, ...
    'HalfPowerFrequency1',oLowCut,'HalfPowerFrequency2', oHighCut, ...
	'SampleRate', opticalTestFs);

[sosBpo,gBpo] = ss2sos(Aobp,Bobp,Cobp,Dobp);
fvtO = fvtool(sosBpo, bandPassOpticalFilter, 'Fs', opticalTestFs);
legend(fvtO,'butter','designfilt');

OPTICS_BANDPASS_FILTER = struct;
OPTICS_BANDPASS_FILTER.filterOrder = bandPassFiltOrder;
OPTICS_BANDPASS_FILTER.lowCut = oLowCut;
OPTICS_BANDPASS_FILTER.highCut = oHighCut;
OPTICS_BANDPASS_FILTER.sampleRate = opticalTestFs;
OPTICS_BANDPASS_FILTER.A = Aobp;
OPTICS_BANDPASS_FILTER.B = Bobp;
OPTICS_BANDPASS_FILTER.C = Cobp;
OPTICS_BANDPASS_FILTER.D = Dobp;
OPTICS_BANDPASS_FILTER.SOS = sosBpo;
OPTICS_BANDPASS_FILTER.G = gBpo;

fprintf('Band-pass butterworth filter for optical cardio energies:\n');
fprintf('\tFilter order: %d\n', bandPassFiltOrder);
fprintf('\tLow cut frequency: %2.2f Hz\n', oLowCut);
fprintf('\tHigh cut frequency: %2.2f Hz\n', oHighCut);
fprintf('\tOptical BPF sample rate (2Hz = normalized frequencies): %d Hz\n', ...
    opticalTestFs);



%% Use this filter to make good optical heart rate signals

filtLed1 = filtfilt( sosBpo, gBpo, oLed1);
filtLed2 = filtfilt( sosBpo, gBpo, oLed2);
filtLed3 = filtfilt( sosBpo, gBpo, oLed3);
filtLed4 = filtfilt( sosBpo, gBpo, oLed4);

meanFiltLed1 = mean(filtLed1);
meanFiltLed2 = mean(filtLed2);
meanFiltLed3 = mean(filtLed3);
meanFiltLed4 = mean(filtLed4);

stdFiltLed1 = std(filtLed1);
stdFiltLed2 = std(filtLed2);
stdFiltLed3 = std(filtLed3);
stdFiltLed4 = std(filtLed4);

r95thFiltLed1 = meanFiltLed1 + stdFiltLed1 * 1.96;
l95thFiltLed1 = meanFiltLed1 - stdFiltLed1 * 1.96;

r95thFiltLed2 = meanFiltLed2 + stdFiltLed2 * 1.96;
l95thFiltLed2 = meanFiltLed2 - stdFiltLed2 * 1.96;

r95thFiltLed3 = meanFiltLed3 + stdFiltLed3 * 1.96;
l95thFiltLed3 = meanFiltLed3 - stdFiltLed3 * 1.96;

r95thFiltLed4 = meanFiltLed4 + stdFiltLed4 * 1.96;
l95thFiltLed4 = meanFiltLed4 - stdFiltLed4 * 1.96;

fprintf('Displaying all band-passed optical channels, and displaying.\n');
fprintf('histograms and summary status for band-passed optical channels.\n');


figOptics = figure('Position',[100 100 1600 1200], 'Color', white);
figOptics.Name = 'Band-pass filtered optical channels';

s1a = subplot(421);
plot(oT, filtLed1, 'Color', black);
% hold on; 
% plot(oT, dFiltLed1, 'Color', [0.5 0.5 0.5]);
% hold off;
xlabel('Time, local');
ylabel('ambient');
grid;

s1b = subplot(423);
plot(oT, filtLed2, 'Color', blue);
% hold on; 
% plot(oT, dFiltLed2, 'Color', [0.0 0.0 0.5]);
% hold off;
xlabel('Time, local');
ylabel('1050 nm');
grid;

s1c = subplot(425);
plot(oT, filtLed3, 'Color', red);
% hold on; 
% plot(oT, dFiltLed3, 'Color', [0.5 0.0 0.0]);
% hold off;
xlabel('Time, local');
ylabel('1200 nm');
grid;

s1d = subplot(427);
plot(oT, filtLed4, 'Color', green);
% hold on; 
% plot(oT, dFiltLed4, 'Color', [0.0 0.5 0.0]);
% hold off;
xlabel('Time, local');
ylabel('amb2 (nonOp 950 nm)');
grid;

fHist1 = subplot(422);
histogram(filtLed1, 'FaceColor', [.7 .7 .7], 'EdgeColor', black);
xline(meanFiltLed1, 'Color', goldenrod, 'LineWidth', 1.2);
xline(r95thFiltLed1, '--', 'Color', goldenrod, 'LineWidth', 1.2);
xline(l95thFiltLed1, '--', 'Color', goldenrod, 'LineWidth', 1.2);
title('Channel 1: ambient');
grid;

fHist2 = subplot(424);
histogram(filtLed2, 'FaceColor', blue, 'EdgeColor', black);
xline(meanFiltLed2, 'Color', goldenrod, 'LineWidth', 1.2);
xline(r95thFiltLed2, '--', 'Color', goldenrod, 'LineWidth', 1.2);
xline(l95thFiltLed2, '--', 'Color', goldenrod, 'LineWidth', 1.2);
title('Channel 2: 1050 nm');
grid;

fHist3 = subplot(426);
histogram(filtLed3, 'FaceColor', red, 'EdgeColor', black);
xline(meanFiltLed3, 'Color', goldenrod, 'LineWidth', 1.2);
xline(r95thFiltLed3, '--', 'Color', goldenrod, 'LineWidth', 1.2);
xline(l95thFiltLed3, '--', 'Color', goldenrod, 'LineWidth', 1.2);
title('Channel 3: 1200 nm');
grid;

fHist4 = subplot(428);
histogram(filtLed4, 'FaceColor', green, 'EdgeColor', black);
xline(meanFiltLed4, 'Color', goldenrod, 'LineWidth', 1.2);
xline(r95thFiltLed4, '--', 'Color', goldenrod, 'LineWidth', 1.2);
xline(l95thFiltLed4, '--', 'Color', goldenrod, 'LineWidth', 1.2);
title('Channel 4: ambient, inoperative 950 nm)');
grid;

linkaxes([s1a s1b s1c s1d],'x');


%% compare how well HR pulses are preserved using filloutliers w/ spline

fprintf('Plotting outlier-replacement with splining for optical channels...\n');

[deOutFiltLed1, deOutIndexLed1, deOutLoThreshLed1, deOutHiThreshLed1 ] = ...
    filloutliers(filtLed1, 'spline');

[deOutFiltLed2, deOutIndexLed2, deOutLoThreshLed2, deOutHiThreshLed2 ] = ...
    filloutliers(filtLed2, 'spline');

[deOutFiltLed3, deOutIndexLed3, deOutLoThreshLed3, deOutHiThreshLed3 ] = ...
    filloutliers(filtLed3, 'spline');

[deOutFiltLed4, deOutIndexLed4, deOutLoThreshLed4, deOutHiThreshLed4 ] = ...
    filloutliers(filtLed4, 'spline');


figDeOut = figure('Position',[100 100 1200 1400 ], 'Color', white);
figDeOut.Name = 'filloutliers with splining for optical channels';

sOut1 = subplot(411);
plot(oT, filtLed1, '--', 'Color', goldenrod);
hold on;
plot(oT, deOutFiltLed1, 'Color', black);
hold off
xlabel('Time, local');
ylabel('AU');
title('Channel 1: ambient');
grid;
legend('pre-outlier signal','outlier-replaced signal');

sOut2 = subplot(412);
plot(oT, filtLed2, '--', 'Color', goldenrod);
hold on;
plot(oT, deOutFiltLed2, 'Color', blue);
hold off
xlabel('Time, local');
ylabel('AU');
title('Channel 2: 1050 nm');
grid;
legend('pre-outlier signal','outlier-replaced signal');

sOut3 = subplot(413);
plot(oT, filtLed3, '--', 'Color', goldenrod);
hold on;
plot(oT, deOutFiltLed3, 'Color', red);
hold off
xlabel('Time, local');
ylabel('AU');
title('Channel 3: 1200 nm');
grid;
legend('pre-outlier signal','outlier-replaced signal');

sOut4 = subplot(414);
plot(oT, filtLed4, '--', 'Color', goldenrod);
hold on;
plot(oT, deOutFiltLed4, 'Color', green);
hold off
xlabel('Time, local');
ylabel('AU');
title('Channel 4: ambient (inoperable 950 nm)');
grid;
legend('pre-outlier signal','outlier-replaced signal');

linkaxes([sOut1,sOut2,sOut3,sOut4],'x');  


%% replace original optical channel data with outlier-smoothed optical data?

replaceTest = lower(input('Use outlier-smoothed optical data (y/n) [y]: ','s'));

switch(replaceTest)
    case 'n'
        
        outliersSmoothed = false;
        fprintf('Leaving original band-passed optical data in place.\n');     
        
    otherwise
        
        outliersSmoothed = true;
        fprintf('Using outlier-smoothed band-passed optical data...\n');
        
        filtLed1 = deOutFiltLed1;
        filtLed2 = deOutFiltLed2;
        filtLed3 = deOutFiltLed3;
        filtLed4 = deOutFiltLed4;

        fprintf('Here are histograms after outlier-smoothing...\n');

end

%   to a look at histograms if outlier-smoothing is used

if (outliersSmoothed)

meanFiltLed1 = mean(filtLed1);
meanFiltLed2 = mean(filtLed2);
meanFiltLed3 = mean(filtLed3);
meanFiltLed4 = mean(filtLed4);

stdFiltLed1 = std(filtLed1);
stdFiltLed2 = std(filtLed2);
stdFiltLed3 = std(filtLed3);
stdFiltLed4 = std(filtLed4);

r95thFiltLed1 = meanFiltLed1 + stdFiltLed1 * 1.96;
l95thFiltLed1 = meanFiltLed1 - stdFiltLed1 * 1.96;

r95thFiltLed2 = meanFiltLed2 + stdFiltLed2 * 1.96;
l95thFiltLed2 = meanFiltLed2 - stdFiltLed2 * 1.96;

r95thFiltLed3 = meanFiltLed3 + stdFiltLed3 * 1.96;
l95thFiltLed3 = meanFiltLed3 - stdFiltLed3 * 1.96;

r95thFiltLed4 = meanFiltLed4 + stdFiltLed4 * 1.96;
l95thFiltLed4 = meanFiltLed4 - stdFiltLed4 * 1.96;

figOptics = figure('Position',[100 100 1600 1200], 'Color', white);
figOptics.Name = 'Band-pass filtered optical channels';

s1a = subplot(421);
plot(oT, filtLed1, 'Color', black);
% hold on; 
% plot(oT, dFiltLed1, 'Color', [0.5 0.5 0.5]);
% hold off;
xlabel('Time, local');
ylabel('ambient');
grid;

s1b = subplot(423);
plot(oT, filtLed2, 'Color', blue);
% hold on; 
% plot(oT, dFiltLed2, 'Color', [0.0 0.0 0.5]);
% hold off;
xlabel('Time, local');
ylabel('1050 nm');
grid;

s1c = subplot(425);
plot(oT, filtLed3, 'Color', red);
% hold on; 
% plot(oT, dFiltLed3, 'Color', [0.5 0.0 0.0]);
% hold off;
xlabel('Time, local');
ylabel('1200 nm');
grid;

s1d = subplot(427);
plot(oT, filtLed4, 'Color', green);
% hold on; 
% plot(oT, dFiltLed4, 'Color', [0.0 0.5 0.0]);
% hold off;
xlabel('Time, local');
ylabel('amb2 (nonOp 950 nm)');
grid;

fHist1 = subplot(422);
histogram(filtLed1, 'FaceColor', [.7 .7 .7], 'EdgeColor', black);
xline(meanFiltLed1, 'Color', goldenrod, 'LineWidth', 1.2);
xline(r95thFiltLed1, '--', 'Color', goldenrod, 'LineWidth', 1.2);
xline(l95thFiltLed1, '--', 'Color', goldenrod, 'LineWidth', 1.2);
title('Channel 1: ambient');
grid;

fHist2 = subplot(424);
histogram(filtLed2, 'FaceColor', blue, 'EdgeColor', black);
xline(meanFiltLed2, 'Color', goldenrod, 'LineWidth', 1.2);
xline(r95thFiltLed2, '--', 'Color', goldenrod, 'LineWidth', 1.2);
xline(l95thFiltLed2, '--', 'Color', goldenrod, 'LineWidth', 1.2);
title('Channel 2: 1050 nm');
grid;

fHist3 = subplot(426);
histogram(filtLed3, 'FaceColor', red, 'EdgeColor', black);
xline(meanFiltLed3, 'Color', goldenrod, 'LineWidth', 1.2);
xline(r95thFiltLed3, '--', 'Color', goldenrod, 'LineWidth', 1.2);
xline(l95thFiltLed3, '--', 'Color', goldenrod, 'LineWidth', 1.2);
title('Channel 3: 1200 nm');
grid;

fHist4 = subplot(428);
histogram(filtLed4, 'FaceColor', green, 'EdgeColor', black);
xline(meanFiltLed4, 'Color', goldenrod, 'LineWidth', 1.2);
xline(r95thFiltLed4, '--', 'Color', goldenrod, 'LineWidth', 1.2);
xline(l95thFiltLed4, '--', 'Color', goldenrod, 'LineWidth', 1.2);
title('Channel 4: ambient, inoperative 950 nm)');
grid;

linkaxes([s1a s1b s1c s1d],'x');

end

%%

dcoFiltOrder        = 3;                % FieldTrip suggested default
dcoPassbandFreqHz   = 0.01;             % DC-removal (cut below 0.01 Hz)

fprintf('DC-removal (high-pass) filter characteristics:\n');
fprintf('\tFilter order: %d\n', dcoFiltOrder);
fprintf('\tPassband frequency: %d Hz\n', dcoPassbandFreqHz);
fprintf('\tNyquist frequency (optics): %d\n', oNyquistFs);
fprintf('\tNyquist frequency (kinematics): %d\n', kNyquistFs);

% filter coefficients for optics high-pass (DC-removal)
[AoDCHP, BoDCHP, CoDCHP, DoDCHP] = ...
    butter(dcoFiltOrder, ( dcoPassbandFreqHz / (oNyquistFs) ), 'high');

[oSOS,oG] = ss2sos(AoDCHP, BoDCHP, CoDCHP, DoDCHP);

% design the high-pass filter. NOTE: using samplingRate = 2, which is the
% designfilt default, results in normalized frequencies

dcoOpticalFilter = designfilt('highpassiir','FilterOrder',dcoFiltOrder, ...
    'PassbandFrequency', dcoPassbandFreqHz, ... 
    'PassbandRipple', 0.2, ...
    'SampleRate', 2);

fvt = fvtool(oSOS, dcoOpticalFilter, 'Fs', oFs);
legend(fvt,'butter','designfilt');


%% Use Gaussian window filtering to smooth optical pulse data

% set a 100 ms windowSize to use for gaussian filtering
%optGaussianWindowSize = oFs / 10;      % = 0.040 s
optGaussianWindowSize = oFs / 25;       % = 0.100 s

gfBo = ones(1, optGaussianWindowSize) / optGaussianWindowSize;
gfA = [1];

gfFiltLed1 = filtfilt(gfBo, gfA, filtLed1);
gfFiltLed2 = filtfilt(gfBo, gfA, filtLed2);
gfFiltLed3 = filtfilt(gfBo, gfA, filtLed3);
gfFiltLed4 = filtfilt(gfBo, gfA, filtLed4);

figGaussOpt = figure('Position',[100 500 1500 600], 'Color', white);

pGf1 = subplot(411);
plot(oT, filtLed1, 'Color', goldenrod);
hold on;
plot(oT, gfFiltLed1, 'Color', black);
hold off;
grid;
legend('bpFilt','gf+bp');

pGf2 = subplot(412);
plot(oT, filtLed2, 'Color', goldenrod);
hold on;
plot(oT, gfFiltLed2, 'Color', blue);
hold off;
grid;
legend('bpFilt','gf+bp');

pGf3 = subplot(413);
plot(oT, filtLed3, 'Color', goldenrod);
hold on;
plot(oT, gfFiltLed3, 'Color', red);
hold off;
grid;
legend('bpFilt','gf+bp');

pGf4 = subplot(414);
plot(oT, filtLed4, 'Color', goldenrod);
hold on;
plot(oT, gfFiltLed4, 'Color', green);
hold off;
grid;
legend('bpFilt','gf+bp');

linkaxes([pGf1 pGf2 pGf3 pGf4],'x');

%% make root sum squared versions of these optical channels

%   NOTE: rss = root sum squared, i.e.: the Euclidian norm of a vector

rssFiltLed1 = sqrt(filtLed1.^2);
rssFiltLed2 = sqrt(filtLed2.^2);
rssFiltLed3 = sqrt(filtLed3.^2);
rssFiltLed4 = sqrt(filtLed4.^2);

rssActiveLeds = sqrt( filtLed2.^2 + filtLed3.^2);
rssAmbientLeds = sqrt( filtLed1.^2 + filtLed4.^2 );

rssGfFiltLed1 = sqrt(gfFiltLed1.^2);
rssGfFiltLed2 = sqrt(gfFiltLed2.^2);
rssGfFiltLed3 = sqrt(gfFiltLed3.^2);
rssGfFiltLed4 = sqrt(gfFiltLed4.^2);

rssGfActiveLeds = sqrt( gfFiltLed2.^2 + gfFiltLed3.^2);
rssGfAmbientLeds = sqrt( gfFiltLed1.^2 + gfFiltLed4.^2 );

figure;

s1a = subplot(321);
plot(oT_s, rssFiltLed2, 'Color', blue);
grid;
s1b = subplot(323);
plot(oT_s, rssFiltLed3, 'Color', red);
grid;
s1c = subplot(325);
plot(oT_s, rssActiveLeds, 'Color', purple);
grid;

s1d = subplot(322);
plot(oT_s, rssGfFiltLed2, 'Color', blue);
grid;
s1d = subplot(324);
plot(oT_s, rssGfFiltLed3, 'Color', red);
grid;
s1e = subplot(326);
plot(oT_s, rssGfActiveLeds, 'Color', purple);
grid;

linkaxes([s1c s1d s1e],'x');



%% Make the optics timetable OTT

Time = oT;

OTT = timetable(Time, oT_s, filtLed1, filtLed2, filtLed3, filtLed4, ...
    gfFiltLed1, gfFiltLed2, gfFiltLed3, gfFiltLed4, ...
    rssFiltLed1, rssFiltLed2, rssFiltLed3, rssFiltLed4, ...
    rssGfFiltLed1, rssGfFiltLed2, rssGfFiltLed3, rssGfFiltLed4, ...
    rssActiveLeds, rssAmbientLeds, rssGfActiveLeds, rssGfAmbientLeds);

%% create a decimated common working timetable

dE = downsample(E, 2);
dO = downsample(OTT, 5);

figure; 
plot(dO.Time, normalize(dO.filtLed2) * -1, 'Color', blue);
hold on;
plot(dO.Time, normalize(dO.filtLed3) * -1, 'Color', red);
plot(dE.Time, normalize(dE.bestImfBpOdba), 'Color', goldenrod);
hold off;
xlabel('Time, local');
ylabel('NormOptics (A.U.)');
title('Optical heart rate peaks with ODBA HR detections');
grid;


%% Run Pearson's correlation test against optical data

figCorrPlotFilt = figure('Position',[100 100 900 900], 'Color', white);
[rFiltLeds, pValsFiltLeds] = ...
    corrplot( [OTT.filtLed2 OTT.filtLed3 OTT.filtLed4, OTT.filtLed1], ...
    'type', 'Pearson', 'tail', 'right', 'varNames', ...
    {'1050nm','1200nm','Amb2','Amb1'} );

figCorrPlotGfFilt = figure('Position',[100 100 900 900], 'Color', white);
[rGfFiltLeds, pValsGfFiltLeds] = ...
    corrplot( [OTT.gfFiltLed2 OTT.gfFiltLed3 OTT.gfFiltLed4, OTT.gfFiltLed1], ...
    'type', 'Pearson', 'tail', 'right', 'varNames', ...
    {'1050nm','1200nm','Amb2','Amb1'} );


%% Okay, everything is in place for heart rate analysis...

%  Pause here and let the user select which analysis to run

fprintf('\n\n');
fprintf('Signal processing is complete on this epoch.\n');
fprintf('Press [enter] to continue analysis, or [q] to go manual.\n');
waitForInput = lower(input('Take your pick... (enter/q) [enter]: ','s'));
switch(waitForInput)
    case 'q'
        clc;
        fprintf('Bailing ftHeartRate_Optics.m around line 800...\n');
        fprintf('Pick it up from there!\n');
        return;
	otherwise
        fprintf('Continuing with analysis...\n');
end


%%

clc;

startCue = 1;
endCue = numel(gfFiltLed3);

oCue(1) = startCue;
oCue(2) = endCue;

markedUp = false;
markIndex = 3;

figMarkUp = figure('Position', ...
    [100 screenSize(4)-100 screenSize(3) screenSize(4) * .7], ...
    'Color', white);
figMarkUp.Name = 'Region selection for optical spectral analysis';

pMark2 = subplot(311);
pMark2.Position(1) = .10;
pMark2.Position(3) = .80;
plot(OTT.filtLed2, 'Color', blue);
hold on;
xline((startCue), 'Color', green', 'LineWidth', 1.4);
xline((endCue), 'Color', red, 'LineWidth', 1.4);
hold off;
xlabel('Samples, 250 Hz');
ylabel('AU');
grid;

pMark3 = subplot(312);
pMark3.Position(1) = .10;
pMark3.Position(3) = .80;
plot(OTT.filtLed3, 'Color', red);
hold on;
xline((startCue), 'Color', green', 'LineWidth', 1.4);
xline((endCue), 'Color', red, 'LineWidth', 1.4);
hold off;
xlabel('Samples, 250 Hz');
ylabel('AU');
grid;

pMark4 = subplot(313);
pMark4.Position(1) = .10;
pMark4.Position(3) = .80;
plot(OTT.filtLed4, 'Color', green);
hold on;
plot(OTT.filtLed1, 'Color', black);
xline((startCue), 'Color', green', 'LineWidth', 1.4);
xline((endCue), 'Color', red, 'LineWidth', 1.4);
hold off;
xlabel('Samples, 250 Hz');
ylabel('AU');
grid;


fprintf('Use this tool to mark-up the optical channels for regional analysis.\n');
fprintf('Commands are:\n');
fprintf('\t[m]ark - use to mark vertical line sample boundaries\n');
fprintf('\t[s]ave - once you have marked all regions, use this to sort and save\n');
fprintf('\t[q]uit - use this to quit without saving. oCue has unsorted boundaries\n');

while (markedUp == false)
    
    fprintf('Choose from one of these: [click] or [m] to mark | [s]ave | [q]uit\n');
    
    [x,~,button] = ginput(1);

    switch(button)
       case 'q'
           fprintf('Caught a [q]uit. Ending.\n');
           markedUp = true;
       case 's'
           fprintf('Caught a [s]ave. Saving and ending.\n');
           oCue = sort(oCue);
           markedUp = true;
        case ('m' || 1)
            
            % xAxis = gca;
            % t1 = num2ruler(x, xAxis.XAxis);
            
            t1 = x;
            
            figure(figMarkUp);
            subplot(pMark2);
            hold on;
            xline(t1, 'Color', maroon, 'LineWidth', 1.2);
            hold off;
            subplot(pMark3);
            hold on;
            xline(t1, 'Color', maroon, 'LineWidth', 1.2);
            hold off;            
            subplot(pMark4);
            hold on;
            xline(t1, 'Color', maroon, 'LineWidth', 1.2);
            hold off;            
            
            oCue(markIndex) = t1;
            markIndex = markIndex + 1;
            
            markedUp = false;
            
        otherwise
            fprintf('Not a valid input choice.\n');
            markedUp = false;
    end
    
end

fprintf('Ending mark-up. See OCUE.sample for region sample numbers.\n');

OCUE = struct;
OCUE.sample     = round(unique(oCue));
OCUE.kSample    = round( OCUE.sample / 2.5 );
OCUE.epochChoice = epochChoice;

%%

% build a model for attempting to auto-classify which ridge is the best fit
% for the consensusHr determined by SCG

useGaussTxt = 'Use Gaussian window smoothed optical data (y/n) [y]: ';
useGaussWindow = lower(input(useGaussTxt,'s'));
switch(useGaussWindow)
    case 'n'
        useGaussSmoothing = false;
    otherwise
        useGaussSmoothing = true;
end


meanConsensusHr = mean(E.consensusHr);
stdConsensusHr = std(E.consensusHr);
m196ConsensusHr = meanConsensusHr + ( -1.96 * stdConsensusHr);
p196ConsensusHr = meanConsensusHr + ( 1.96 * stdConsensusHr);

consensusHr_95percent = [m196ConsensusHr p196ConsensusHr];

regions = numel(OCUE.sample) - 1;

figRidgeExtract = figure('Position',[100 100 1400 800], 'Color', white);
figRidgeExtract.Name = 'Optical ridge extraction';

for i = 1:regions
    
   thisStart = round(OCUE.sample(i));
   thisEnd   = round(OCUE.sample(i+1))-1; 
   
   kStart    = round(OCUE.kSample(i));
   kEnd      = round(OCUE.kSample(i+1))-1;
   
   if (kStart == 0)
      kStart = 1; 
   end
   
   fprintf('Processing region %d: samples %d to %d...\n', ...
       i, thisStart, thisEnd);
   
   if (useGaussSmoothing)
       sampLed1 = gfFiltLed1(thisStart:thisEnd,1);
       sampLed2 = gfFiltLed2(thisStart:thisEnd,1);
       sampLed3 = gfFiltLed3(thisStart:thisEnd,1);
       sampLed4 = gfFiltLed4(thisStart:thisEnd,1);       
   else
       sampLed1 = filtLed1(thisStart:thisEnd,1);
       sampLed2 = filtLed2(thisStart:thisEnd,1);
       sampLed3 = filtLed3(thisStart:thisEnd,1);
       sampLed4 = filtLed4(thisStart:thisEnd,1);       
   end

   
   fprintf('\tsamples = %d | samples represent: %2.2f seconds\n', ...
       numel(sampLed3), numel(sampLed3)/oFs );
   
   [wsstThis1, fThis1] = wsst(sampLed1, oFs);
   [wsstThis2, fThis2] = wsst(sampLed2, oFs);
   [wsstThis3, fThis3] = wsst(sampLed3, oFs);
   [wsstThis4, fThis4] = wsst(sampLed4, oFs);
   
   numRidges = 2;

   [thisRidge1, iThisRidge1] = wsstridge(wsstThis1, 50, fThis1, ...
       'NumRidges', numRidges);
   [thisRidge2, iThisRidge2] = wsstridge(wsstThis2, 50, fThis2, ...
       'NumRidges', numRidges);
   [thisRidge3, iThisRidge3] = wsstridge(wsstThis3, 50, fThis3, ...
       'NumRidges', numRidges);
   [thisRidge4, iThisRidge4] = wsstridge(wsstThis4, 50, fThis4, ...
       'NumRidges', numRidges);   

   meanRidge1 = mean(thisRidge1(:,1:numRidges));
   stdRidge1  = std(thisRidge1(:,1:numRidges));
   p196Ridge1 = meanRidge1 + (1.96 * stdRidge1);
   m196Ridge1 = meanRidge1 + (-1.96 * stdRidge1);
   
   meanRidge2 = mean(thisRidge2(:,1:numRidges));
   stdRidge2  = std(thisRidge2(:,1:numRidges));
   p196Ridge2 = meanRidge2 + (1.96 * stdRidge2);
   m196Ridge2 = meanRidge2 + (-1.96 * stdRidge2);   
   
   meanRidge3 = mean(thisRidge3(:,1:numRidges));
   stdRidge3  = std(thisRidge3(:,1:numRidges));
   p196Ridge3 = meanRidge3 + (1.96 * stdRidge3);
   m196Ridge3 = meanRidge3 + (-1.96 * stdRidge3);
   
   meanRidge4 = mean(thisRidge4(:,1:numRidges));
   stdRidge4  = std(thisRidge4(:,1:numRidges));
   p196Ridge4 = meanRidge4 + (1.96 * stdRidge4);
   m196Ridge4 = meanRidge4 + (-1.96 * stdRidge4);
   
   ridge2_95percent = [m196Ridge2(1) p196Ridge2(1) ; ...
       m196Ridge2(2) p196Ridge2(2) ];

    % extract consensusHr from E

    thisConsensusHr = E.consensusHr(kStart:kEnd);

    figure(figRidgeExtract);

    sSCG = subplot(411);
    plot(thisConsensusHr, 'Color', maroon);
    xlabel('Samples, 100 Hz');
    ylabel('Frequency, Hz');
    title('Consensus SCG Heart Rate');
    grid;
    
    sO2 = subplot(412);
    plot(thisRidge2);
    xlabel('Samples, 250 Hz');
    ylabel('Frequency, Hz');
    title('1050 nm WSST ridges');
    grid;

    sO3 = subplot(413);
    plot(thisRidge3);
    xlabel('Samples, 250 Hz');
    ylabel('Frequency, Hz');
    title('1200 nm WSST ridges');
    grid;   

    sAmb = subplot(414);
    plot(thisRidge1);
    hold on;
    plot(thisRidge4);
    hold off;
    xlabel('Samples, 250 Hz');
    ylabel('Frequency, Hz');
    title('ambient 1 & 2 WSST ridges');
    grid;   
   
    linkaxes([sSCG sO2 sO3 sAmb],'x');
   
   
end


%% Generate some VMD IMFs from LED3 / optical channel 3

numImfs = 10;

numRows = numImfs / 2;
numCols = numImfs / numRows;

[imfFiltLed3, resFiltLed3, infoFiltLed3] = vmd(gfFiltLed3, 'NumIMFs', numImfs);

t = tiledlayout(numRows,numCols,'TileSpacing','compact','Padding','compact');
for n = 1:numImfs
    ax(n) = nexttile(t);
    plot(oT,imfFiltLed3(:,n)')
    xlim([oT(1) oT(end)])
    txt = ['IMF',num2str(n)];
    title(txt)
    xlabel('Time, local')
end
title(t,'Variational Mode Decomposition - 1200 nm optical channel (filtLed3)');


%% Generate some VMD IMFs from active optical channels

numImfs = 10;

numRows = numImfs / 2;
numCols = numImfs / numRows;

[imfActiveLeds, resActiveLeds, infoActiveLeds] = vmd(rssActiveLeds, ...
    'NumIMFs', numImfs);

t = tiledlayout(numRows,numCols,'TileSpacing','compact','Padding','compact');
for n = 1:numImfs
    ax(n) = nexttile(t);
    plot(oT,imfActiveLeds(:,n)')
    xlim([oT(1) oT(end)])
    txt = ['IMF',num2str(n)];
    title(txt)
    xlabel('Time, local')
end
title(t,'Variational Mode Decomposition')

%% Do peak-finding for bestImf and sumImf activeLeds (LED2 and LED3)

bestImfActiveLeds = imfActiveLeds(:,9);

% show this graphically in the time domain

figure;
findpeaks(bestImfActiveLeds, oFs, 'MinPeakDistance',0.5);

% ... but get results for computation in the sample domain

[pkBestImfActiveLeds, locBestImfActiveLeds, widBestImfActiveLeds, ...
    promBestImfActiveLeds] = findpeaks(bestImfActiveLeds, ...
    'MinPeakDistance', 125);

% sanity check this plot to compare time and sample domain results

figure;
plot(locBestImfActiveLeds, pkBestImfActiveLeds, 'bv');
grid;

% ... then do the same for sumImfOptics

% peaks for sumImf
sumImfActiveLeds = sum(imfActiveLeds(:,2:9),2);

figure;
findpeaks(sumImfActiveLeds, oFs, 'MinPeakDistance', 0.5);

[pkSumImfActiveLeds, locSumImfActiveLeds, widSumImfActiveLeds, ...
    promSumImfActiveLeds] = ...
    findpeaks(sumImfActiveLeds, 'MinPeakDistance', 125);

% sanity check this plot to compare time and sample domain results
figure;
plot(locSumImfActiveLeds, pkSumImfActiveLeds, 'bv');
grid;


%% Create timetables of these two peak estimations


HR_BESTIMF_ACTIVELEDS = timetable( oT(locBestImfActiveLeds), ...
    pkBestImfActiveLeds, widBestImfActiveLeds, promBestImfActiveLeds, ...
    'VariableNames', {'peaks','width','prominence'});

HR_SUMIMF_ACTIVELEDS = timetable( oT(locSumImfActiveLeds), ...
    pkSumImfActiveLeds, widSumImfActiveLeds, promSumImfActiveLeds, ...
    'VariableNames', {'peaks','width','prominence'});

% now spline interpolate in context of kT

bestImfActiveLedsSpline = interp1(HR_BESTIMF_ACTIVELEDS.Time, ...
    HR_BESTIMF_ACTIVELEDS.peaks, oT, 'spline');
bestImfActiveLedsRespCycle = csaps( oT_s(locBestImfActiveLeds), ...
    pkBestImfActiveLeds, 0.5, oT_s);

sumImfActiveLedsSpline = interp1(HR_SUMIMF_ACTIVELEDS.Time, ...
    HR_SUMIMF_ACTIVELEDS.peaks, oT, 'spline');
sumImfActiveLedsRespCycle = csaps(oT_s(locSumImfActiveLeds), ...
    pkSumImfActiveLeds, 0.5, oT_s);

% generate CUEs within the range of kT

cueTimes = CUE.Time(CUE.Time < oT(end) & CUE.Time > oT(1));
cueTypes = CUE.type(unique(cueTimes));

f1 = figure;
f1.Name = 'Respiration test based on IMF4 BpOdba';
figSize = f1.Position;
figSize(3) = screenSize(3);
figSize(4) = screenSize(4) * 0.5;
f1.Position = figSize;  

pPeaksBestImf = subplot(211);
plot(oT, bestImfActiveLeds, 'Color', red, 'LineWidth', 1.1);
hold on;
% yline(thresholdBestImfOdba, '--', 'Color', red);
plot(oT(locBestImfActiveLeds), pkBestImfActiveLeds, 'd', 'Color', blue, ...
    'MarkerFace', blue);
plot(oT, bestImfActiveLedsSpline, 'Color', blue);
plot(oT, bestImfActiveLedsRespCycle, 'Color', maroon, 'LineWidth', 1.1);
yPos = pPeaksBestImf.YLim(2) - (pPeaksBestImf.YLim(2) * 0.1);
numCuePlots = ftPlotCues(cueTimes, cueTypes, yPos);
hold off;
xlabel('Time, local');
ylabel('IHR bestImf');
grid;
legend('bestImfActiveLeds','Peaks','peakSpline','respCycle');

pPeaksSumImf = subplot(212);
plot(oT, sumImfActiveLeds, 'Color', green, 'LineWidth', 1.1);
hold on;
% yline(thresholdBestImfOdav, '--', 'Color', green);
plot(oT(locSumImfActiveLeds), pkSumImfActiveLeds, 'd', 'Color', blue, ...
    'MarkerFace', blue);
plot(oT, sumImfActiveLedsSpline, 'Color', green);
plot(oT, sumImfActiveLedsRespCycle, 'Color', maroon, 'LineWidth', 1.1);
yPos = pPeaksSumImf.YLim(2) - (pPeaksSumImf.YLim(2) * 0.1);
numCuePlots = ftPlotCues(cueTimes, cueTypes, yPos);
hold off;
xlabel('Time, local');
ylabel('IHR sumImf');
grid;
legend('sumImfActiveLeds','Peaks','peakSpline','respCycle');

linkaxes([pPeaksBestImf pPeaksSumImf],'x');

%% Try to make some frequency conversions of these peaks

bestImfActiveLedsLength = numel(bestImfActiveLeds);
sumImfActiveLedsLength = numel(sumImfActiveLeds);

peakBestImfActiveLedsElements = numel(locBestImfActiveLeds);
peakSumImfActiveLedsElements = numel(locSumImfActiveLeds);

peaksBestImfActiveLeds = zeros(bestImfActiveLedsLength,1);
peaksSumImfActiveLeds = zeros(sumImfActiveLedsLength,1);

for i = 1:peakBestImfActiveLedsElements

    thisStart = locBestImfActiveLeds(i);
    thisEnd = locBestImfActiveLeds(i) + round( widBestImfActiveLeds(i));
    if (thisEnd > bestImfActiveLedsLength)
        peaksBestImfActiveLeds(thisStart:bestImfActiveLedsLength,1) = 1;
    else
        peaksBestImfActiveLeds(thisStart:thisEnd,1) = 1;
    end
 
end

for i = 1:peakSumImfActiveLedsElements
   
    thisStart = locSumImfActiveLeds(i);
    thisEnd = locSumImfActiveLeds(i) + round( widSumImfActiveLeds(i));
    if (thisEnd > sumImfActiveLedsLength)
        peaksSumImfActiveLeds(thisStart:sumImfActiveLedsLength,1) = 1;    
    else
        peaksSumImfActiveLeds(thisStart:thisEnd,1) = 1;    
    end
    
end


figPeakFreq = figure('Position',[100 500 1400 1000]);

pPeak1 = subplot(211);
plot(oT, peaksBestImfActiveLeds, 'Color', red);
hold off;
xlabel('Time, seconds');
ylabel('HR Pulses');
title('BestImfActiveLeds - binary peak heart rate estimate');

pPeak2 = subplot(212);
plot(oT, peaksSumImfActiveLeds, 'Color', green);
hold off;
xlabel('Time, seconds');
ylabel('HR pulses');
title('SumImfActiveLeds - binary peak heart rate estimate');

linkaxes([pPeak1 pPeak2],'x');

%% compute WSSTs

% static code bestImfNumActiveLeds = 9 for now

bestImfNumActiveLeds = 9;

fprintf('Starting WSST computations...\n');

tic

fprintf('Creating WSSTs for band-passed active LED channels...\n');
[wsstActiveLeds, fActiveLeds] = wsst(rssActiveLeds, oFs);


fprintf('Creating WSSTs for sumImfActiveLeds...\n');
[wsstSumImfActiveLeds, fSumImfActiveLeds] = wsst(sumImfActiveLeds, oFs);


fprintf('Creating WSSTs for bestImf (%d) active LED channels...\n', ...
    bestImfNumActiveLeds);
[wsstBestImfActiveLeds, fBestImfActiveLeds] = wsst(bestImfActiveLeds, oFs);


fprintf('Creating WSSTs for peak-based sqwaves, sumImf & bestImf active LEDs ...\n');

[wsstPeakSqwaveSumImfActiveLeds, fPeakSqwaveSumImfActiveLeds] = ...
    wsst(peaksSumImfActiveLeds, oFs);

[wsstPeakSqwaveBestImfActiveLeds, fPeakSqwaveBestImfActiveLeds] = ...
    wsst(peaksBestImfActiveLeds, oFs);


fprintf('Completed WSST computations for all signals of interest...\n');
toc;


%% Plot all the bpOdba WSSTs

tic

fprintf('Plotting all the WSSTs...BE PATIENT! This can take a minute!\n');

figWsstPlots = figure('Position',[50 50 800 1200 ]);

p1wsst1 = subplot(511);
pcolor(oT, fActiveLeds, abs(wsstActiveLeds));
shading interp;
ylim([0 4]);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('WSST for band-passed active LED channels');

p1wsst2 = subplot(512);
pcolor(oT, fSumImfActiveLeds, abs(wsstSumImfActiveLeds));
shading interp;
ylim([0 4]);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('WSST for SumImfActiveLeds');

p1wsst3 = subplot(513);
pcolor(oT, fBestImfActiveLeds, abs(wsstBestImfActiveLeds));
shading interp;
ylim([0 4]);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('WSST for bestImfActiveLeds');

p1wsst4 = subplot(514);
pcolor(oT, fPeakSqwaveSumImfActiveLeds, abs(wsstPeakSqwaveSumImfActiveLeds));
shading interp;
ylim([0 4]);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('WSST for peak-based square-wave of sumImfActiveLeds');

p1wsst5 = subplot(515);
pcolor(oT, fPeakSqwaveBestImfActiveLeds, abs(wsstPeakSqwaveBestImfActiveLeds));
shading interp;
ylim([0 4]);
xlabel('Time, seconds');
ylabel('Frequency, Hz');
title('WSST for peak-based square-wave of bestImfActiveLeds');

fprintf('Completed plots of WSSTs for active LED channels.\n');

toc


%%

fprintf('Starting WSST ridge computations for ridge 1 and 2...\n');
fprintf('\t*** BE VERY PATIENT! This always takes 60+ seconds!\n');

tic

% first, do the band-passed active LED channels

fprintf('Computing WSST ridges for band-passed active optical channels...\n');

[ridgeActiveLeds, iRidgeActiveLeds] = wsstridge(wsstActiveLeds, 30, ...
    fActiveLeds, 'NumRidges', 2);    

% second, do sumImf for activeLeds

fprintf('Computing WSST ridges for sumImf of active LED channels)...\n');
[ridgeSumImfActiveLeds, iRidgeSumImfActiveLeds] = ...
    wsstridge(wsstSumImfActiveLeds, 30, fSumImfActiveLeds, 'NumRidges', 2);       

% third, do the bestImf activeLeds

fprintf('Computing WSST ridges for sumImf of active LED channels)...\n');
[ridgeBestImfActiveLeds, iRidgeBestImfActiveLeds] = ...
    wsstridge(wsstBestImfActiveLeds, 30, fBestImfActiveLeds, 'NumRidges', 2);       

% fourth, do peak-based square-wave of sumImf of activeLeds

fprintf('Computing WSST ridges for peak-based sumImfActiveLeds square wave...\n');
[ridgePeakSqwaveSumImfActiveLeds, iRidgePeakSqwaveSumImfActiveLeds] = ...
    wsstridge(wsstPeakSqwaveSumImfActiveLeds, 30, ...
    fPeakSqwaveSumImfActiveLeds, 'NumRidges', 2);


% fifth, and finally, do peak-based square-wave bestImf bpOdba & bpOdav

fprintf('Computing WSST ridges for peak-based bestImfActiveLeds square wave...\n');
[ridgePeakSqwaveBestImfActiveLeds, iRidgePeakSqwaveBestImfActiveLeds] = ...
    wsstridge(wsstPeakSqwaveBestImfActiveLeds, 30, ...
    fPeakSqwaveBestImfActiveLeds, 'NumRidges', 2);
fprintf('Completed WSST ridge computations for all signals of interest...\n');

toc;   


%%


%% examine first and second ridges for these respective predictors

figRidges = figure('Position', [ 0 0 screenSize(3)/2 screenSize(4)-100 ]);

% wsstridges for band-passed active optical channels

pRidge1 = subplot(511);
plot(oT, ridgeActiveLeds(:,1), 'Color', blue);
hold on;
plot(oT, ridgeActiveLeds(:,2), 'Color', red);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
ylim([0 4]);
title('WSST ridges for band-passed active optical channels');
grid;
legend('r1', 'r2');

% wsstridges for sumImf of band-passed active optical channels

pRidge2 = subplot(512);
plot(oT, ridgeSumImfActiveLeds(:,1), 'Color', blue);
hold on;
plot(oT, ridgeSumImfActiveLeds(:,2), 'Color', red);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
ylim([0 4]);
title('WSST ridges for sumImf active optical channels');
grid;
legend('r1', 'r2');

% wsstridges for bestImf of band-passed active optical channels

pRidge3 = subplot(513);
plot(oT, ridgeBestImfActiveLeds(:,1), 'Color', blue);
hold on;
plot(oT, ridgeBestImfActiveLeds(:,2), 'Color', red);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
ylim([0 4]);
title('WSST ridges for bestImf active optical channels');
grid;
legend('r1', 'r2');

% wsstridges for square wave of sumImf of band-passed active optical channels

pRidge4 = subplot(514);
plot(oT, ridgePeakSqwaveSumImfActiveLeds(:,1), 'Color', blue);
hold on;
plot(oT, ridgePeakSqwaveSumImfActiveLeds(:,2), 'Color', red);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
ylim([0 4]);
title('WSST ridges for peak-based square wave of sumImf active optical channels');
grid;
legend('r1', 'r2');

% wsstridges for square wave of bestImf of band-passed active optical channels

pRidge5 = subplot(515);
plot(oT, ridgePeakSqwaveBestImfActiveLeds(:,1), 'Color', blue);
hold on;
plot(oT, ridgePeakSqwaveBestImfActiveLeds(:,2), 'Color', red);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
ylim([0 4]);
title('WSST ridges for peak-based square wave of bestImf active optical channels');
grid;
legend('r1', 'r2');

linkaxes([pRidge1 pRidge2 pRidge3 pRidge4 pRidge5],'x');


%% Have the user confirm which ridge corresponds with HR ridges

%   This step will ideally be automated eventually, but for now it's tricky
%   for the ridge detection associated with heart rate to be automatically
%   detected in a correct way. So prompt user to specify which ridge is the
%   HR ridge, if any, to use for consensus HR estimation


confirmed = false;

while(~confirmed)

bestHrRidgeActiveLeds                   = 100;
bestHrRidgeSumImfActiveLeds             = 100;
bestHrRidgeBestImfActiveLeds            = 100;
bestHrRidgePeakSqwaveSumImfActiveLeds	= 100;
bestHrRidgePeakSqwaveBestImfActiveLeds  = 100;


fprintf('For each of the following, choose the ridge that best fits HR...\n');
fprintf('Valid ridge values are:\n');
fprintf('\t0 = neither ridge\n');
fprintf('\t1 = ridge 1\n');
fprintf('\t2 = ridge 2\n');

% band-passed active optical channels

while (bestHrRidgeActiveLeds == 100)
    thisTxt = 'Best HR ridge for active optical channels: ';
    bestHrRidgeActiveLeds = str2double(input(thisTxt,'s'));
    switch(bestHrRidgeActiveLeds)
        case 0
            fprintf('\tYou selected neither ridge for activeLeds.\n');
        case 1
            fprintf('\tYou selected ridge 1 for activeLeds.\n')
        case 2
            fprintf('\tYou selected ridge 2 for activeLeds.\n')
        otherwise
            fprintf('\tYou selected an invalid option. Try again.\n');
            bestHrRidgeActiveLeds = 100;
    end
end

% sumImfActiveLeds

while (bestHrRidgeSumImfActiveLeds == 100)
    bestHrRidgeSumImfActiveLeds = str2double(input('Best HR ridge for sumImfActiveLeds: ','s'));
    switch(bestHrRidgeSumImfActiveLeds)
        case 0
            fprintf('\tYou selected neither ridge for sumImfActiveLeds.\n');
        case 1
            fprintf('\tYou selected ridge 1 for sumImfActiveLeds.\n')
        case 2
            fprintf('\tYou selected ridge 2 for sumImfActiveLeds.\n')
        otherwise
            fprintf('\tYou selected an invalid option. Try again.\n');
            bestHrRidgeSumImfActiveLeds = 100;
    end
end


% bestImfActiveLeds

while (bestHrRidgeBestImfActiveLeds == 100)
    bestHrRidgeBestImfActiveLeds = str2double(input('Best HR ridge for bestImfActiveLeds: ','s'));
    switch(bestHrRidgeBestImfActiveLeds)
        case 0
            fprintf('\tYou selected neither ridge for bestImfActiveLeds.\n');
        case 1
            fprintf('\tYou selected ridge 1 for bestImfActiveLeds.\n')
        case 2
            fprintf('\tYou selected ridge 2 for bestImfActiveLeds.\n')
        otherwise
            fprintf('\tYou selected an invalid option. Try again.\n');
            bestHrRidgeBestImfActiveLeds = 100;
    end
end

% peakSqwaveSumImfActiveLeds

while (bestHrRidgePeakSqwaveSumImfActiveLeds == 100)
    bestHrRidgePeakSqwaveSumImfActiveLeds = ...
        str2double(input('Best HR ridge for peakSqwaveBestImfBpOdba: ','s'));
    switch(bestHrRidgePeakSqwaveSumImfActiveLeds)
        case 0
            fprintf('\tYou selected neither ridge for peakSqwaveBestImfBpOdba.\n');
        case 1
            fprintf('\tYou selected ridge 1 for peakSqwaveBestImfBpOdba.\n')
        case 2
            fprintf('\tYou selected ridge 2 for peakSqwaveBestImfBpOdba.\n')
        otherwise
            fprintf('\tYou selected an invalid option. Try again.\n');
            bestHrRidgePeakSqwaveSumImfActiveLeds = 100;
    end
end

% peakSqwaveBestImfActiveLeds

while (bestHrRidgePeakSqwaveBestImfActiveLeds == 100)
    bestHrRidgePeakSqwaveBestImfActiveLeds = ...
        str2double(input('Best HR ridge for peakSqwaveBestImfBpOdav: ','s'));
    switch(bestHrRidgePeakSqwaveBestImfActiveLeds)
        case 0
            fprintf('\tYou selected neither ridge for peakSqwaveBestImfBpOdav.\n');
        case 1
            fprintf('\tYou selected ridge 1 for peakSqwaveBestImfBpOdav.\n')
        case 2
            fprintf('\tYou selected ridge 2 for peakSqwaveBestImfBpOdav.\n')
        otherwise
            fprintf('\tYou selected an invalid option. Try again.\n');
            bestHrRidgePeakSqwaveBestImfActiveLeds = 100;
    end
end

fprintf('Best HR ridge summary:\n');

fprintf('\tBand-passed active optical channels: ridge %d\n', ...
    bestHrRidgeActiveLeds);

fprintf('\tSumImf active optical channels: ridge %d\n', ...
    bestHrRidgeSumImfActiveLeds);

fprintf('\tSumImf active optical channels: ridge %d\n', ...
    bestHrRidgeBestImfActiveLeds);

fprintf('\tPeak-based square wave sumImfActiveLeds: ridge %d\n', ...
    bestHrRidgePeakSqwaveSumImfActiveLeds);

fprintf('\tPeak-based square wave bestImfActiveLeds: ridge %d\n', ...
    bestHrRidgePeakSqwaveBestImfActiveLeds);


confirmStr = lower(input('That all look good? (y/n): ','s'));
switch(confirmStr)
    case 'y'
        fprintf('Good. Continuing!\n');
        confirmed = true;
    otherwise
        fprintf('Okay... redoing best ridge selection.\n');
        confirmed = false;
end

end

%% Create some summary statistics

meanActiveLedsRidges = [ mean(ridgeActiveLeds(:,1)), mean(ridgeActiveLeds(:,2)) ];
stdActiveLedsRidges = [ std(ridgeActiveLeds(:,1)), std(ridgeActiveLeds(:,2)) ];

meanSumImfActiveLedsRidges = [ mean(ridgeSumImfActiveLeds(:,1)), ...
    mean(ridgeSumImfActiveLeds(:,2)) ];
stdSumImfActiveLedsRidges = [ std(ridgeSumImfActiveLeds(:,1)), ...
    std(ridgeSumImfActiveLeds(:,2)) ];

meanBestImfActiveLedsRidges = [ mean(ridgeBestImfActiveLeds(:,1)), ...
    mean(ridgeBestImfActiveLeds(:,2)) ];
stdBestImfActiveLedsRidges = [ std(ridgeBestImfActiveLeds(:,1)), ...
    std(ridgeBestImfActiveLeds(:,2)) ];

meanPeakSqwaveSumImfActiveLedsRidges = [ mean(ridgePeakSqwaveSumImfActiveLeds(:,1)), ...
    mean(ridgePeakSqwaveSumImfActiveLeds(:,2)) ];
stdPeakSqwaveSumImfActiveLedsRidges = [ std(ridgePeakSqwaveSumImfActiveLeds(:,1)), ...
    std(ridgePeakSqwaveSumImfActiveLeds(:,2)) ];

meanPeakSqwaveBestImfActiveLedsRidges = [ mean(ridgePeakSqwaveBestImfActiveLeds(:,1)), ...
    mean(ridgePeakSqwaveBestImfActiveLeds(:,2)) ];
stdPeakSqwaveBestImfActiveLedsRidges = [ std(ridgePeakSqwaveBestImfActiveLeds(:,1)), ...
    std(ridgePeakSqwaveBestImfActiveLeds(:,2)) ];

%% Make histograms

figHrFreq = figure('Position',[0 0 screenSize(3)/3 screenSize(4) ], ...
    'Color', white, 'NumberTitle', 'off');

sHrF1 = subplot(511);
if (bestHrRidgeActiveLeds ~=0)
    histogram( ridgeActiveLeds(:,bestHrRidgeActiveLeds), 'Normalization', 'pdf', ...
        'FaceColor', red, 'EdgeColor', black, 'Orientation', 'horizontal');
    hold on;
    yline(meanActiveLedsRidges(bestHrRidgeActiveLeds), 'Color', blue, ...
        'LineWidth', 1.2);
    yline(meanActiveLedsRidges(bestHrRidgeActiveLeds) + ...
        1.96 * stdActiveLedsRidges(bestHrRidgeActiveLeds), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    yline(meanActiveLedsRidges(bestHrRidgeActiveLeds) + ...
        -1.96 * stdActiveLedsRidges(bestHrRidgeActiveLeds), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    xMean = mean(sHrF1.XAxis.Limits);
    xLoc = xMean + xMean * 0.5;
    yMean = mean(sHrF1.YAxis.Limits);
    yLoc = yMean + yMean * 0.375;
	ySpacing = range(sHrF1.YAxis.Limits) / 12;
    meanStr = sprintf('%1.2f Hz', meanActiveLedsRidges(bestHrRidgeActiveLeds));
    text(xLoc, yLoc, ['{\mu}: ' meanStr],'FontSize',12);
    stdStr = sprintf('%1.2f Hz', stdActiveLedsRidges(bestHrRidgeActiveLeds));
    text(xLoc, yLoc - ySpacing, ['{\sigma}: ' stdStr],'FontSize',12);
    hold off;
    ylabel('Frequency, Hz');
    grid;
else
    xTextLoc = mean(sHrF1.XAxis.Limits) / 2;
    yTextLoc = mean(sHrF1.YAxis.Limits);
    text(xTextLoc, yTextLoc, 'This ridge was not selected for ensemble averaging.', 'FontSize', 12);
end

sHrF2 = subplot(512);
if (bestHrRidgeSumImfActiveLeds ~=0)
    histogram( ridgeSumImfActiveLeds(:,bestHrRidgeSumImfActiveLeds), 'Normalization', 'pdf', ...
        'FaceColor', red, 'EdgeColor', black, 'Orientation', 'horizontal');
    hold on;
    yline(meanSumImfActiveLedsRidges(bestHrRidgeSumImfActiveLeds), 'Color', blue, ...
        'LineWidth', 1.2);
    yline(meanSumfImfActiveLedsRidges(bestHrRidgeSumImfActiveLeds) + ...
        1.96 * stdSumImfActiveLedsRidges(bestHrRidgeSumImfActiveLeds), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    yline(meanSumImfActiveLedsRidges(bestHrRidgeSumImfActiveLeds) + ...
        -1.96 * stdSumImfActiveLedsRidges(bestHrRidgeSumImfActiveLeds), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    xMean = mean(sHrF2.XAxis.Limits);
    xLoc = xMean + xMean * 0.5;
    yMean = mean(sHrF2.YAxis.Limits);
    yLoc = yMean + yMean * 0.375;
	ySpacing = range(sHrF2.YAxis.Limits) / 12;    
    meanStr = sprintf('%1.2f Hz', meanSumImfActiveLedsRidges(bestHrRidgeSumImfActiveLeds));
    text(xLoc, yLoc, ['{\mu}: ' meanStr],'FontSize',12);
    stdStr = sprintf('%1.2f Hz', stdSumImfActiveLedsRidges(bestHrRidgeSumImfActiveLeds));
    text(xLoc, yLoc - ySpacing, ['{\sigma}: ' stdStr],'FontSize',12);
    hold off;
    ylabel('Frequency, Hz');
    grid;
else
    xTextLoc = mean(sHrF2.XAxis.Limits) / 2;
    yTextLoc = mean(sHrF2.YAxis.Limits);
    text(xTextLoc, yTextLoc, 'This ridge was not selected for ensemble averaging.', 'FontSize', 12);
end

sHrF3 = subplot(513);
if (bestHrRidgeBestImfActiveLeds ~=0)
    histogram( ridgeSumImfActiveLeds(:,bestHrRidgeBestImfActiveLeds), 'Normalization', 'pdf', ...
        'FaceColor', red, 'EdgeColor', black, 'Orientation', 'horizontal');
    hold on;
    yline(meanBestImfActiveLedsRidges(bestHrRidgeBestImfActiveLeds), 'Color', blue, ...
        'LineWidth', 1.2);
    yline(meanBestImfActiveLedsRidges(bestHrRidgeBestImfActiveLeds) + ...
        1.96 * stdBestImfActiveLedsRidges(bestHrRidgeBestImfActiveLeds), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    yline(meanBestImfActiveLedsRidges(bestHrRidgeBestImfActiveLeds) + ...
        -1.96 * stdBestImfActiveLedsRidges(bestHrRidgeBestImfActiveLeds), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    xMean = mean(sHrF3.XAxis.Limits);
    xLoc = xMean + xMean * 0.5;
    yMean = mean(sHrF3.YAxis.Limits);
    yLoc = yMean + yMean * 0.375;
	ySpacing = range(sHrF3.YAxis.Limits) / 12;  
    meanStr = sprintf('%1.2f Hz', meanBestImfActiveLedsRidges(bestHrRidgeBestImfActiveLeds));
    text(xLoc, yLoc, ['{\mu}: ' meanStr],'FontSize',12);
    stdStr = sprintf('%1.2f Hz', stdBestImfActiveLedsRidges(bestHrRidgeBestImfActiveLeds));
    text(xLoc, yLoc - ySpacing, ['{\sigma}: ' stdStr],'FontSize',12);
    hold off;
    ylabel('Frequency, Hz');
    grid;
else
    xTextLoc = mean(sHrF3.XAxis.Limits) / 2;
    yTextLoc = mean(sHrF3.YAxis.Limits);
    text(xTextLoc, yTextLoc, 'This ridge was not selected for ensemble averaging.', 'FontSize', 12);
end

sHrF4 = subplot(514);
if (bestHrRidgePeakSqwaveSumImfActiveLeds ~=0)
    histogram( ridgePeakSqwaveSumImfActiveLeds(:,bestHrRidgePeakSqwaveSumImfActiveLeds), ...
        'Normalization', 'pdf', 'FaceColor', red, ...
        'EdgeColor', black, 'Orientation', 'horizontal');
    hold on;
    yline(meanPeakSqwaveSumImfActiveLedsRidges(bestHrRidgePeakSqwaveSumImfActiveLeds), ...
        'Color', blue, 'LineWidth', 1.2);
    yline(meanPeakSqwaveSumImfActiveLedsRidges(bestHrRidgePeakSqwaveSumImfActiveLeds) + ...
        1.96 * stdSumImfActiveLedsRidges(bestHrRidgePeakSqwaveSumImfActiveLeds), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    yline(meanPeakSqwaveSumImfActiveLedsRidges(bestHrRidgePeakSqwaveSumImfActiveLeds) + ...
        -1.96 * stdSumImfActiveLedsRidges(bestHrRidgePeakSqwaveSumImfActiveLeds), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    
    xMean = mean(sHrF4.XAxis.Limits);
    xLoc = xMean + xMean * 0.5;
    yMean = mean(sHrF4.YAxis.Limits);
    yLoc = yMean + yMean * 0.375;   
    ySpacing = range(sHrF4.YAxis.Limits) / 12;
    meanStr = sprintf('%1.2f Hz', meanPeakSqwaveSumImfActiveLedsRidges(bestHrRidgePeakSqwaveSumImfActiveLeds));
    text(xLoc, yLoc, ['{\mu}: ' meanStr],'FontSize',12);
    stdStr = sprintf('%1.2f Hz', stdPeakSqwaveSumImfActiveLedsRidges(bestHrRidgePeakSqwaveSumImfActiveLeds));
    text(xLoc, yLoc - ySpacing, ['{\sigma}: ' stdStr],'FontSize',12);
    hold off;
    ylabel('Frequency, Hz');
    grid;
else
    xTextLoc = mean(sHrF4.XAxis.Limits) / 2;
    yTextLoc = mean(sHrF4.YAxis.Limits);
    text(xTextLoc, yTextLoc, 'This ridge was not selected for ensemble averaging.', 'FontSize', 12);
end

sHrF5 = subplot(515);
if (bestHrRidgePeakSqwaveBestImfActiveLeds ~=0)
    histogram( ridgePeakSqwaveBestImfActiveLeds(:,bestHrRidgePeakSqwaveBestImfActiveLeds), ...
        'Normalization', 'pdf', 'FaceColor', red, ...
        'EdgeColor', black, 'Orientation', 'horizontal');
    hold on;
    yline(meanPeakSqwaveBestImfActiveLedsRidges(bestHrRidgePeakSqwaveBestImfActiveLeds), ...
        'Color', blue, 'LineWidth', 1.2);
    yline(meanPeakSqwaveBestImfActiveLedsRidges(bestHrRidgePeakSqwaveBestImfActiveLeds) + ...
        1.96 * stdBestImfActiveLedsRidges(bestHrRidgePeakSqwaveBestImfActiveLeds), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    yline(meanPeakSqwaveBestImfActiveLedsRidges(bestHrRidgePeakSqwaveBestImfActiveLeds) + ...
        -1.96 * stdPeakSqwaveBestImfActiveLedsRidges(bestHrRidgePeakSqwaveBestImfActiveLeds), ...
        '--', 'Color', blue, 'LineWidth', 1.1);

    xMean = mean(sHrF5.XAxis.Limits);
    xLoc = xMean + xMean * 0.5;
    yMean = mean(sHrF5.YAxis.Limits);
    yLoc = yMean + yMean * 0.375;
	ySpacing = range(sHrF5.YAxis.Limits) / 12;
    meanStr = sprintf('%1.2f Hz', meanPeakSqwaveBestImfActiveLedsRidges(bestHrRidgePeakSqwaveBestImfActiveLeds));
    text(xLoc, yLoc, ['{\mu}: ' meanStr],'FontSize',12);
    stdStr = sprintf('%1.2f Hz', stdPeakSqwaveBestImfActiveLedsRidges(bestHrRidgePeakSqwaveBestImfActiveLeds));
    text(xLoc, yLoc - ySpacing, ['{\sigma}: ' stdStr],'FontSize',12);
    hold off;
    ylabel('Frequency, Hz');
    grid;
else
    xTextLoc = mean(sHrF5.XAxis.Limits) / 2;
    yTextLoc = mean(sHrF5.YAxis.Limits);
    text(xTextLoc, yTextLoc, 'This ridge was not selected for ensemble averaging.', 'FontSize', 12);
end


%% create smoothed versions of the peak sqwave ridges

peakSqwaveSumImfRidge = smooth( ...
    ridgePeakSqwaveSumImfActiveLeds(:,bestHrRidgePeakSqwaveSumImfActiveLeds), ...
    125);

peakSqwaveBestImfRidge = smooth( ...,
    ridgePeakSqwaveBestImfActiveLeds(:,bestHrRidgePeakSqwaveBestImfActiveLeds), ...
    125);


figure;

p5a = subplot(211);
plot(oT, peakSqwaveSumImfRidge, 'Color', green);
hold on;
plot(oT, peakSqwaveBestImfRidge, 'Color', blue);
plot(E.Time, E.consensusHr, 'Color', red, 'LineWidth', 1.1);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
grid;
legend('sumImfPeaks','bestImfPeaks','kinematics');

p5b = subplot(212);
plot(oT, peakSqwaveBestImfRidge, 'Color', blue);
hold on;
plot(E.Time, E.consensusHr, 'Color', red, 'LineWidth', 1.1);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
grid;
legend('optics','kinematics');

linkaxes([p5a p5b],'x');

%% decimate to compare optics and kinematic "truth" with Bland-Altman

% this should achieve a common 25 Hz 

dPeakSumImfRidge = decimate(peakSqwaveSumImfRidge, 10);
dPeakBestImfRidge = decimate(peakSqwaveBestImfRidge, 10);
dConsensusHr = decimate(E.consensusHr, 4);
dT = downsample(oT, 10);
dT_S = decimate(oT_s, 10);

cueTimes = CUE.Time(CUE.Time < dT(end) & CUE.Time > dT(1));
cueTypes = CUE.type(unique(cueTimes));

figAll = figure('Position', [100 100 800 500], 'Color', white);

plot(dT, dPeakSumImfRidge, 'g--');
hold on;
plot(dT, dPeakBestImfRidge, 'b-');
plot(dT, dConsensusHr, 'r-');
figInfo = gca;
yLims = figInfo.YAxis.Limits;
yRange = range(figInfo.YAxis.Limits);
yPos = yLims(2) - yRange/20;
numCuePlots = ftPlotCues(cueTimes, cueTypes, yPos);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Comparison of optical heart rate estimate with SCG ground-truth');
grid;
legend('BestImf optics','SCG');





%%

figSizeBA = [ 100 100 600 600 ];

fig1 = figure('Position', figSizeBA, 'Color', white, 'NumberTitle', 'off');
[ rpc1, fig1, stats1 ]  = BlandAltman(fig1, ...
    dPeakSumImfRidge, dConsensusHr, ...
    {'peakSqwaveSumImf','kinematicHR'}, ...
    'Bland-Altman test for peakSqwave sum vs. SCG', ...
    {'peakSumIMF','SCG'}, ...
    'corrInfo', {'eq','r','r2','p','SSE','RMSE','n'}, ...
    'baStatsMode','Non-parametric');

fig2 = figure('Position', figSizeBA, 'Color', white, 'NumberTitle', 'off');
[ rpc2, fig2, stats2 ]  = BlandAltman(fig2, ...
    (dPeakBestImfRidge).^2, (dConsensusHr).^2, ...
    {'peakSqwaveBestImf','kinematicHR'}, ...
    'Bland-Altman test for peakSqwave best IMF vs. SCG', ...
    {'peakBestIMF','SCG'}, ...
    'corrInfo', {'eq','r','r2','p','SSE','RMSE','n'}, ...
    'baStatsMode','Non-parametric');



%%


meanConsensusHr = mean(E.consensusHr) * 60;
stdConsensusHr = std(E.consensusHr) * 60;

meanOpticsHr = mean(dPeakBestImfRidge) * 60;
stdOpticsHr = std(dPeakBestImfRidge) * 60;

figHist = figure('Position', [100 100 500 500], 'Color', white);

histogram(dPeakBestImfRidge.*60, 'FaceColor', blue, ...
    'EdgeColor', black, 'Orientation', 'horizontal');
hold on;
histogram(dConsensusHr.*60, 'FaceColor', red, ...
    'EdgeColor', black, 'Orientation', 'horizontal');

yline(meanConsensusHr, 'Color', red, 'LineWidth', 1.2);
yline(meanConsensusHr + 1.96 * stdConsensusHr, '--', 'Color', red);
yline(meanConsensusHr + -1.96 * stdConsensusHr, '--', 'Color', red);

yline(meanOpticsHr, 'Color', blue, 'LineWidth', 1.2);
yline(meanOpticsHr + 1.96 * stdOpticsHr, '--', 'Color', blue);
yline(meanOpticsHr + -1.96 * stdOpticsHr, '--', 'Color', blue);

figInfo = gca;

xLims = figInfo.XAxis.Limits;
yLims = figInfo.YAxis.Limits;
xRange = range(figInfo.XAxis.Limits);
yRange = range(figInfo.YAxis.Limits);

xLoc1 = mean(xLims) + xRange * 0.2; 
yLoc1 = yLims(2) - yLims(2) * .05;
ySpacing = yRange / 30;

xLoc2 = mean(xLims) + xRange * 0.2;
yLoc2 = yLims(1) + yLims(1) * 0.05;

meanStr = sprintf('%1.2f BPM', meanConsensusHr);
text(xLoc1, yLoc1, ['{\mu} SCG: ' meanStr],'FontSize',12);
stdStr = sprintf('%1.2f BPM', stdConsensusHr);
text(xLoc1, yLoc1 - ySpacing, ['{\sigma} SCG: ' stdStr],'FontSize',12);

meanStr = sprintf('%1.2f BPM', meanOpticsHr);
text(xLoc2, yLoc2 + ySpacing, ['{\mu} OCG: ' meanStr],'FontSize',12);
stdStr = sprintf('%1.2f BPM', stdOpticsHr);
text(xLoc2, yLoc2, ['{\sigma} OCG: ' stdStr],'FontSize',12);

hold off;
ylabel('Frequency, Hz');
grid;
legend('','');


%% Ask user if they are satisfied with analysis results for save to RAW

prepToSave = lower(input('Make optics HR timetables and peak structs? (y/n): ','s'));

switch(prepToSave)
    
    case 'y'

        fprintf('Building HR timetable and peak structs...\n');
        
        Time = oT;

        % Build a generic timetable with all necessary constituent parts

        TT = timetable(Time, filtLed1, filtLed2, filtLed3, filtLed4);

        % Build the correct timetable and structs based on epochChoice

        switch(epochChoice)

            case 1

                fprintf('Creating E1 timetable and peak structures...\n');
                E1_OPTICS = TT;
                E1_PEAKS_SUMIMF_ACTIVELEDS = HR_SUMIMF_ACTIVELEDS;
                E1_PEAKS_BESTIMF_ACTIVELEDS = HR_BESTIMF_ACTIVELEDS;


            case 2

                fprintf('Creating E2 timetable and peak structures...\n');
                E2_OPTICS = TT;
                E2_PEAKS_SUMIMF_ACTIVELEDS = HR_SUMIMF_ACTIVELEDS;
                E2_PEAKS_BESTIMF_ACTIVELEDS = HR_BESTIMF_ACTIVELEDS;


            case 3

                fprintf('Creating E3 timetable and peak structures...\n');
                E3_OPTICS = TT;
                E3_PEAKS_SUMIMF_ACTIVELEDS = HR_SUMIMF_ACTIVELEDS;
                E3_PEAKS_BESTIMF_ACTIVELEDS = HR_BESTIMF_ACTIVELEDS;

            case 4

                fprintf('Creating E4 timetable and peak structures...\n');
                E4_OPTICS = TT;
                E4_PEAKS_SUMIMF_ACTIVELEDS = HR_SUMIMF_ACTIVELEDS;
                E4_PEAKS_BESTIMF_ACTIVELEDS = HR_BESTIMF_ACTIVELEDS;


            case 5

                fprintf('Creating E5 timetable and peak structures...\n');
                E5_OPTICS = TT;
                E5_PEAKS_SUMIMF_ACTIVELEDS = HR_SUMIMF_ACTIVELEDS;
                E5_PEAKS_BESTIMF_ACTIVELEDS = HR_BESTIMF_ACTIVELEDS;

            otherwise
                fprintf('Epoch choice %d is invalid. Bailing!\n');
                return;

        end
        
        fprintf('Preparing to append HR timetable and peak structs to RAW file:\n');
        fprintf('\t --> %s\n', fullRawFileName );
        
        if (exist(fullRawFileName, 'file'))

            fprintf('RAW file exists... ');

            prompt = 'OK to append new data? y/n [n]: ';
            str = input(prompt,'s');
            if isempty(str)
                str = 'n';
            end

            if (strcmp(str, 'y'))
               fprintf('\nOver-writing and appending new data to RAW file...\n');
               switch(epochChoice)
                   case 1
                       save(fullRawFileName, 'E1_KINEMATICS', ...
                           'E1_PEAKS_ODBA', 'E1_PEAKS_ODAV', '-append');
                   case 2
                       save(fullRawFileName, 'E2_KINEMATICS', ...
                           'E2_PEAKS_ODBA', 'E2_PEAKS_ODAV', '-append');
                   case 3
                       save(fullRawFileName, 'E3_KINEMATICS', ...
                           'E3_PEAKS_ODBA', 'E3_PEAKS_ODAV', '-append');
                   case 4
                       save(fullRawFileName, 'E4_KINEMATICS', ...
                           'E4_PEAKS_ODBA', 'E4_PEAKS_ODAV', '-append');
                   case 5
                       save(fullRawFileName, 'E5_KINEMATICS', ...
                           'E5_PEAKS_ODBA', 'E5_PEAKS_ODAV', '-append');
                   otherwise
                       fprintf('There is epochChoice specified. Bailing!\n');
                       return;
               end

            end 

        else

            fprintf('Raw file does not exist... this should never happen... Bail!\n');
            return;

        end        
        
    otherwise
        
        fprintf('Skipping append+save to RAW.\n');

end


%% if E1-E5 are all present, make a stitch plot of consensusHr 

if ( exist('E1_KINEMATICS','var') & exist('E2_KINEMATICS','var') & ...
    exist('E3_KINEMATICS','var') & exist('E4_KINEMATICS','var') & ...
    exist('E5_KINEMATICS','var') )

figHrStitch = figure('Position', [100 100 1200 500], 'Color', white); 

plot(E1_KINEMATICS.Time, E1_KINEMATICS.consensusHr, 'Color', red, ...
    'LineWidth', 1.2); 
hold on;
plot(E2_KINEMATICS.Time, E2_KINEMATICS.consensusHr, 'Color', black, ...
    'LineWidth', 0.8, 'LineStyle', '--');
plot(E3_KINEMATICS.Time, E3_KINEMATICS.consensusHr, 'Color', blue, ...
    'LineWidth', 1.2);
plot(E4_KINEMATICS.Time, E4_KINEMATICS.consensusHr, 'Color', black, ...
    'LineWidth', 0.8, 'LineStyle', '--');
plot(E5_KINEMATICS.Time, E5_KINEMATICS.consensusHr, 'Color', green, ...
    'LineWidth', 1.2);
hold off; 
xlabel('Time, local');
ylabel('Heart rate, Hz');

grid;

figHrStitch = figure('Position', [100 100 1200 500], 'Color', white); 

plot(E1_KINEMATICS.Time, E1_KINEMATICS.consensusHr .* 60, 'Color', red, ...
    'LineWidth', 1.2); 
hold on;
plot(E2_KINEMATICS.Time, E2_KINEMATICS.consensusHr .* 60, 'Color', purple, ...
    'LineWidth', 0.8, 'LineStyle', '--');
plot(E3_KINEMATICS.Time, E3_KINEMATICS.consensusHr .* 60, 'Color', blue, ...
    'LineWidth', 1.2);
plot(E4_KINEMATICS.Time, E4_KINEMATICS.consensusHr .* 60, 'Color', cyan, ...
    'LineWidth', 0.8, 'LineStyle', '--');
plot(E5_KINEMATICS.Time, E5_KINEMATICS.consensusHr .* 60, 'Color', green, ...
    'LineWidth', 1.2);
hold off; 
xlabel('Time, local');
ylabel('Heart rate, BPM');
titleStr = sprintf('Heart rate profile for %s - %s', ...
    CUE.Properties.CustomProperties.SUBJECT.id, ...
    CUE.Properties.CustomProperties.tag);
title(titleStr,'Interpreter','none');
grid;
legend('pre-apnea','transition1','apnea phase',...
    'transition2','post-apnea phase');

end








