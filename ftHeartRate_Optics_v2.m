%%  ftHeartRate_Optics_v2.m

%   Written by Dave Haas between 9 July and 14 August 2021
%
%   Last updated by Dave Haas on 16 August 2021, to make wsstridges from
%   findpeaks data.

clc;
clear;

%% Load globals and commonly-used color defintions

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

%% I. Make introductions...

fprintf('Welcome to the FaunaTag optical heart rate analysis app.\n');

makePlots = lower(input('Do you want to see diagnostic plots (y/n) [n]: ','s'));
switch(makePlots)
    case 'y'
        makePlots = true;
    otherwise
        makePlots = false;
end

%% II. select a tag for analysis

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

%% III. Select which of the five ranges to use for analysis

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
                E_PEAKS_ODBA = E1_PEAKS_ODBA;
                E_PEAKS_ODAV = E1_PEAKS_ODAV;
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
                E_PEAKS_ODBA = E2_PEAKS_ODBA;
                E_PEAKS_ODAV = E2_PEAKS_ODAV;                
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
                E_PEAKS_ODBA = E3_PEAKS_ODBA;
                E_PEAKS_ODAV = E3_PEAKS_ODAV;                
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
                E_PEAKS_ODBA = E4_PEAKS_ODBA;
                E_PEAKS_ODAV = E4_PEAKS_ODAV;
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
                E_PEAKS_ODBA = E5_PEAKS_ODBA;
                E_PEAKS_ODAV = E5_PEAKS_ODAV;                
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


%% IV. Display epoch data

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


%% V. Confirm this epoch is the one you want

% load this EPOCH's optical data into led1-4

led1 = oLed1;
led2 = oLed2;
led3 = oLed3;
led4 = oLed4;

% perform the simple reflectance inversion for ease of interpretation
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

% Make sure that sizes of optics and kinematics are roughly 2.5:1 

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


%% VI. Decimate optical and SCG data from the selected epoch

% specify a common sample frequency after downsample / decimation

dFs = 50;

% specify downsample / decimation factors to use for 50 Hz down/dec

o50Hz = 5;
k50Hz = 2;

% build a common 50 Hz time basis

dT = downsample(oT, o50Hz);
dT_s = downsample(oT_s, o50Hz);

% decimate raw optical channel data using a 3rd order FIR

dLed1 = decimate(oLed1, o50Hz, 3, 'FIR');
dLed2 = decimate(oLed2, o50Hz, 3, 'FIR');
dLed3 = decimate(oLed3, o50Hz, 3, 'FIR');
dLed4 = decimate(oLed4, o50Hz, 3, 'FIR');

% finally, decimate key values from E (i.e., E = E1 for EPOCH 1)

dBpOdba = decimate(E.bpOdba, k50Hz);
dBpOdav = decimate(E.bpOdav, k50Hz);
dBestImfBpOdba = decimate(E.bestImfBpOdba, k50Hz);
dBestImfBpOdav = decimate(E.bestImfBpOdav, k50Hz);
dConsensusHr = decimate(E.consensusHr, k50Hz);

if ( ( numel(dBestImfBpOdav) ~= numel(dLed2) ) )
   fprintf('Decimated K and O length mismatch. Fixing...\n'); 
   koDiff = numel(dLed2) - numel(dBestImfBpOdba);
   if (koDiff > 0)
        fprintf('Adding %d elements to k sigs.\n', koDiff);
        dBpOdba(end+koDiff) = dBpOdba(end);
        dBpOdav(end+koDiff) = dBpOdav(end);
        dBestImfBpOdba(end+koDiff) = dBestImfBpOdba(end);
        dBestImfBpOdav(end+koDiff) = dBestImfBpOdav(end);
        dConsensusHr(end+koDiff) = dConsensusHr(end);
   else
        fprintf('Removing %d elements from k sigs.\n', koDiff);
        dBestImfBpOdbaScale = dBestImfBpOdbaScale(1:end-koDiff);
        dBestImfBpOdavScale = dBestImfBpOdavScale(1:end-koDiff);
   end
   
else
   fprintf('Decimated K and O lengths match. Proceeding...\n');   
end

fprintf('Decimated optical channels and key EPOCH SCG data to 50 Hz.\n');


%% VII. Bandpass filter led1-4 into filtLed1-4

switch(epochChoice)
    case 3                  % apnea
        bandPassFiltOrder  = 3;      % 3rd order butterworth filter
        oLowCut            = 0.5;    % ~60 beats per minute
        oHighCut           = 7.5;    % ~90 beats per minute
        % oHighCut           = 60;
        opticalTestFs      = oHighCut * 2;
        %opticalTestFs      = 8;      % 8 produces good pulses on optics apnea        
    otherwise               % pre- or post-apnea (E1, E2, E4, or E5)
        bandPassFiltOrder  = 3;      % 6th order butterworth filter
        oLowCut            = 2.5;    % ~30 beats per minute
        %oHighCut           = 3.5;   % ~210 beats per minute <- 3.5 test 12.5
        oHighCut           = 12.5;   % ~210 beats per minute <- 3.5 test 12.5
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

if (makePlots)
    fvtO = fvtool(sosBpo, bandPassOpticalFilter, 'Fs', opticalTestFs);
    legend(fvtO,'butter','designfilt');
end

OPTICS = struct;
OPTICS.BANDPASS_FILTER.filterOrder = bandPassFiltOrder;
OPTICS.BANDPASS_FILTER.lowCut = oLowCut;
OPTICS.BANDPASS_FILTER.highCut = oHighCut;
OPTICS.BANDPASS_FILTER.sampleRate = opticalTestFs;
OPTICS.BANDPASS_FILTER.A = Aobp;
OPTICS.BANDPASS_FILTER.B = Bobp;
OPTICS.BANDPASS_FILTER.C = Cobp;
OPTICS.BANDPASS_FILTER.D = Dobp;
OPTICS.BANDPASS_FILTER.SOS = sosBpo;
OPTICS.BANDPASS_FILTER.G = gBpo;

fprintf('Band-pass butterworth filter for optical cardio energies:\n');
fprintf('\tFilter order: %d\n', bandPassFiltOrder);
fprintf('\tLow cut frequency: %2.2f Hz\n', oLowCut);
fprintf('\tHigh cut frequency: %2.2f Hz\n', oHighCut);
fprintf('\tOptical BPF sample rate (2Hz = normalized frequencies): %d Hz\n', ...
    opticalTestFs);

% Use this filter to make good optical heart rate signals

filtLed1 = filtfilt(sosBpo, gBpo, oLed1);
filtLed2 = filtfilt(sosBpo, gBpo, oLed2);
filtLed3 = filtfilt(sosBpo, gBpo, oLed3);
filtLed4 = filtfilt(sosBpo, gBpo, oLed4);

dFiltLed1 = filtfilt(sosBpo, gBpo, dLed1);
dFiltLed2 = filtfilt(sosBpo, gBpo, dLed2);
dFiltLed3 = filtfilt(sosBpo, gBpo, dLed3);
dFiltLed4 = filtfilt(sosBpo, gBpo, dLed4);

% Generate some stats for histogram displays in next section (VII)

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


mean_dFiltLed1 = mean(dFiltLed1);
mean_dFiltLed2 = mean(dFiltLed2);
mean_dFiltLed3 = mean(dFiltLed3);
mean_dFiltLed4 = mean(dFiltLed4);

std_dFiltLed1 = std(dFiltLed1);
std_dFiltLed2 = std(dFiltLed2);
std_dFiltLed3 = std(dFiltLed3);
std_dFiltLed4 = std(dFiltLed4);

r95th_dFiltLed1 = mean_dFiltLed1 + std_dFiltLed1 * 1.96;
l95th_dFiltLed1 = mean_dFiltLed1 - std_dFiltLed1 * 1.96;

r95th_dFiltLed2 = mean_dFiltLed2 + std_dFiltLed2 * 1.96;
l95th_dFiltLed2 = mean_dFiltLed2 - std_dFiltLed2 * 1.96;

r95th_dFiltLed3 = mean_dFiltLed3 + std_dFiltLed3 * 1.96;
l95th_dFiltLed3 = mean_dFiltLed3 - std_dFiltLed3 * 1.96;

r95th_dFiltLed4 = mean_dFiltLed4 + std_dFiltLed4 * 1.96;
l95th_dFiltLed4 = mean_dFiltLed4 - std_dFiltLed4 * 1.96;



%% VIII. Make band-passed optical channel plots and histograms

fprintf('Displaying all band-passed optical channels, and displaying.\n');
fprintf('histograms and summary status for band-passed optical channels.\n');


figOptics = figure('Position',[100 100 1600 1200], 'Color', white);
figOptics.Name = 'Band-pass filtered optical channels';

s1a = subplot(421);
plot(dT, dFiltLed1, 'Color', black);
% plot(oT, filtLed1, 'Color', black);
% hold on; 
% plot(oT, dFiltLed1, 'Color', [0.5 0.5 0.5]);
% hold off;
xlabel('Time, local');
ylabel('Reflectance (A.U.)');
title('Ambient optical channel');
grid;

s1b = subplot(423);
plot(dT, dFiltLed2, 'Color', blue);
% plot(oT, filtLed2, 'Color', blue);
xlabel('Time, local');
ylabel('Reflectance (A.U.)');
title('1050 nm optical channel');
grid;

s1c = subplot(425);
plot(dT, dFiltLed3, 'Color', red);
% plot(oT, filtLed3, 'Color', red);
xlabel('Time, local');
ylabel('Reflectance (A.U.)');
title('1200 nm optical channel');
grid;

s1d = subplot(427);
plot(dT, dFiltLed4, 'Color', green);
% plot(oT, filtLed4, 'Color', green);
xlabel('Time, local');
ylabel('amb2 (nonOp 950 nm)');
title('Ambient (inoperative 950 nm) optical channel');  
grid;

fHist1 = subplot(422);
histogram(dFiltLed1, 'FaceColor', [.7 .7 .7], 'EdgeColor', black);
xline(mean_dFiltLed1, 'Color', goldenrod, 'LineWidth', 1.2);
xline(r95th_dFiltLed1, '--', 'Color', goldenrod, 'LineWidth', 1.2);
xline(l95th_dFiltLed1, '--', 'Color', goldenrod, 'LineWidth', 1.2);
title('Channel 1: ambient');
grid;

fHist2 = subplot(424);
histogram(dFiltLed2, 'FaceColor', blue, 'EdgeColor', black);
xline(mean_dFiltLed2, 'Color', goldenrod, 'LineWidth', 1.2);
xline(r95th_dFiltLed2, '--', 'Color', goldenrod, 'LineWidth', 1.2);
xline(l95th_dFiltLed2, '--', 'Color', goldenrod, 'LineWidth', 1.2);
title('Channel 2: 1050 nm');
grid;

fHist3 = subplot(426);
histogram(dFiltLed3, 'FaceColor', red, 'EdgeColor', black);
xline(mean_dFiltLed3, 'Color', goldenrod, 'LineWidth', 1.2);
xline(r95th_dFiltLed3, '--', 'Color', goldenrod, 'LineWidth', 1.2);
xline(l95th_dFiltLed3, '--', 'Color', goldenrod, 'LineWidth', 1.2);
title('Channel 3: 1200 nm');
grid;

fHist4 = subplot(428);
histogram(dFiltLed4, 'FaceColor', green, 'EdgeColor', black);
xline(mean_dFiltLed4, 'Color', goldenrod, 'LineWidth', 1.2);
xline(r95th_dFiltLed4, '--', 'Color', goldenrod, 'LineWidth', 1.2);
xline(l95th_dFiltLed4, '--', 'Color', goldenrod, 'LineWidth', 1.2);
title('Channel 4: ambient, inoperative 950 nm)');
grid;

linkaxes([s1a s1b s1c s1d],'x');


%% IX. Compare how well HR pulses are preserved using filloutliers w/ spline

fprintf('Plotting outlier-replacement with splining for optical channels...\n');

[deOutFiltLed1, deOutIndexLed1, deOutLoThreshLed1, deOutHiThreshLed1 ] = ...
    filloutliers(dFiltLed1, 'spline');

[deOutFiltLed2, deOutIndexLed2, deOutLoThreshLed2, deOutHiThreshLed2 ] = ...
    filloutliers(dFiltLed2, 'spline');

[deOutFiltLed3, deOutIndexLed3, deOutLoThreshLed3, deOutHiThreshLed3 ] = ...
    filloutliers(dFiltLed3, 'spline');

[deOutFiltLed4, deOutIndexLed4, deOutLoThreshLed4, deOutHiThreshLed4 ] = ...
    filloutliers(dFiltLed4, 'spline');


figDeOut = figure('Position',[100 100 1200 1400 ], 'Color', white);
figDeOut.Name = 'filloutliers with splining for optical channels';

sOut1 = subplot(411);
plot(dT, dFiltLed1, '--', 'Color', goldenrod);
hold on;
plot(dT, deOutFiltLed1, 'Color', black);
hold off
xlabel('Time, local');
ylabel('AU');
title('Channel 1: ambient');
grid;
legend('pre-outlier signal','outlier-replaced signal');

sOut2 = subplot(412);
plot(dT, dFiltLed2, '--', 'Color', goldenrod);
hold on;
plot(dT, deOutFiltLed2, 'Color', blue);
hold off
xlabel('Time, local');
ylabel('AU');
title('Channel 2: 1050 nm');
grid;
legend('pre-outlier signal','outlier-replaced signal');

sOut3 = subplot(413);
plot(dT, dFiltLed3, '--', 'Color', goldenrod);
hold on;
plot(dT, deOutFiltLed3, 'Color', red);
hold off
xlabel('Time, local');
ylabel('AU');
title('Channel 3: 1200 nm');
grid;
legend('pre-outlier signal','outlier-replaced signal');

sOut4 = subplot(414);
plot(dT, dFiltLed4, '--', 'Color', goldenrod);
hold on;
plot(dT, deOutFiltLed4, 'Color', green);
hold off
xlabel('Time, local');
ylabel('AU');
title('Channel 4: ambient (inoperable 950 nm)');
grid;
legend('pre-outlier signal','outlier-replaced signal');

linkaxes([sOut1,sOut2,sOut3,sOut4],'x');  


%% X. Just say yes to using outlier-smoothed optical data.

replaceTest = lower(input('Use outlier-smoothed optical data (y/n) [y]: ','s'));

switch(replaceTest)
    case 'n'
        
        outliersSmoothed = false;
        OPTICS.outlierSmoothed = false;
        fprintf('Leaving original band-passed optical data in place.\n');     
        
    otherwise
        
        outliersSmoothed = true;
        OPTICS.outlierSmoothed = true;
        fprintf('Using outlier-smoothed band-passed optical data...\n');
        
        dFiltLed1 = deOutFiltLed1;
        dFiltLed2 = deOutFiltLed2;
        dFiltLed3 = deOutFiltLed3;
        dFiltLed4 = deOutFiltLed4;

        fprintf('Here are histograms after outlier-smoothing...\n');

end

%  take a look at histograms if outlier-smoothing is used

if (outliersSmoothed)

mean_dFiltLed1 = mean(dFiltLed1);
mean_dFiltLed2 = mean(dFiltLed2);
mean_dFiltLed3 = mean(dFiltLed3);
mean_dFiltLed4 = mean(dFiltLed4);

std_dFiltLed1 = std(filtLed1);
std_dFiltLed2 = std(filtLed2);
std_dFiltLed3 = std(filtLed3);
std_dFiltLed4 = std(filtLed4);

r95th_dFiltLed1 = meanFiltLed1 + stdFiltLed1 * 1.96;
l95th_dFiltLed1 = meanFiltLed1 - stdFiltLed1 * 1.96;

r95th_dFiltLed2 = meanFiltLed2 + stdFiltLed2 * 1.96;
l95th_dFiltLed2 = meanFiltLed2 - stdFiltLed2 * 1.96;

r95th_dFiltLed3 = meanFiltLed3 + stdFiltLed3 * 1.96;
l95th_dFiltLed3 = meanFiltLed3 - stdFiltLed3 * 1.96;

r95th_dFiltLed4 = meanFiltLed4 + stdFiltLed4 * 1.96;
l95th_dFiltLed4 = meanFiltLed4 - stdFiltLed4 * 1.96;

figOptics = figure('Position',[100 100 1600 1200], 'Color', white);
figOptics.Name = 'Band-pass filtered optical channels';

s1a = subplot(421);
plot(dT, dFiltLed1, 'Color', black);
xlabel('Time, local');
ylabel('Reflectance (A.U.)');
title('Ambient optical channel');
grid;

s1b = subplot(423);
plot(dT, dFiltLed2, 'Color', blue);
xlabel('Time, local');
ylabel('Reflectance (A.U.)');
title('1050 nm optical channel');
grid;

s1c = subplot(425);
plot(dT, dFiltLed3, 'Color', red);
xlabel('Time, local');
ylabel('Reflectance (A.U.)');
title('1200 nm optical channel');
grid;

s1d = subplot(427);
plot(dT, dFiltLed4, 'Color', green);
xlabel('Time, local');
ylabel('amb2 (nonOp 950 nm)');
title('Ambient (inoperative 950 nm) optical channel');  
grid;

fHist1 = subplot(422);
histogram(dFiltLed1, 'FaceColor', [.7 .7 .7], 'EdgeColor', black);
xline(mean_dFiltLed1, 'Color', goldenrod, 'LineWidth', 1.2);
xline(r95th_dFiltLed1, '--', 'Color', goldenrod, 'LineWidth', 1.2);
xline(l95th_dFiltLed1, '--', 'Color', goldenrod, 'LineWidth', 1.2);
title('Channel 1: ambient');
grid;

fHist2 = subplot(424);
histogram(dFiltLed2, 'FaceColor', blue, 'EdgeColor', black);
xline(mean_dFiltLed2, 'Color', goldenrod, 'LineWidth', 1.2);
xline(r95th_dFiltLed2, '--', 'Color', goldenrod, 'LineWidth', 1.2);
xline(l95th_dFiltLed2, '--', 'Color', goldenrod, 'LineWidth', 1.2);
title('Channel 2: 1050 nm');
grid;

fHist3 = subplot(426);
histogram(dFiltLed3, 'FaceColor', red, 'EdgeColor', black);
xline(mean_dFiltLed3, 'Color', goldenrod, 'LineWidth', 1.2);
xline(r95th_dFiltLed3, '--', 'Color', goldenrod, 'LineWidth', 1.2);
xline(l95th_dFiltLed3, '--', 'Color', goldenrod, 'LineWidth', 1.2);
title('Channel 3: 1200 nm');
grid;

fHist4 = subplot(428);
histogram(dFiltLed4, 'FaceColor', green, 'EdgeColor', black);
xline(mean_dFiltLed4, 'Color', goldenrod, 'LineWidth', 1.2);
xline(r95th_dFiltLed4, '--', 'Color', goldenrod, 'LineWidth', 1.2);
xline(l95th_dFiltLed4, '--', 'Color', goldenrod, 'LineWidth', 1.2);
title('Channel 4: ambient, inoperative 950 nm)');
grid;

linkaxes([s1a s1b s1c s1d],'x');

end

%% XI. DC-removal filter - NOTE: not significantly different from filtLed*

dcoFiltOrder        = 3;                % FieldTrip suggested default
dcoPassbandFreqHz   = 0.25;             % DC-removal (cut below 0.01 Hz)

fprintf('DC-removal (high-pass) filter characteristics:\n');
fprintf('\tFilter order: %d\n', dcoFiltOrder);
fprintf('\tPassband frequency: %d Hz\n', dcoPassbandFreqHz);
fprintf('\tNyquist frequency (optics): %d\n', oNyquistFs);
fprintf('\tNyquist frequency (kinematics): %d\n', kNyquistFs);

% filter coefficients for optics high-pass (DC-removal)
[AoDCHP, BoDCHP, CoDCHP, DoDCHP] = ...
    butter(dcoFiltOrder, ( dcoPassbandFreqHz / (oNyquistFs) ), 'high');

[sosDco,gDco] = ss2sos(AoDCHP, BoDCHP, CoDCHP, DoDCHP);

% design the high-pass filter. NOTE: using samplingRate = 2, which is the
% designfilt default, results in normalized frequencies

dcoOpticalFilter = designfilt('highpassiir','FilterOrder',dcoFiltOrder, ...
    'PassbandFrequency', dcoPassbandFreqHz, ... 
    'PassbandRipple', 0.2, ...
    'SampleRate', 2);

if (makePlots)
    fvt = fvtool(sosDco, dcoOpticalFilter, 'Fs', oFs);
    legend(fvt,'butter','designfilt');
end

dDcoLed2 = filtfilt(sosDco, gDco, dLed2);
dDcoLed3 = filtfilt(sosDco, gDco, dLed3);

dDcoFiltLed1 = filtfilt(sosDco, gDco, dFiltLed1);
dDcoFiltLed2 = filtfilt(sosDco, gDco, dFiltLed2);
dDcoFiltLed3 = filtfilt(sosDco, gDco, dFiltLed3);
dDcoFiltLed4 = filtfilt(sosDco, gDco, dFiltLed4);

%% XII. Gaussian window filter to smooth optical data - not useful yet

% set a 100 ms windowSize to use for gaussian filtering
%optGaussianWindowSize = oFs / 10;      % = 0.040 s
optGaussianWindowSize = oFs / 25;       % = 0.100 s

gfBo = ones(1, optGaussianWindowSize) / optGaussianWindowSize;
gfA = [1];

dGfFiltLed1 = filtfilt(gfBo, gfA, dFiltLed1);
dGfFiltLed2 = filtfilt(gfBo, gfA, dFiltLed2);
dGfFiltLed3 = filtfilt(gfBo, gfA, dFiltLed3);
dGfFiltLed4 = filtfilt(gfBo, gfA, dFiltLed4);

figGaussOpt = figure('Position',[100 500 1500 600], 'Color', white);

pGf1 = subplot(411);
plot(dT, dFiltLed1, 'Color', goldenrod);
hold on;
plot(dT, dGfFiltLed1, 'Color', black);
hold off;
grid;
legend('bpFilt','gf+bp');

pGf2 = subplot(412);
plot(dT, dFiltLed2, 'Color', goldenrod);
hold on;
plot(dT, dGfFiltLed2, 'Color', blue);
hold off;
grid;
legend('bpFilt','gf+bp');

pGf3 = subplot(413);
plot(dT, dFiltLed3, 'Color', goldenrod);
hold on;
plot(dT, dGfFiltLed3, 'Color', red);
hold off;
grid;
legend('bpFilt','gf+bp');

pGf4 = subplot(414);
plot(dT, dFiltLed4, 'Color', goldenrod);
hold on;
plot(dT, dGfFiltLed4, 'Color', green);
hold off;
grid;
legend('bpFilt','gf+bp');

linkaxes([pGf1 pGf2 pGf3 pGf4],'x');


%% XIII. Create first derivatives of optical channels (extra peak-y!)

% NOTE: all diff optical data is padded with a zero for its last sample, 
% due to reduction of the data set by one sample by the diff() command.

% take first derivative of decimated dFiltLed*

dDiffFiltLed1 = diff(dFiltLed1 * -1);
dDiffFiltLed2 = diff(dFiltLed2 * -1);
dDiffFiltLed3 = diff(dFiltLed3 * -1);
dDiffFiltLed4 = diff(dFiltLed4 * -1);

dDiffFiltLed1(end+1) = 0;
dDiffFiltLed2(end+1) = 0;
dDiffFiltLed3(end+1) = 0;
dDiffFiltLed4(end+1) = 0;

% and make the absolute value versions of these too

dAbsLed1 = abs(dLed1);
dAbsLed2 = abs(dLed2);
dAbsLed3 = abs(dLed3);
dAbsLed4 = abs(dLed4);

dAbsDiffFiltLed1 = abs(dDiffFiltLed1);
dAbsDiffFiltLed2 = abs(dDiffFiltLed2);
dAbsDiffFiltLed3 = abs(dDiffFiltLed3);
dAbsDiffFiltLed4 = abs(dDiffFiltLed4);

dAbsDiffFiltActiveOchans = dAbsDiffFiltLed2 + dAbsDiffFiltLed3;
dAbsDiffFiltAmbientOchans = dAbsDiffFiltLed1 + dAbsDiffFiltLed4;

% have a look at the dLed* (decimated filtLed*) vs. diffLed*

figDiff = figure('Position', [100 100 1200 1400], 'Color', white);

pDiff1 = subplot(211);
plot(dT, dFiltLed2, 'Color', goldenrod);
hold on;
plot(dT, dDiffFiltLed2, 'Color', blue);
hold on;
xlabel('Time, local');
ylabel('Reflectance (A.U.)');
title('1050 nm - Band-passed reflectance and its first derivative');
grid;
legend('filtLed2','diffLed2');

pDiff2 = subplot(212);
plot(dT, dFiltLed3, 'Color', goldenrod);
hold on;
plot(dT, dDiffFiltLed3, 'Color', red);
hold on;
xlabel('Time, local');
ylabel('Reflectance (A.U.)');
title('1200 nm - Band-passed reflectance and its first derivative');
grid;
legend('filtLed3','diffLed3');

linkaxes([pDiff1 pDiff2],'x');


%% XIV. Create root sum squared versions of these optical channels

%   NOTE: rss = root sum squared, i.e.: the Euclidian norm of a vector

% start with dFiltLed*

dRssFiltLed1 = sqrt( dFiltLed1 .^ 2);
dRssFiltLed2 = sqrt( dFiltLed2 .^ 2);
dRssFiltLed3 = sqrt( dFiltLed3 .^ 2);
dRssFiltLed4 = sqrt( dFiltLed4 .^ 2);

dRssFiltActiveOchans  = sqrt( dFiltLed2 .^ 2 + dFiltLed3 .^ 2 );
dRssFiltAmbientOchans = sqrt( dFiltLed1 .^ 2 + dFiltLed4 .^ 2 );

% then do dDiffFiltLed*

dRssDiffFiltLed1 = sqrt( dDiffFiltLed1 .^ 2 );
dRssDiffFiltLed2 = sqrt( dDiffFiltLed2 .^ 2 );
dRssDiffFiltLed3 = sqrt( dDiffFiltLed3 .^ 2 );
dRssDiffFiltLed4 = sqrt( dDiffFiltLed4 .^ 2 );

dRssDiffFiltActiveOchans  = sqrt( dDiffFiltLed2 .^ 2 + dDiffFiltLed3 .^ 2 );
dRssDiffFiltAmbientOchans = sqrt( dDiffFiltLed1 .^ 2 + dDiffFiltLed4 .^ 2 );

fprintf('Be patient. Doing an exploratory 20 IMF VMD of dRssFiltActiveOchans.\n');
tic
[dImfRssFiltActiveOchans] = vmd(dRssFiltActiveOchans,'NumIMFs',20);
toc

figure;

s1a = subplot(411);
plot(dT, dRssFiltLed2, 'Color', blue);
hold on;
plot(dT, dRssDiffFiltLed2, 'Color', cyan);
hold off;
xlabel('Time, local');
ylabel('RSS (A.U.)');
title('RSS dFiltLed2 and dDiffFiltLed2');
grid;
legend('dRssFiltLed2','dRssDiffFiltLed2');

s1b = subplot(412);
plot(dT, dRssFiltLed3, 'Color', red);
hold on;
plot(dT, dRssDiffFiltLed3, 'Color', purple);
hold off;
xlabel('Time, local');
ylabel('RSS (A.U.)');
title('RSS dFiltLed3 and dDiffFiltLed3');
grid;
legend('dRssFiltLed3','dRssDiffFiltLed3');

s1c = subplot(413);
plot(dT, dRssFiltActiveOchans, 'Color', maroon);
hold on;
plot(dT, dRssDiffFiltActiveOchans, 'Color', goldenrod);
plot(dT, dImfRssFiltActiveOchans(:,17), 'Color', red);
hold off;
xlabel('Time, local');
ylabel('RSS (A.U.)');
title('RSS Active optical channels (Filt and DiffFilt at 1050 nm and 1200 nm)');
grid;
legend('dRssFiltActive','dRssDiffFiltActive');

s1d = subplot(414);
plot(dT, dRssFiltAmbientOchans, 'Color', black);
hold on;
plot(dT, dRssDiffFiltAmbientOchans, 'Color', [0.5 0.5 0.5] );
hold off;
xlabel('Time, local');
ylabel('RSS (A.U.)');
title('RSS Ambient optical channels (1050 nm and 1200 nm)');
grid; 
legend('dRssFiltAmbient','dRssDiffFiltAmbient');

linkaxes([s1a s1b s1c s1d],'x');

%% XV. Make some correlation plots to ensure this hasn't gone amok

figLed2 = figure('Color', white, 'Name', 'dFiltLed2 : dFiltLed3' ); 
[rValLed2, pValLed2, figLed2] = corrplot( [dFiltLed2, dFiltLed3] ,  ...
    'type', 'Pearson', 'tail', 'both', 'testR', 'on', 'alpha', 0.01, ...
    'varNames', {'dFO2','dFO3'} );

figLed23 = figure('Color', white, 'Name', 'diffFiltLed2:diffFiltLed3');
[rValLed23, pValLed23, figLed23] = corrplot( [dDiffFiltLed2, dDiffFiltLed3] ,  ...
    'type', 'Pearson', 'tail', 'both', 'testR', 'on', 'alpha', 0.01, ...
    'varNames', {'d1FO2','d1FO3'} );

figLed2 = figure('Color', white, 'Name', 'dFiltLed2 : dDiffFiltLed2' ); 
[rValLed2, pValLed2, figLed2] = corrplot( [dFiltLed2, dDiffFiltLed2] ,  ...
    'type', 'Pearson', 'tail', 'both', 'testR', 'on', 'alpha', 0.01, ...
    'varNames', {'O2','d1O2'} );

figLed3 = figure('Color', white, 'Name', 'dFiltLed3 : dDiffFiltLed3' ); 
[rValLed3, pValLed3, figLed3] = corrplot( [dFiltLed3, dDiffFiltLed3] ,  ...
    'type', 'Pearson', 'tail', 'both', 'testR', 'on', 'alpha', 0.01, ...
    'varNames', {'O3','d1O3'} );


%% XVI. Plot wsst and ridge of dImfRssFiltActiveOchans

% NOTE: use this to assess the need for regionalizing

fprintf('Be patient... doing an exploratory WSST and ridge-find...\n');

tic

[wsstTest, fTest] = wsst(dImfRssFiltActiveOchans(:,17), dFs);

[testRidge10, testIridge] = wsstridge(wsstTest, 10, fTest);
[testRidge20, testIridge] = wsstridge(wsstTest, 20, fTest);
[testRidge30, testIridge] = wsstridge(wsstTest, 30, fTest);
[testRidge40, testIridge] = wsstridge(wsstTest, 40, fTest);

meanRidge = mean([ testRidge10' testRidge20' testRidge30' testRidge40' ], 2);

toc

figure;

subplot(211);
pcolor(dT, fTest, abs(wsstTest));
shading interp;
xlabel('Time, local');
ylabel('Frequency, Hz');
title('WSST for IMF 17 of dRssFiltActiveOchans');
ylim([0 4]);

subplot(212);
plot(dT, dConsensusHr, 'Color', red, 'LineWidth', 2);
hold on;
plot(dT, meanRidge, 'Color', blue);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Single wsstridge find from IMF 17 of dRssFiltActiveOchans');
ylim([0 4]);
grid;

figSizeBA = [ 100 100 600 600 ];

fig1 = figure('Position', figSizeBA, 'Color', white, 'NumberTitle', 'off');
[ rpc1, fig1, stats1 ]  = BlandAltman(fig1, ...
    dConsensusHr, meanRidge, ...
    {'consensusSCG','meanRidge'}, ...
    'Bland-Altman test for SCG vs. IMF17 dRssFiltActiveOchans', ...
    {'SCG','dRSS'}, ...
    'corrInfo', {'eq','r','r2','p','SSE','RMSE','n'}, ...
    'baStatsMode','Non-parametric');


%% XV. Regionalize analysis to partition artifacts

regionalizeTxt = 'Do you want to regionalize optical analysis? (y/n) [n]: ';
regionalize = lower(input(regionalizeTxt,'s'));
switch(regionalize)
    case 'y'
        doRegionalAnalysis = true;
    otherwise
        doRegionalAnalysis = false;
end

if (doRegionalAnalysis)         % start doRegionalAnalysis

clc;

if (exist('REGION','var'))
    clear REGION;
end

if (exist('oCue','var'))
    clear oCue;
end

if (exist('regionRange','var'))
    clear regionRange;
end

startCue = 1;
endCue = numel(dLed3);

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
plot(dDiffFiltLed2, 'Color', blue);
hold on;
plot(dRssFiltActiveOchans);
xline((startCue), 'Color', green', 'LineWidth', 2);
xline((endCue), 'Color', red, 'LineWidth', 2);
hold off;
xlabel('Samples, 50 Hz');
ylabel('AU');
grid;

pMark3 = subplot(312);
pMark3.Position(1) = .10;
pMark3.Position(3) = .80;
plot(dDiffFiltLed3, 'Color', red);
hold on;
xline((startCue), 'Color', green', 'LineWidth', 2);
xline((endCue), 'Color', red, 'LineWidth', 2);
hold off;
xlabel('Samples, 50 Hz');
ylabel('AU');
grid;

pMark4 = subplot(313);
pMark4.Position(1) = .10;
pMark4.Position(3) = .80;
plot(dDiffFiltLed4, 'Color', green);
hold on;
plot(dDiffFiltLed1, 'Color', black);
xline((startCue), 'Color', green', 'LineWidth', 2);
xline((endCue), 'Color', red, 'LineWidth', 2);
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
            fprintf('Caught a [s]ave. Creating REGION structure...\n');
            oCue = sort(oCue);
            markedUp = true;
           
            fprintf('Ending mark-up. See REGION.sample for region sample numbers.\n');

            REGION = struct;
            REGION.sample      = round(unique(oCue));
            REGION.epochChoice = epochChoice;

            numRegions = numel(REGION.sample) - 1;

            regionRange(numRegions,2) = zeros;

            for i = 1:numRegions

                regionRange(i,:) = [REGION.sample(i) REGION.sample(i+1) ];

            end

            REGION.ranges = regionRange;
            
            disp(REGION)
           
           
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
    
end         % end while(markedUp == false)         
     
end         % end doRegionalAnalysis


%% XVI. Make a timetable of relevant decimated data for later use

Time = dT;

DTT = timetable(Time, dLed1, dLed2, dLed3, dLed4, ...
    dFiltLed1, dFiltLed2, dFiltLed3, dFiltLed4, ...
    dDiffFiltLed1, dDiffFiltLed2, dDiffFiltLed4, dDiffFiltLed4, ...
    dRssFiltLed1, dRssFiltLed2, dRssFiltLed3, dRssFiltLed4, ...
    dRssDiffFiltLed1, dRssDiffFiltLed2, dRssDiffFiltLed3, dRssDiffFiltLed4, ...
    dRssFiltActiveOchans, dRssFiltAmbientOchans, ...
    dRssDiffFiltActiveOchans, dRssDiffFiltAmbientOchans, ...
    dBpOdba, dBpOdav, dBestImfBpOdba, dBestImfBpOdav, dConsensusHr);

fprintf('Created DTT timetable of common 50 Hz optical and SCG data.\n');

%% XVII. Everything is in place for heart rate analysis... proceed?

%  Pause here and let the user select which analysis to run

fprintf('\n\n');
fprintf('Signal processing is complete on this epoch.\n');
inputTxt = 'Press [enter] to continue analysis, or [q] to go manual: ';
waitForInput = lower(input(inputTxt,'s'));
switch(waitForInput)
    case 'q'
        clc;
        fprintf('Bailing ftHeartRate_Optics.m around line 1400...\n');
        fprintf('Pick it up from there!\n');
        return;
	otherwise
        fprintf('Continuing with analysis...\n');
end


%% XVIII. Create numImf (usually 10) VMD IMFs of active absDiffFilt optics

% get 10 IMFs of VMDs for absDiffFiltLed*

numImfs = 20;
numRows = numImfs / 2;
numCols = numImfs / numRows;

fprintf('Computing VMDs for absDiffFiltLed(2&3) and dAbsDiffFiltLed(2&3)...\n');

tic

fprintf('\t... first for decimated rssFiltActiveOchans...\n');

[ imfRssFiltActive, residualsRssFiltActive, infoRssFiltActive ] = ...
    vmd( dRssFiltActiveOchans, 'NumIMFs', numImfs);

fprintf('\t... and now for decimated rssDiffActiveFiltOchans...\n');

[ imfRssDiffFiltActive, residualsRssDiffFiltActive, infoRssDiffFiltActive] = ...
    vmd( dRssDiffFiltActiveOchans, 'NumIMFs', numImfs);

toc

fprintf('Completed compute of VMDs for undecimated and decimated absDiffFiltLed(2&3).\n');

%% XIX. Make some plots of these IMFs...

figRssFiltActiveImfs = figure('Position', [50 50 1400 800], 'Color', white);

numRows = 5;
numCols = 2;

t = tiledlayout(numRows,numCols,'TileSpacing','compact','Padding','compact');
for n = 11:numImfs
    ax(n) = nexttile(t);
    plot(dT,imfRssFiltActive(:,n)')
    xlim([dT(1) dT(end)])
    txt = ['IMF',num2str(n)];
    title(txt)
    xlabel('Time, local')
end
title(t,'Variational Mode Decomposition - dRssFiltActiveOchans');


figRssDiffFiltActiveImfs = figure('Position', [50 500 1400 800], 'Color', white);

t = tiledlayout(numRows,numCols,'TileSpacing','compact','Padding','compact');
for n = 11:numImfs
    ax(n) = nexttile(t);
    plot(dT,imfRssDiffFiltActive(:,n)');
    xlim([dT(1) dT(end)]);
    txt = ['IMF',num2str(n)];
    title(txt)
    xlabel('Time, local')
end
title(t,'Variational Mode Decomposition - dRssDiffFiltActiveOchans)');



%% XX. Create sumImf and bestImf versions of these

sumImfRssFiltActive = sum(imfRssFiltActive(:,2:19),2);
bestImfRssFiltActive = imfRssFiltActive(:,17);

sumImfRssDiffFiltActive = sum(imfRssDiffFiltActive(:,2:19),2);
bestImfRssDiffFiltActive = imfRssDiffFiltActive(:,19);

   
figSumAndBestImf = figure('Color', white);

p5a = subplot(211);
plot(dT, sumImfRssFiltActive, 'Color', blue);
hold on;
plot(dT, bestImfRssFiltActive, 'Color', green, 'LineWidth', 1.2);
hold off;
xlabel('Time, local');
ylabel('A.U.');
title('sumImf and bestImf - rssFiltActiveOchans');
grid;
legend('sumImf','bestImf');

p5b = subplot(212);
plot(dT, sumImfRssDiffFiltActive, 'Color', red);
hold on;
plot(dT, bestImfRssDiffFiltActive, 'Color', goldenrod, 'LineWidth', 1.2);
hold off;
xlabel('Time, local');
ylabel('A.U.');
title('sumImf and bestImf - rssDiffFiltActiveOchans');
grid;
legend('sumImf','bestImf');

linkaxes([p5a p5b],'x');


%% XXI. Create threshold-based square waves

% build target square waves for bestImfs & sumImfs of rssFilt & rssDiffFilt

% bestImf rssFiltActiveOchans

fprintf('Building target square wave for bestImfRssFiltActiveOchans...\n');

meanBestImfRssFiltActive = mean(bestImfRssFiltActive);
stdBestImfRssFiltActive  = std(bestImfRssFiltActive);

fprintf('\tMean bestImf RssFiltActiveOchan: %.1f\n', ...
    meanBestImfRssFiltActive);
fprintf('\tSD bestImf RssFiltActiveOchan: %.1f\n', ...
    stdBestImfRssFiltActive);

threshBestImfRssFiltActive = meanBestImfRssFiltActive;

sqwThreshBestImfRssFiltActive = double( ...
    (bestImfRssFiltActive <= threshBestImfRssFiltActive) );

% bestImf rssDiffFiltActiveOchans

fprintf('Building target square wave for bestImfRssDiffFiltActiveOchans...\n');

meanBestImfRssDiffFiltActive = mean(bestImfRssDiffFiltActive);
stdBestImfRssDiffFiltActive  = std(bestImfRssDiffFiltActive);

fprintf('\tMean bestImf RssDiffFiltActiveOchan: %.1f\n', ...
    meanBestImfRssDiffFiltActive);
fprintf('\tSD bestImf RssDiffFiltActiveOchan: %.1f\n', ...
    stdBestImfRssDiffFiltActive);

threshBestImfRssDiffFiltActive = meanBestImfRssDiffFiltActive;

sqwThreshBestImfRssDiffFiltActive = double( ...
    (bestImfRssDiffFiltActive <= threshBestImfRssDiffFiltActive) );

% sumImf rssFiltActiveOchans

fprintf('Building target square wave for sumImfRssFiltActiveOchans...\n');

meanSumImfRssFiltActive = mean(sumImfRssFiltActive);
stdSumImfRssFiltActive  = std(sumImfRssFiltActive);

fprintf('\tMean sumImf RssFiltActiveOchan: %.1f\n', ...
    meanSumImfRssFiltActive);
fprintf('\tSD sumImf RssFiltActiveOchan: %.1f\n', ...
    stdSumImfRssFiltActive);

threshSumImfRssFiltActive = meanSumImfRssFiltActive;

sqwThreshSumImfRssFiltActive = double( ...
    (sumImfRssFiltActive <= threshSumImfRssFiltActive) );


% sumImf rssDiffFiltActiveOchans

fprintf('Building target square wave for sumImfRssDiffFiltActiveOchans...\n');

meanSumImfRssDiffFiltActive = mean(sumImfRssDiffFiltActive);
stdSumImfRssDiffFiltActive  = std(sumImfRssDiffFiltActive);

fprintf('\tMean sumImf RssDiffFiltActiveOchan: %.1f\n', ...
    meanSumImfRssDiffFiltActive);
fprintf('\tSD sumImf RssDiffFiltActiveOchan: %.1f\n', ...
    stdSumImfRssDiffFiltActive);

threshSumImfRssDiffFiltActive = meanSumImfRssDiffFiltActive;

sqwThreshSumImfRssDiffFiltActive = double( ...
    (sumImfRssDiffFiltActive <= threshSumImfRssDiffFiltActive) );


% % compute proposedHr and proposedAntiHr threshSqwave signals
% 
% fprintf('Creating threshold-based bestImf Ochan2 and Ochan3 square waves...\n');
% threshSqwaveOchan2 = double(~bwareaopen(targetBestImfOchan2, (dFs/2)));
% 
% fprintf('Creating best IMF bpOdav-based AutC square wave...\n');
% threshSqwaveOchan3 = double(~bwareaopen(targetBestImfOchan3, (dFs/5)));

% plot proposed Hr square waves

figSqwaves = figure('Position',[100 300 1400 800]);

pSq1 = subplot(221);
plot(dT, sqwThreshBestImfRssFiltActive, 'Color', blue);
xlabel('Time, local');
title('Square wave of bestIMF RssFiltActive optical channels');
grid;

pSq2 = subplot(222);
plot(dT, sqwThreshBestImfRssDiffFiltActive, 'Color', red);
xlabel('Time, local');
title('Square wave of bestIMF RssDiffFiltActive optical channels');
grid;

pSq3 = subplot(223);
plot(dT, sqwThreshSumImfRssFiltActive, 'Color', green);
xlabel('Time, local');
title('Square wave of sumIMF RssFiltActive optical channels');
grid;

pSq4 = subplot(224);
plot(dT, sqwThreshSumImfRssDiffFiltActive, 'Color', purple);
xlabel('Time, local');
title('Square wave of sumIMF RssDiffFiltActive optical channels');
grid;

linkaxes([ pSq1 pSq2 pSq3 pSq4 ],'x');


%% XXIII. Do peak findings

if (epochChoice == 3)
    minPeakDistance = 50; 
else
    minPeakDistance = 50;
end


[pkBestImfRssFiltActive, locBestImfRssFiltActive, widBestImfRssFiltActive, ...
    promBestImfRssFiltActive] = findpeaks(bestImfRssFiltActive, ...
    'MinPeakHeight', meanBestImfRssFiltActive, ...
    'MinPeakDistance', minPeakDistance);

[pkBestImfRssDiffFiltActive, locBestImfRssDiffFiltActive, ...
    widBestImfRssDiffFiltActive, promBestImfRssDiffFiltActive] = ...
    findpeaks(bestImfRssDiffFiltActive, 'MinPeakHeight', ...
    meanBestImfRssDiffFiltActive,  ...
    'MinPeakDistance', minPeakDistance);
    
% create timetables from these findpeak results

HR_RssFiltActive = timetable(dT(locBestImfRssFiltActive), ...
    pkBestImfRssFiltActive, widBestImfRssFiltActive, ...
    promBestImfRssFiltActive, 'VariableNames', ...
    {'peaks','width','prominence'});

HR_RssDiffFiltActive = timetable(dT(locBestImfRssDiffFiltActive), ...
    pkBestImfRssDiffFiltActive, widBestImfRssDiffFiltActive, ... 
    promBestImfRssDiffFiltActive, 'VariableNames', ...
    {'peaks','width','prominence'});

% now spline interpolate in context of dT

bestImfRssFiltActiveSpline = interp1(HR_RssFiltActive.Time, ...
    HR_RssFiltActive.peaks, dT, 'spline');
bestImfRssFiltActiveRespCycle = csaps(dT_s(locBestImfRssDiffFiltActive), ...
    pkBestImfRssDiffFiltActive, 0.5, dT_s);

bestImfRssDiffFiltActiveSpline = interp1(HR_RssDiffFiltActive.Time, ...
    HR_RssDiffFiltActive.peaks, dT, 'spline');
bestImfRssDiffFiltActiveRespCycle = csaps(dT_s(locBestImfRssFiltActive), ...
    pkBestImfRssFiltActive, 0.5, dT_s);

cueTimes = CUE.Time(CUE.Time < dT(end) & CUE.Time > dT(1));
cueTypes = CUE.type(unique(cueTimes));

f1 = figure;
f1.Name = 'Respiration test based on bestImfOchan2 and bestImfOchan3';
figSize = f1.Position;
figSize(3) = screenSize(3);
figSize(4) = screenSize(4) * 0.5;
f1.Position = figSize;  

pPeaksRssFiltActive = subplot(211);
plot(dT, bestImfRssFiltActive, 'Color', blue, 'LineWidth', 1.1);
hold on;
plot(dT(locBestImfRssFiltActive), pkBestImfRssFiltActive, 'd', ...
    'Color', blue, 'MarkerFace', blue);
plot(dT, bestImfRssFiltActiveSpline, 'Color', blue);
plot(dT, bestImfRssFiltActiveRespCycle, 'Color', maroon, 'LineWidth', 1.1);
yPos = pPeaksRssFiltActive.YLim(2) - (pPeaksRssFiltActive.YLim(2) * 0.1);
numCuePlots = ftPlotCues(cueTimes, cueTypes, yPos);
hold off;
xlabel('Time, local');
ylabel('Cycle');
title('Best IMF of rssFiltActive optical channels');
grid;
legend('bestImf','Peaks','peakSpline','respCycle');

pPeaksRssDiffFiltActive = subplot(212);
plot(dT, bestImfRssDiffFiltActive, 'Color', red, 'LineWidth', 1.1);
hold on;
plot(dT(locBestImfRssDiffFiltActive), pkBestImfRssDiffFiltActive, 'd', ...
    'Color', blue, 'MarkerFace', blue);
plot(dT, bestImfRssDiffFiltActiveSpline, 'Color', green);
plot(dT, bestImfRssDiffFiltActiveRespCycle, 'Color', maroon, 'LineWidth', 1.1);
yPos = pPeaksRssDiffFiltActive.YLim(2) - (pPeaksRssDiffFiltActive.YLim(2) * 0.1);
numCuePlots = ftPlotCues(cueTimes, cueTypes, yPos);
hold off;
xlabel('Time, local');
ylabel('Cycle');
title('Best IMF of rssDiffFiltActive optical channels');
grid;
legend('bestImf','Peaks','peakSpline','respCycle');

linkaxes([pPeaksRssFiltActive pPeaksRssDiffFiltActive],'x');

%% XXIV. Try to make some frequency conversions of these peaks

bestImfRssFiltActiveLength = numel(bestImfRssFiltActive);
bestImfRssDiffFiltActiveLength = numel(bestImfRssDiffFiltActive);

peakBestImfRssFiltActiveElements = numel(locBestImfRssFiltActive);
peakBestImfRssDiffFiltActiveElements = numel(locBestImfRssDiffFiltActive);

sqwPeakBestImfRssFiltActive = zeros(bestImfRssFiltActiveLength,1);
sqwPeakBestImfRssDiffFiltActive = zeros(bestImfRssDiffFiltActiveLength,1);

for i = 1:peakBestImfRssDiffFiltActiveElements

    thisStart = locBestImfRssDiffFiltActive(i);
    thisEnd = locBestImfRssDiffFiltActive(i) + round( widBestImfRssDiffFiltActive(i));
    if (thisEnd > bestImfRssFiltActiveLength)
        sqwPeakBestImfRssFiltActive(thisStart:bestImfRssFiltActiveLength,1) = 1;
    else
        sqwPeakBestImfRssFiltActive(thisStart:thisEnd,1) = 1;
    end
    
end

for i = 1:peakBestImfRssFiltActiveElements
   
    thisStart = locBestImfRssFiltActive(i);
    thisEnd = locBestImfRssFiltActive(i) + round( widBestImfRssFiltActive(i));
    if (thisEnd > bestImfRssDiffFiltActiveLength)
        sqwPeakBestImfRssDiffFiltActive(thisStart:bestImfRssDiffFiltActiveLength,1) = 1;    
    else
        sqwPeakBestImfRssDiffFiltActive(thisStart:thisEnd,1) = 1;    
    end
    
end


figPeakFreq = figure('Position',[100 500 1400 1000]);

pPeak1 = subplot(211);
plot(dT, sqwPeakBestImfRssFiltActive, 'Color', red);
hold off;
xlabel('Time, seconds');
ylabel('HR Pulses');
title('BestImfOchan2 - binary peak heart rate estimate');

pPeak2 = subplot(212);
plot(dT, sqwPeakBestImfRssDiffFiltActive, 'Color', green);
hold off;
xlabel('Time, seconds');
ylabel('HR pulses');
title('BestImfBpOchan3 - binary peak heart rate estimate');

linkaxes([pPeak1 pPeak2],'x');


%% XXV. Don't forget about this: engine for looping through regional analysis 

if (doRegionalAnalysis)
    
for i = 1:numRegions
    
    this = REGION.ranges(i,:);
    thisT = timerange(dT(this(1)), dT(this(2)) );
        
    numOsamples = numel(dRssFiltActiveOchans(this(1):this(2)) );
    numSamples = this(2) - this(1) + 1;
    
    fprintf('Processing region %d with %d samples...\n', i, numSamples);
    
    % good wsst and wsstridge extraction for leds1-4 during breath-hold on
    % tt21_142c

    [wsstThis, fThis] = wsst(abs(dRssFiltActiveOchans(this(1):this(2))), dFs); 
    [ridge, iridge] = wsstridge(wsstThis, 40, fThis, 'NumRidges', 2);

    figure;
    subplot(211);
    pcolor(dT(this(1):this(2)), fThis, abs(wsstThis) );
    xlabel('Time, seconds');
    ylabel('Frequency, Hz');
    title('abs(dDiffFiltLed2)');
    shading interp;
    ylim([0 5]);

    subplot(212);
    plot(dT(this(1):this(2)), ridge);
    xlabel('Time, seconds');
    ylabel('Frequency, Hz');
    title('abs(dDiffFiltLed2)');
    grid;


    [wsstThis, fThis] = wsst(abs(dRssDiffFiltActiveOchans(this(1):this(2))), dFs); 
    [ridge, iridge] = wsstridge(wsstThis, 40, fThis, 'NumRidges', 2);
    figure;
    subplot(211);
    pcolor(dT(this(1):this(2)), fThis, abs(wsstThis) );
    xlabel('Time, seconds');
    ylabel('Frequency, Hz');
    title('abs(dDiffFiltLed3)');
    shading interp;
    ylim([0 5]);

    subplot(212);
    plot(dT(this(1):this(2)), ridge);
    xlabel('Time, seconds');
    ylabel('Frequency, Hz');
    title('abs(dDiffFiltLed3)');
    grid;    
    
end    

end             % end  doRegionalAnalysis

%% XXVI. WSST and ridge-finding against those WSSTs

fprintf('Starting synchrosqueeze transforms...\n');

tic

% build WSSTs for:
%   1. dRssFiltActiveOchans
%   2. sumImf of dRssFiltActiveOchans
%   3. bestImf of dRssFiltActiveOchans
%   4. sqwThresh of bestImf of dRssFiltActiveOchans
%   5. sqwPeak of bestImf of dRssDiffFiltActiveOchans

fprintf('working on WSSTs for...\n');
fprintf('\t1. dRssFiltActiveOchans\n');
fprintf('\t2. sumImf of dRssFiltActiveOchans\n');
fprintf('\t3. bestImf of dRssFiltActiveOchans\n');
fprintf('\t4. sqwThresh of bestImf of dRssFiltActiveOchans\n');
fprintf('\t5. sqwPeak of bestImf of dRssDiffFiltActiveOchans\n');

[ wsstRssFiltActive, fRssFiltActive ] = ...
    wsst( dRssDiffFiltActiveOchans, dFs);

[ wsstSumImfRssFiltActive, fSumImfRssFiltActive ] = ...
    wsst( sumImfRssFiltActive, dFs );

[ wsstBestImfRssFiltActive, fBestImfRssFiltActive ] = ...
    wsst( bestImfRssFiltActive, dFs );

[ wsstSqwThreshBestImfFiltActive, fSqwThreshBestImfFiltActive ] = ...
    wsst( sqwThreshBestImfRssFiltActive, dFs );

[ wsstSqwPeakBestImfDiffFiltActive, fSqwPeakBestImfDiffFiltActive ] = ...
    wsst( sqwPeakBestImfRssDiffFiltActive, dFs );
    
toc

fprintf('Completed wavelet synchrosqueeze transforms for optical HR.\n');


%% XXVI. Look at some plots of these WSSTs...

figWsstLed2 = figure('Position',[100 100 900 1600], 'Color', white);
figWsstLed.Name = 'Wavelet synchrosqueeze transforms for optical HR estimates';

subplot(511);
pcolor(dT, fRssFiltActive, abs(wsstRssFiltActive));
shading interp;
xlabel('Time, local');
ylabel('Frequency, Hz');
title('RSS band-pass filtered active optical channels');
ylim([0 4]);

subplot(512);
pcolor(dT, fSumImfRssFiltActive, abs(wsstSumImfRssFiltActive ));
shading interp;
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Sum of IMFs 2:19 of RSS band-pass filtered active optical channels');
ylim([0 4]);

subplot(513);
pcolor(dT, fBestImfRssFiltActive, abs(wsstBestImfRssFiltActive ));
shading interp;
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Best IMF of RSS band-passed filtered active optical channels');
ylim([0 4]);

subplot(514);
pcolor(dT, fSqwThreshBestImfFiltActive, abs(wsstSqwThreshBestImfFiltActive ));
shading interp;
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Threshold-based square wave approximation of bestImfRssFiltActive');
ylim([0 4]);

subplot(515);
pcolor(dT, fSqwPeakBestImfDiffFiltActive, abs(wsstSqwPeakBestImfDiffFiltActive ));
shading interp;
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Peak-based square wave approximation of bestImfRssDiffFiltActive');
ylim([0 4]);



%% XXVII. Ensure perfect correlation among fSum, fBest, fThresh & fPeaks

fCorr = corr ( [ fRssFiltActive', fSumImfRssFiltActive', ...
    fBestImfRssFiltActive', fSqwThreshBestImfFiltActive' , ...
    fSqwPeakBestImfDiffFiltActive' ] );

if (fCorr < 1.0001 & fCorr > 0.9999) 
    
    useWholeSpectrum = false;
    fprintf('There is perfect correlation among the wsst fOutputs.\n'); 
    fprintf('Winnowing the fOutput fields to select for 0.5 - 3.5 Hz...\n');
    
    fLimit = find(fBestImfRssFiltActive < 3.5 & ...
        fBestImfRssFiltActive > 0.5);
    
    fStart = fLimit(1);
    fEnd   = fLimit(end);
    
    fprintf('\tDefining a new fRange.\n');
    
    fRange = fBestImfRssFiltActive(1, fStart:fEnd);

    wsstRssFiltActive = wsstRssFiltActive(fStart:fEnd, :);
    wsstSumImfRssFiltActive = wsstSumImfRssFiltActive(fStart:fEnd, :);
    wsstBestImfRssFiltActive = wsstBestImfRssFiltActive(fStart:fEnd, :);
    wsstSqwThreshBestImfFiltActive = ...
        wsstSqwThreshBestImfFiltActive(fStart:fEnd,:);
    wsstSqwPeakBestImfDiffFiltActive = ...
        wsstSqwPeakBestImfDiffFiltActive(fStart:fEnd,:);

    fprintf('Min frequency for ridging: %d -> %2.2f Hz (%3.1f BPM).\n', ...
        fStart, fBestImfRssFiltActive(fLimit(1)), ...
        fBestImfRssFiltActive(fLimit(1)) * 60 );
    fprintf('Max frequency for ridging: %d -> %2.2f Hz (%3.1f BPM)\n', ...
        fEnd, fBestImfRssFiltActive(fLimit(end)), ...
        fBestImfRssFiltActive(fLimit(end)) * 60 );
    
else 
    
    useWholeSpectrum = true;
	fprintf('There is not perfect correlation among the wsst fOutputs.\n'); 
	fprintf('Using complete frequency space for wsstridging...\n');
    
end


%% XXVIII. find ridges for sumImf, bestImf, threshSqave & peaksSqwave WSSTs

fprintf('Finding ridges for these WSSTs... \n');
fprintf('... this might take a minute if using whole spectrum!\n');

tic

numRidges = 2;
penalty   = 20;

fprintf('Finding %d ridges with penality %d for wsstRssFiltActive...\n', ...
    numRidges, penalty);

[ridgeRssFiltActive, iRidgeRssFiltActive] = wsstridge(wsstRssFiltActive, ...
    penalty, fRange, 'NumRidges', numRidges);

% ---===---

fprintf('Finding %d ridges with penality %d for wsstSumImfRssFiltActive...\n', ...
    numRidges, penalty);

[ridgeSumImfRssFiltActive, iRidgeSumImfRssFiltActive] = ...
    wsstridge(wsstSumImfRssFiltActive, ...
    penalty, fRange, 'NumRidges', numRidges);

% ---===---

fprintf('Finding 1 ridge with penalty %d for wsstBestImfRssFiltActive...\n', ...
    20);

[ridgeBestImfRssFiltActive, iRidgeBestImfRssFiltActive] = ...
    wsstridge(wsstBestImfRssFiltActive, ...
    20, fRange, 'NumRidges', 1);

% ---===---

fprintf('Finding 1 ridge with penalty %d for wsstSqwThreshBestImfRssFiltActive...\n', ...
    20);

[ridgeSqwThreshBestImfRssFiltActive, iRidgeSqwThreshBestImfRssFiltActive] = ...
    wsstridge(wsstSqwThreshBestImfFiltActive, ...
    20, fRange, 'NumRidges', 1);


% ---===---

fprintf('Finding 1 ridge with penality %d for wsstSqwPeakBestImfRssDiffFiltActive...\n', ...
    20);

[ridgeSqwPeakBestImfRssDiffFiltActive, iRidgeSqwPeakBestImfRssDiffFiltActive] = ...
    wsstridge(wsstSqwThreshBestImfFiltActive, ...
    20, fRange, 'NumRidges', 1);

toc

fprintf('Completed wsstridge-finding for optical signals of interest.\n');



%% XXIX. Make some WSST ridge plots

figRidges = figure('Position',[100 600 screenSize(3)*0.5 screenSize(4) ], ...
    'Color', white);

pRidge1 = subplot(511);
plot(dT, ridgeRssFiltActive);
hold on; 
plot(dT, dConsensusHr, 'Color', goldenrod, 'LineWidth', 1.5);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Root-sum-squared band-passed active optics');
grid;
legend('r1','r2','SCG');

pRidge2 = subplot(512);
plot(dT, ridgeSumImfRssFiltActive);
hold on; 
plot(dT, dConsensusHr, 'Color', goldenrod, 'LineWidth', 1.1);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Sum IMF of root-sum-squared band-passed active optics');
grid;
legend('r1','r2','SCG');

pRidge3 = subplot(513);
plot(dT, ridgeBestImfRssFiltActive);
hold on; 
plot(dT, dConsensusHr, 'Color', goldenrod, 'LineWidth', 1.1);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Best IMF of root-sum-squared band-passed active optics');
grid;
legend('r1','SCG');

pRidge4 = subplot(514);
plot(dT, ridgeSqwThreshBestImfRssFiltActive);
hold on; 
plot(dT, dConsensusHr, 'Color', goldenrod, 'LineWidth', 1.1);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Threshold-based square wave: best IMF RSS band-passed active optics');
grid;
legend('r1','SCG');

pRidge5 = subplot(515);
plot(dT, ridgeSqwPeakBestImfRssDiffFiltActive);
hold on; 
plot(dT, dConsensusHr, 'Color', goldenrod, 'LineWidth', 1.1);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Peak-based square wave of RSS 1° derivative band-passed active optics');
grid;
legend('r1','SCG');

linkaxes([pRidge1 pRidge2 pRidge3 pRidge4 pRidge5], 'x');

%% XXX. Take a look at a correlation plot for these...


figCorrR1 = figure('Position', [100 100 600 600], 'Color', white, ...
    'Name', 'Consensus SCG HR Frequency vs. Ridge 1 Ochan2');

[rR1, pValR1, figCorrR1 ] = corrplot( ...
    [ dConsensusHr, ridgeRssFiltActive(:,1), ...
     ridgeSumImfRssFiltActive(:,1) ], ...
    'type', 'Pearson', 'alpha', 0.01, 'testR', 'on', 'varName', ...
    {'SCG','RFA','SUMRFA'} );


% ---===---

figCorrR2 = figure('Position', [100 100 600 600], 'Color', white, ...
    'Name', 'Consensus SCG HR Frequency vs. Ridge 1 Ochan2');

[rR2, pValR2, figCorrR2 ] = corrplot( ...
    [ dConsensusHr, ridgeRssFiltActive(:,2), ...
     ridgeSumImfRssFiltActive(:,2) ], ...
    'type', 'Pearson', 'alpha', 0.01, 'testR', 'on', 'varName', ...
    {'SCG','RFA','SUMRFA'} );

% ---===---

figCorrRidges = figure('Position', [100 100 600 600], 'Color', white, ...
    'Name', 'Consensus SCG HR Frequency vs. Useful Ridges');

[rRidges, pValRidges, figCorrRidges ] = corrplot( ...
    [ dConsensusHr, ridgeRssFiltActive(:,1), ridgeBestImfRssFiltActive', ...
    ridgeSqwThreshBestImfRssFiltActive', ridgeSqwPeakBestImfRssDiffFiltActive' ], ...
    'type', 'Pearson', 'alpha', 0.01, 'testR', 'on', 'varName', ...
    {'SCG','RFA','BEST','THRSH','PEAKS'} );

% ---===---


%% XXXI. Try some regressions for all instances of ridge1

fprintf('-------------------------------------------------------\n\n');
fprintf('-------------------------------------------------------\n\n');

fprintf('Creating a stepwise linear model...\n');
fprintf('\tr1 = ridgeRssFiltActive(:,1)\n');
fprintf('\tr2 = ridgeBestImfRssFiltActive(1,:)\n');
fprintf('\tr3 = ridgeSqwThreshBestImfRssFiltActive(1,:)\n');
fprintf('\tr4 = ridgeSqwPeakBestImfRssDiffFiltActive(1,:)\n');
fprintf('\ty  = dConsensusHr\n');

y = dConsensusHr;
r1 = ridgeRssFiltActive(:,1);
r2 = ridgeBestImfRssFiltActive(1,:)';
r3 = ridgeSqwThreshBestImfRssFiltActive(1,:)';
r4 = ridgeSqwPeakBestImfRssDiffFiltActive(1,:)';

modelOCG = stepwiselm([ r4, r3, r2, r1 ], y , ...
    'Criterion', 'aic', 'VarNames',{'r4','r3','r2','r1','scg'} );

disp(modelOCG);

%%

runRidge1 = lower(input('Run ridge1 Ochan2 stepwise GLM (y/n) [y]: ','s'));
switch(runRidge1)
    case 'n'
        fprintf('Okay, skipping ridge1 Ochan2 GLM.\n');
    otherwise

fprintf('Running regressions for ridge1 Ochan2...\n');
r1_Ochan2_GLM = stepwiselm([ ridgeSumImfOchan2(:,1), ridgeBestImfOchan2(:,1), ...
    ridgeThreshOchan2(:,1), ridgePeaksOchan2(:,1) ], dConsensusHr, ...
    'Criterion', 'aic', 'VarNames',{'sum','best','thresh','peaks','scg'} );
disp(r1_Ochan2_GLM);
fprintf('-------------------------------------------------------\n\n');

end

runRidge1 = lower(input('Run ridge1 Ochan3 stepwise GLM (y/n) [y]: ','s'));
switch(runRidge1)
    case 'n'
        fprintf('Okay, skipping ridge1 Ochan3 GLM.\n');
    otherwise

fprintf('Running regressions for ridge1 Ochan3...\n');
r1_Ochan3_GLM = stepwiselm([ ridgeSumImfOchan3(:,1), ridgeBestImfOchan3(:,1), ...
    ridgeThreshOchan3(:,1), ridgePeaksOchan3(:,1) ], dConsensusHr, ...
    'Criterion', 'aic', 'VarNames',{'sum','best','thresh','peaks','scg'} );
disp(r1_Ochan3_GLM);
fprintf('-------------------------------------------------------\n\n');

end

runRidge2 = lower(input('Run ridge2 Ochan2 stepwise GLM (y/n) [n]: ','s'));
switch(runRidge2)
    case 'y'

fprintf('Running regressions for ridge2 Ochan3...\n');
r2_Ochan2_GLM = stepwiselm([ ridgeSumImfOchan2(:,2), ridgeBestImfOchan2(:,2), ...
    ridgeThreshOchan2(:,2), ridgePeaksOchan2(:,2) ], dConsensusHr, ...
    'Criterion', 'aic', 'VarNames',{'sum','best','thresh','peaks','scg'} );
disp(r2_Ochan2_GLM);
fprintf('-------------------------------------------------------\n\n');

    otherwise
        fprintf('Okay, skipping ridge2 Ochan3 GLM.\n');
    
end

runRidge2 = lower(input('Run ridge2 Ochan3 stepwise GLM (y/n) [n]: ','s'));
switch(runRidge2)
    case 'y'
        
fprintf('Running regressions for ridge2 Ochan2...\n');
r2_Ochan3_GLM = stepwiselm([ ridgeSumImfOchan3(:,2), ridgeBestImfOchan3(:,2), ...
    ridgeThreshOchan3(:,2), ridgePeaksOchan3(:,2) ], dConsensusHr, ...
    'Criterion', 'aic', 'VarNames',{'sum','best','thresh','peaks','scg'} );
disp(r2_Ochan3_GLM);
fprintf('-------------------------------------------------------\n\n');

    otherwise
        fprintf('Okay, skipping ridge2 Ochan3 GLM.\n');
    
end

%% THIS IS END OF NEW / START OF OLD. BRANCH TO CONTINUE WITH OLD IF NEEDED



holdUp = lower(input('You have reached the OLD stuff. Continue(y/n) [n]: ','s'));
switch(holdUp)
    case 'y'
        fprintf('Okay... continuing...\n');
    otherwise
        fprintf('Okay, stopping here. Make sure to implement saving.\n');
        return;
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

numRegions = numel(OCUE.sample) - 1;

figRidgeExtract = figure('Position',[100 100 1400 800], 'Color', white);
figRidgeExtract.Name = 'Optical ridge extraction';

for i = 1:numRegions
    
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
f1.Name = 'Respiration test based on IMF4 Ochan2';
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








