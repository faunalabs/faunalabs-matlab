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

%makePlots = lower(input('Do you want to see diagnostic plots (y/n) [n]: ','s'));
makePlots = 'n';
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
        
        if (exist('OPTICS_METADATA.fs','var'))
            oFs = OPTICS_METADATA.Properties.SampleRate;
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
    %goodRegionTxt = 'Is this region good for your analysis? (y/n) [y]: ';
    %thatsGood = lower(input(goodRegionTxt,'s'));        

    thatsGood = 'y';
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
   
   %okToProceed = lower(input('Okay to proceed (y/n) [y]: ','s'));
   okToProceed = 'y';
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

if (~exist('makePlots','var'))
    makePlots = false;
end

switch(epochChoice)
    case 3                  % apnea
        bandPassFiltOrder  = 3;      % 3rd order butterworth filter
        oLowCut            = 0.5;    % ~60 beats per minute
        oHighCut           = 7.5;    % ~90 beats per minute
        opticalTestFs      = oHighCut * 2;
        %opticalTestFs      = 8;      % 8 produces good pulses on optics apnea        
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

if (makePlots)
    fvtO = fvtool(sosBpo, bandPassOpticalFilter, 'Fs', opticalTestFs);
    legend(fvtO,'butter','designfilt');
end

OPTICS_METADATA = struct;
OPTICS_METADATA.BANDPASS_FILTER.filterOrder = bandPassFiltOrder;
OPTICS_METADATA.BANDPASS_FILTER.lowCut = oLowCut;
OPTICS_METADATA.BANDPASS_FILTER.highCut = oHighCut;
OPTICS_METADATA.BANDPASS_FILTER.sampleRate = opticalTestFs;
OPTICS_METADATA.BANDPASS_FILTER.A = Aobp;
OPTICS_METADATA.BANDPASS_FILTER.B = Bobp;
OPTICS_METADATA.BANDPASS_FILTER.C = Cobp;
OPTICS_METADATA.BANDPASS_FILTER.D = Dobp;
OPTICS_METADATA.BANDPASS_FILTER.SOS = sosBpo;
OPTICS_METADATA.BANDPASS_FILTER.G = gBpo;

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

%replaceTest = lower(input('Use outlier-smoothed optical data (y/n) [y]: ','s'));
replaceTest = 'y';
switch(replaceTest)
    case 'n'
        
        outliersSmoothed = false;
        OPTICS_METADATA.outlierSmoothed = false;
        fprintf('Leaving original band-passed optical data in place.\n');     
        
    otherwise
        
        outliersSmoothed = true;
        OPTICS_METADATA.outlierSmoothed = true;
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

%%


% set a 100 ms windowSize to use for gaussian filtering
%optGaussianWindowSize = oFs / 10;      % = 0.040 s | 25
%optGaussianWindowSize = oFs / 25;       % = 0.100 s | 10
optGaussianWindowSize = 7;  

gfBo = ones(1, optGaussianWindowSize) / optGaussianWindowSize;
gfA = [1];

gf_dRssDiffFiltLed2 = filtfilt(gfBo, gfA, dRssDiffFiltLed2);
gf_dRssDiffFiltLed3 = filtfilt(gfBo, gfA, dRssDiffFiltLed3);
gf_sumActiveOptics = gf_dRssDiffFiltLed2 + gf_dRssDiffFiltLed3;
log10GfActiveOptics = real(log10(gf_sumActiveOptics));
log10GfActiveOptics(end) = log10GfActiveOptics(end-1);
normGfActiveOptics = normalize(gf_sumActiveOptics);

[wsstThis, fThis] = wsst(normGfActiveOptics, 50);
[ridgeThis, iridgeThis] = wsstridge(wsstThis, 20, fThis);

switch(epochChoice)
    case 1
        thisCC          = E1_CC;
        thisCC_ENERGY   = E1_CC_ENERGY;
        thisK           = E1_KINEMATICS;
        thisVENT        = E1_VENT_DATA;
    case 3
        thisCC          = E3_CC;
        thisCC_ENERGY   = E3_CC_ENERGY;
        thisK           = E3_KINEMATICS;
        
        
        appendVent = timetable(thisCC.Time(1), {'v'}, ...
            'VariableNames', {'cueType'} )
        thisCC = [thisCC; appendVent]
        thisCC = sortrows(thisCC);
        
        
    case 5
        thisCC          = E5_CC;
        thisCC_ENERGY   = E5_CC_ENERGY;
        thisK           = E5_KINEMATICS;
        thisVENT        = E5_VENT_DATA;
end

%
figure('Color', white);

s3a = subplot(211);
plot(dT, normGfActiveOptics, 'Color', black);
hold on;
xline(thisCC.Time((string(thisCC.cueType)=='ao')), '-', ...
    'Color', goldenrod, 'LineWidth', 2);
xline(thisCC.Time((string(thisCC.cueType)=='x')), '--', ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
xline(thisCC.Time((string(thisCC.cueType)=='v')), '-', ...
    'Color', green, 'LineWidth', 2);
yline(mean(normGfActiveOptics), '--', 'Color', purple);
hold off;
grid;

s3b = subplot(212);
plot(dT, ridgeThis, '-', 'Color', blue);
hold on;
plot(thisK.Time, thisK.consensusHr, '-', 'Color', red);
hold off;
grid;
legend('optics','ensWsstHr');


%%

figure; 

r1a = subplot(411);
plot(dT, gf_dRssDiffFiltLed2, 'Color', blue);
hold on;
xline(thisCC.Time((string(thisCC.cueType)=='ao')), '-', ...
    'Color', goldenrod, 'LineWidth', 2);
xline(thisCC.Time((string(thisCC.cueType)=='x')), '--', ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
xline(thisCC.Time((string(thisCC.cueType)=='v')), '-', ...
    'Color', green, 'LineWidth', 2);
hold off;
grid;

r1b = subplot(412);
plot(dT, gf_dRssDiffFiltLed3, 'Color', red);
hold on;
xline(thisCC.Time((string(thisCC.cueType)=='ao')), '-', ...
    'Color', goldenrod, 'LineWidth', 2);
xline(thisCC.Time((string(thisCC.cueType)=='x')), '--', ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
xline(thisCC.Time((string(thisCC.cueType)=='v')), '-', ...
    'Color', green, 'LineWidth', 2);
hold off;
grid;

r1c = subplot(413);
plot(dT, normGfActiveOptics, 'Color', black);
hold on;
xline(thisCC.Time((string(thisCC.cueType)=='ao')), '-', ...
    'Color', goldenrod, 'LineWidth', 2);
xline(thisCC.Time((string(thisCC.cueType)=='x')), '--', ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
xline(thisCC.Time((string(thisCC.cueType)=='v')), '-', ...
    'Color', green, 'LineWidth', 2);
yline(mean(normGfActiveOptics), '--', 'Color', purple);
hold off;
grid;

r1d = subplot(414);
bar(thisCC_ENERGY.Time, [thisCC_ENERGY.iLinKE thisCC_ENERGY.iRotKE], ...
    'stacked' );
hold on;
xline(thisCC.Time((string(thisCC.cueType)=='ao')), '-', ...
    'Color', goldenrod, 'LineWidth', 2);
xline(thisCC.Time((string(thisCC.cueType)=='v')), '-', ...
    'Color', green, 'LineWidth', 2);
hold off;
xlabel('Time, seconds');
ylabel('J');
title('Cardiac cycle kinetic energy'); 
grid;

linkaxes([r1a r1b r1c r1d],'x');

%%

figure;

s1 = subplot(311);
plot(dT, dRssDiffFiltLed2, 'Color', blue);
hold on;
xline(thisCC.Time((string(thisCC.cueType)=='ao')), '-', ...
    'Color', goldenrod, 'LineWidth', 2);
xline(thisCC.Time((string(thisCC.cueType)=='v')), '-', ...
    'Color', green, 'LineWidth', 2);
hold off;
xlabel('Time, seconds');
ylabel('A.U.');
title('Decimated 1st Derivative Band pass-filtered 1050 nm optical trace');
grid;
legend('1050 nm','SCG+GCG CC','ventilation');

s2 = subplot(312);
plot(dT, dRssDiffFiltLed3, 'Color', red);
hold on;
xline(thisCC.Time((string(thisCC.cueType)=='ao')), '-', ...
    'Color', goldenrod, 'LineWidth', 2);
xline(thisCC.Time((string(thisCC.cueType)=='v')), '-', ...
    'Color', green, 'LineWidth', 2);
hold off;
xlabel('Time, seconds');
ylabel('A.U.');
title('Decimated 1st Derivative Band pass-filtered 1200 nm optical trace');
grid;
legend('1200 nm','SCG+GCG CC','ventilation');

s3 = subplot(313);
bar(thisCC_ENERGY.Time, [thisCC_ENERGY.iLinKE thisCC_ENERGY.iRotKE], ...
    'stacked' );
hold on;
xline(thisCC.Time((string(thisCC.cueType)=='ao')), '-', ...
    'Color', goldenrod, 'LineWidth', 2);
xline(thisCC.Time((string(thisCC.cueType)=='v')), '-', ...
    'Color', green, 'LineWidth', 2);
hold off;
xlabel('Time, seconds');
ylabel('J');
title('Cardiac cycle kinetic energy'); 
grid;

linkaxes([s1 s2 s3],'x');

%%

figure;

s2a = subplot(211);
plot(dT, normalize(real(log10GfActiveOptics)), 'Color', blue);
hold on;
xline(thisCC.Time((string(thisCC.cueType)=='ao')), '-', ...
    'Color', goldenrod, 'LineWidth', 2);
xline(thisCC.Time((string(thisCC.cueType)=='v')), '-', ...
    'Color', green, 'LineWidth', 2);
hold off;
xlabel('Time, local');
ylabel('A.U.');
grid;
s2b = subplot(212);
bar(thisCC_ENERGY.Time, thisCC_ENERGY.iTotalKE, 'FaceColor', red);
hold on;
xline(thisCC.Time((string(thisCC.cueType)=='ao')), '-', ...
    'Color', goldenrod, 'LineWidth', 1);
xline(thisCC.Time((string(thisCC.cueType)=='v')), '-', ...
    'Color', green, 'LineWidth', 2);
hold off;
xlabel('Time, local');
ylabel('Kinetic energy, J');
grid;
linkaxes([s2a s2b],'x');