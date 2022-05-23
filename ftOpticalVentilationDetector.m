%%  ftOpticalVentilationDetector.m

%   Written by Dave Haas between 9 July and 3 October 2021

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
                        if (exist('E1_VENT_DATA','var'))
                            thisVent = E1_VENT_DATA;
                            validVent = true;
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
                        if (exist('E2_VENT_DATA','var'))
                            thisVent = E2_VENT_DATA;
                            validVent = true;
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
                        if (exist('E3_VENT_DATA','var'))
                            thisVent = E3_VENT_DATA;
                            validVent = true;
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
                        if (exist('E4_VENT_DATA','var'))
                            thisVent = E4_VENT_DATA;
                            validVent = true;
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
                        if (exist('E5_VENT_DATA','var'))
                            thisVent = E5_VENT_DATA;
                            validVent = true;
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
   bpOdba = thisEpoch.bpOdba;
   bpOdav = thisEpoch.bpOdav;
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

%% prepare optical data in this analysisRange

oT = OPTICS.Time(analysisRange);
oT_s = OPTICS.time_s(analysisRange);

oLed1 = OPTICS.led1(analysisRange);
oLed2 = OPTICS.led2(analysisRange);
oLed3 = OPTICS.led3(analysisRange);
oLed4 = OPTICS.led4(analysisRange);

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
kLength = numel(bpAx);
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

dsLed2 = downsample(oLed2, o50Hz);
dsLed3 = downsample(oLed3, o50Hz);

% finally, decimate key values from E (i.e., E = E1 for EPOCH 1)

dBpOdba = decimate(thisEpoch.bpOdba, k50Hz);
dBpOdav = decimate(thisEpoch.bpOdav, k50Hz);
dBestImfBpOdba = decimate(thisEpoch.bestImfBpOdba, k50Hz);
dBestImfBpOdav = decimate(thisEpoch.bestImfBpOdav, k50Hz);
dConsensusHr = decimate(thisEpoch.consensusHr, k50Hz);

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
        oHighCut           = 3.5;
        % oHighCut           = 7.5;    % ~90 beats per minute
        % oHighCut           = 60;
        opticalTestFs      = oHighCut * 2;
        %opticalTestFs      = 8;      % 8 produces good pulses on optics apnea        
    otherwise               % pre- or post-apnea (E1, E2, E4, or E5)
        bandPassFiltOrder  = 3;      % 6th order butterworth filter
        oLowCut            = 0.5;    % ~30 beats per minute
        oHighCut           = 3.5;   % ~210 beats per minute <- 3.5 test 12.5
        %oHighCut           = 7.5;   % ~210 beats per minute <- 3.5 test 12.5
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

if (makeFilterPlots)
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

%%

sym4_5_deOutFiltLed2 = cmddenoise(deOutFiltLed2, 'sym4', 5);
sym4_5_deOutFiltLed3 = cmddenoise(deOutFiltLed3, 'sym4', 5);

figure('Color', white);

s1 = subplot(411);
% plot(dT, deOutFiltLed2, 'Color', blue);
plot(dT, real(log10(dLed2)), 'Color', blue);
hold on;
% plot(dT, dFiltLed2, 'Color', cyan);
% plot(dT, log10(sym4_5_deOutFiltLed2), 'Color', purple, 'LineWidth', 2);
xline(thisCC.Time(string(thisCC.cueType)=='v'), ...
    'Color', green, 'LineWidth', 1.5);

hold off;
grid;

s2 = subplot(412);
plot(dT, real(log10(dLed3)), 'Color', red);
hold on;
% plot(dT, dFiltLed3, 'Color', maroon);
% plot(dT, sym4_5_deOutFiltLed3, 'Color', maroon, 'LineWidth', 2);
xline(thisCC.Time(string(thisCC.cueType)=='v'), ...
    'Color', green, 'LineWidth', 1.5);
% xline(thisVent.Time, 'Color', green, 'LineWidth', 1.5);
hold off;
grid;

s3 = subplot(413);
plot(kT, bpAx, 'Color', red);
grid;

s4 = subplot(414);
plot(kT, bpGy, 'Color', green);
grid;

linkaxes([s1 s2 s3 s4],'x');


%%



sym4_5_deOutFiltLed2 = cmddenoise(deOutFiltLed2, 'sym4' ,5);
sym4_5_deOutFiltLed3 = cmddenoise(deOutFiltLed3, 'sym4' ,5);

