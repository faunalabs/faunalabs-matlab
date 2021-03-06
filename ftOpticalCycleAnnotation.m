%%  ftOpticalCycleAnnotation.m

%   Written by Dave Haas between 9 July & 16 October 2021

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

% define makeFilterPlots if you want to see them during early steps...

makeFilterPlots = false;

% store screenSize now for later use in auto-sized plots

screenSize = get(0,'screensize');

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
                        if (exist('E1_KE','var'))
                           thisKE = E1_KE;
                           validKE = true;
                        else 
                            validKE = false;
                        end
                        if (exist('E1_OCC','var'))
                            thisOCC = E1_OCC;
                            validOCC = true;
                        else
                            validOCC = false;
                        end
                     
                        
                    otherwise
                        analysisRange = EPOCH.preApnea;
                        thisEpoch.Time = KINEMATICS.Time(analysisRange);
                        validEpoch = true;
                        %oT = OPTICS.Time(analysisRange);
                        %oT_s = OPTICS.time_s(analysisRange);
                        
                end                
            else
                analysisRange = EPOCH.preApnea;
                thisEpoch.Time = KINEMATICS.Time(analysisRange);
                fprintf('EPOCH not found. Analyze EPOCH first or choose a new epoch!\n');
                validEpoch = true;       
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
                        if (exist('E2_OCC','var'))
                            thisOCC = E2_OCC;
                            validOCC = true;
                        else
                            validOCC = false;
                        end                        
                    otherwise
                        analysisRange = EPOCH.t1;
                        %fprintf('Skipping re-analysis. Choose a new epoch!\n');
                        %validEpoch = false;
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
                        if (exist('E3_OCC','var'))
                            thisOCC = E3_OCC;
                            validOCC = true;
                        else
                            validOCC = false;
                        end
                        
                    otherwise
                        analysisRange = EPOCH.apnea;
                        thisEpoch.Time = KINEMATICS.Time(analysisRange);
                        %fprintf('Choose a new epoch!\n');
                        validEpoch = true;
                end                
            else
                analysisRange = EPOCH.apnea;
                thisEpoch.Time = KINEMATICS.Time(analysisRange);
                validEpoch= true;
                %fprintf('EPOCH not found. Analyze EPOCH first or choose a new epoch!\n');
                %validEpoch = false;                     
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
                        if (exist('E4_OCC','var'))
                            thisOCC = E4_OCC;
                            validOCC = true;
                        else
                            validOCC = false;
                        end
                        
                    otherwise
                        analysisRange = EPOCH.t2;
                        % fprintf('Choose a new epoch!\n');
                        % validEpoch = false;
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
                        if (exist('E5_OCC','var'))
                            thisOCC = E5_OCC;
                            validOCC = true;
                        else
                            validOCC = false;
                        end
                        
                    otherwise
                        analysisRange = EPOCH.postApnea;
                        thisEpoch.Time = KINEMATICS.Time(analysisRange);
                        %fprintf('Choose a new epoch!\n');
                        validEpoch = true;
                end                
            else
                analysisRange = EPOCH.postApnea;
                thisEpoch.Time = KINEMATICS.Time(analysisRange);
                validEpoch = true;
                %fprintf('EPOCH not found. Analyze EPOCH first or choose a new epoch!\n');
                %validEpoch = false;               
            end
            
        otherwise
            fprintf('That was not a valid epoch choice. Try again.\n');
    end    
end

% fix seconds rounding issue in thisCC while we're here...

if (validCC)
    thisCC.Time.Second = round(thisCC.Time.Second, 2);
end




%% make band-pass kinematics for this region

kFs = 100;

Time = thisEpoch.Time;
kT_s = KINEMATICS.time_s(analysisRange);

if (exist('thisKE','var'))
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

% also add in depth and temp...

depth100 = PRESSURE.depth(analysisRange);
temp100  = PRESSURE.temperature(analysisRange);

%% prepare optical data in this analysisRange

fprintf('Preparing optical data...\n');

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

% dsLed1 = downsample(oLed1, o50Hz);
% dsLed2 = downsample(oLed2, o50Hz);
% dsLed3 = downsample(oLed3, o50Hz);
% dsLed4 = downsample(oLed4, o50Hz);


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
        oLowCut            = 0.2;    % ~30 beats per minute
        % oLowCut            = 0.5;    % ~30 beats per minute
        oHighCut           = 3.5;   % ~210 beats per minute <- 3.5 test 12.5
        %oHighCut           = 5.5;   % ~210 beats per minute <- 3.5 test 12.5
        %oHighCut           = 7.5;   % ~210 beats per minute <- 3.5 test 12.5
        %oHighCut = 60;
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


figOptics = figure('Color', white);
figSize = figOptics.Position;
figSize(3) = screenSize(3);
figSize(4) = screenSize(4) * 0.7;
figOptics.Position = figSize;  
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


cleanLed1 = deOutFiltLed1;
cleanLed2 = deOutFiltLed2;
cleanLed3 = deOutFiltLed3;
cleanLed4 = deOutFiltLed4;

d1FiltLed1 = diff(dFiltLed1);
d1FiltLed2 = diff(dFiltLed2);
d1FiltLed3 = diff(dFiltLed3);
d1FiltLed4 = diff(dFiltLed4);

d1FiltLed1(end+1) = d1FiltLed1(end);
d1FiltLed2(end+1) = d1FiltLed2(end);
d1FiltLed3(end+1) = d1FiltLed3(end);
d1FiltLed4(end+1) = d1FiltLed4(end);

figDeOut = figure('Color', white);
figSize = figDeOut.Position;
figSize(3) = screenSize(3);
figSize(4) = screenSize(4) * 0.7;
figDeOut.Position = figSize;  
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

%% Sovitzky-Golay

sgLed1_30_401 = sgolayfilt(dLed1, 30, 401);
sgLed2_30_401 = sgolayfilt(dLed2, 30, 401);
sgLed3_30_401 = sgolayfilt(dLed3, 30, 401);
sgLed4_30_401 = sgolayfilt(dLed4, 30, 401);

sgFiltLed1 = sgLed1_30_401 - dLed1;
sgFiltLed2 = sgLed2_30_401 - dLed2;
sgFiltLed3 = sgLed3_30_401 - dLed3;
sgFiltLed4 = sgLed4_30_401 - dLed4;

%% decimate kinematics and depth and temp

dbpAx = downsample(bpAx, k50Hz);
dbpAy = downsample(bpAy, k50Hz);
dbpAz = downsample(bpAz, k50Hz);
dbpGx = downsample(bpGx, k50Hz);
dbpGy = downsample(bpGy, k50Hz);
dbpGz = downsample(bpGz, k50Hz);
dbpOdba = downsample(bpOdba, k50Hz);
dbpOdav = downsample(bpOdav, k50Hz);

dPitch = downsample(KINEMATICS.pitch(analysisRange), k50Hz);
dRoll = downsample(KINEMATICS.roll(analysisRange), k50Hz);
dHeading = downsample(KINEMATICS.heading(analysisRange), k50Hz);

dBestImfBpOdba = downsample(thisEpoch.bestImfBpOdba, k50Hz);
dBestImfBpOdav = downsample(thisEpoch.bestImfBpOdav, k50Hz);
dConsensusHr = decimate(thisEpoch.consensusHr, k50Hz);

dDepth = decimate(depth100, k50Hz);
dTemp  = decimate(temp100, k50Hz);

koDiff = numel(dLed2) - numel(dbpAx);

if ( koDiff > 0) 
    
    fprintf('Decimated K and O length mismatch. Fixing...\n'); 
	fprintf('Adding %d elements to k sigs.\n', koDiff);
    
    dbpAx(end+koDiff) = dbpAx(end);
    dbpAy(end+koDiff) = dbpAy(end);
    dbpAz(end+koDiff) = dbpAz(end);
    dbpGx(end+koDiff) = dbpGx(end);
    dbpGy(end+koDiff) = dbpGy(end);
    dbpGz(end+koDiff) = dbpGz(end);
    dbpOdba(end+koDiff) = dbpOdba(end);
    dbpOdav(end+koDiff) = dbpOdav(end);
    dPitch(end+koDiff) = dPitch(end);
    dRoll(end+koDiff) = dRoll(end);
    dHeading(end+koDiff) = dHeading(end);
    dBestImfBpOdba(end+koDiff) = dBestImfBpOdba(end);
    dBestImfBpOdav(end+koDiff) = dBestImfBpOdav(end);
    dConsensusHr(end+koDiff) = dConsensusHr(end);
    dDepth(end+koDiff) = dDepth(end);
    dTemp(end+koDiff) = dTemp(end);
    
elseif ( koDiff < 0 )
    
    fprintf('Removing %d elements from k sigs.\n', koDiff);
    
    dbpAx = dbpAx(1:end-koDiff); 
    dbpAy = dbpAy(1:end-koDiff); 
    dbpAz = dbpAz(1:end-koDiff); 
    dbpGx = dbpGx(1:end-koDiff); 
    dbpGy = dbpGy(1:end-koDiff); 
    dbpGz = dbpGz(1:end-koDiff); 
    dbpOdba = dbpOdba(1:end-koDiff);
    dbpOdav = dbpOdav(1:end-koDiff);
    dPitch = dPitch(1:end-koDiff);
    dRoll = dRoll(1:end-koDiff);
    dHeading = dHeading(1:end-koDiff);
    dBestImfBpOdba = dBestImfBpOdba(1:end-koDiff);
    dBestImfBpOdav = dBestImfBpOdav(1:end-koDiff);
    dConsensusHr = dConsensusHr(1:end-koDiff);
    dDepth = dDepth(1:end-koDiff);
    dTemp = dTemp(1:end-koDiff);

else
    
    fprintf('Decimated K and O lengths match! Proceding...\n'); 
    
end
   
fprintf('Decimated optical channels and key EPOCH SCG data to 50 Hz.\n');

sampleNum(1:length(dT),1) = 1:length(dT);


%% make a decimated timetable of EVERYTHING!



fprintf('Creating DATA timetable of all data decimated to 50 Hz...\n');

Time = dT;
time_s = dT_s;

DATA = timetable(Time, dbpAx, dbpAy, dbpAz, dbpGx, dbpGy, dbpGz, ...
    dbpOdba, dbpOdav, dBestImfBpOdba, dBestImfBpOdav, ...
    dConsensusHr, dPitch, dRoll, dHeading, ...
    dDepth, dTemp, dLed1, dLed2, dLed3, dLed4, ...
    dFiltLed1, dFiltLed2, dFiltLed3, dFiltLed4, ...
    d1FiltLed1, d1FiltLed2, d1FiltLed3, d1FiltLed4, ...
    cleanLed1, cleanLed2, cleanLed3, cleanLed4, ...
    sgFiltLed1, sgFiltLed2, sgFiltLed3, sgFiltLed4, ...
    time_s, sampleNum, ...
    'VariableNames', ...
    {'ax','ay','az','gx','gy','gz', ...
    'odba','odav','bestImfOdba', 'bestImfOdav', ...
    'wsstHr', 'pitch', 'roll', 'heading', ...
    'depth', 'temp', 'led1','led2','led3','led4', ...
    'filtLed1','filtLed2','filtLed3','filtLed4', ...
    'd1FiltLed1','d1FiltLed2','d1FiltLed3','d1FiltLed4', ...
    'cleanLed1','cleanLed2','cleanLed3','cleanLed4', ...
    'sgLed1', 'sgLed2', 'sgLed3', 'sgLed4', ...
    'time_s', 'sampleNum'} );

DATA.Time.Second = round(DATA.Time.Second, 2);

head(DATA)


%% 

clear -regexp ^bp ^clean ^d1 ^dBest ^dbp ^dConsensus ^deOut ...
    ^dDepth ^dFilt ^dHead ^dLed ^dPitch ^dRoll ...
    ^fHist ^figDeOut ^figOptics ^filtLed ^iLed ^KINEMATICS...
    ^kT ^l95th ^led ^MAG ^mean ^oLed ^OPTICS ^oT ...
    ^POWER ^PRESSURE ^r95 ^rLed ^std 
     

%% 
% 
% figOptics = figure('Color', white);
% figSize = figOptics.Position;
% figSize(3) = screenSize(3);
% figSize(4) = screenSize(4) * 0.7;
% figOptics.Position = figSize;  
% 
% s1 = subplot(411);
% plot(DATA.Time, DATA.d1FiltLed2, 'Color', blue);
% 
% hold on;
% if (exist('thisCC','var'))
%     xline(thisCC.Time(string(thisCC.cueType)=='ao'), ...
%         'Color', goldenrod, 'LineWidth', 1.2);    
% end
% hold off;
% xlabel('Time, local');
% ylabel('A.U.');
% title('Optical cardiography (\lambda = 1050 nm)');
% grid;
% legend('1050 nm','cardiac cycle');
% 
% s2 = subplot(412);
% plot(DATA.Time, DATA.d1FiltLed3, 'Color', red);
% hold on;
% if (exist('thisCC','var'))
%     xline(thisCC.Time(string(thisCC.cueType)=='ao'), ...
%         'Color', goldenrod, 'LineWidth', 1.2);    
% end
% hold off;
% xlabel('Time, local');
% ylabel('A.U.');
% title('Optical cardiography (\lambda = 1200 nm)');
% grid;
% legend('1200 nm','cardiac cycle');
% 
% s3 = subplot(413);
% plot(DATA.Time, DATA.ax, 'Color', red);
% hold on;
% if (exist('thisCC','var'))
%     xline(thisCC.Time(string(thisCC.cueType)=='ao'), ...
%         'Color', goldenrod, 'LineWidth', 1.2);    
% end
% hold off;
% xlabel('Time, local');
% ylabel('m ?? s^-^2');
% title('Seismocardiography (band-passed x-axis accelerometer');
% grid;
% legend('bpAx','cardiac cycle');
% 
% s4 = subplot(414);
% plot(DATA.Time, DATA.gy, 'Color', green);
% hold on;
% if (exist('thisCC','var'))
%     xline(thisCC.Time(string(thisCC.cueType)=='ao'), ...
%         'Color', goldenrod, 'LineWidth', 1.2);    
% end
% hold off;
% xlabel('Time, local');
% ylabel('degrees ?? s^-^1');
% title('Gyrocardiography (band-passed y-axis gyroscope');
% grid;
% legend('bpGy','cardiac cycle');
% 
% linkaxes([s1 s2 s3 s4],'x');



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



% now go to work...

cardiacCyclesDefined = true;

overRideYscale = false;

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
                    startTime = 1;
                    timeWindow   = 5;                   % five seconds
                    sampleWindow = timeWindow * dFs;    % number of samples per window
                    lastSample = numel(DATA.led2);                     
            end

        otherwise       % deal with pre-started CC
            
            startHere = thisCC.Time(end);
            secondsJump = seconds(startHere - datetime(thisEpoch.Time(1)));
            sampleJump = dFs * round(secondsJump, 2);
            fprintf('Okay, starting at sample %d | time: %s\n', ...
                sampleJump, startHere);

    end                 % end startWhere switch-case
end                     % end pre-existing case of CC timetable

% Create the epoch sample window for cardiac cycle annotation

timeWindow   = 5;                   % five seconds
sampleWindow = timeWindow * dFs;    % number of samples per window
lastSample = numel(DATA.led2);    
    
if (exist('sampleJump','var'))
    startTime = round(sampleJump,0);
else 
    startTime = 1;
	timeWindow   = 5;                   % five seconds
    sampleWindow = timeWindow * dFs;    % number of samples per window
    lastSample = numel(DATA.led2); 
end

endTime = startTime + sampleWindow;
if ( endTime > numel(dT) )
    endTime = numel(dT); 
    startTime = endTime - sampleWindow;
end

thisRange = timerange(DATA.Time(startTime), DATA.Time(endTime) );

if (exist('thisCC','var'))
    useCC = true;
    ccTimes = thisCC.Time(thisRange);
    ccCues = thisCC.cueType(thisRange);
    aoTimes =  ccTimes(string(ccCues)=='ao');
else
    useCC = false;
end

if (exist('thisOCC','var'))
   useOCC = true;
   occTimes = thisOOC.Time(thisRange);
   occCues  = thisOCC.cueType(thisRange);
else
    useOCC = false;
end

figAnalysisRegion = figure('Color', white);
figSize = figAnalysisRegion.Position;
figSize(1) = screenSize(3) / 2;
figSize(2) = screenSize(4) * 0.3;
figSize(3) = screenSize(3) * 0.5;
figSize(4) = screenSize(4) * 0.8;
figAnalysisRegion.Position = figSize;  

s1 = subplot(411);
plot(DATA.Time(thisRange), DATA.d1FiltLed2(thisRange), ...
    'Color', blue, 'LineWidth', 1.2);
ylim([ -rms(DATA.d1FiltLed2(thisRange)) * 3 ...
    rms(DATA.d1FiltLed2(thisRange)) * 3 ]);
hold on;
if (exist('thisCC','var'))
    yPos = s1.YLim(2) - ((s1.YLim(2) - s1.YLim(1)) / 10);
    ftPlotCues(ccTimes, ccCues, yPos);
    xline(aoTimes, 'Color', goldenrod, 'LineWidth', 1.2);    
end
hold off;
xlabel('Time, local');
ylabel('A.U.');
title('Optical cardiography (\lambda = 1050 nm)');

grid;
legend('1050 nm','cardiac cycle', 'Location', 'southeast');

s2 = subplot(412);
plot(DATA.Time(thisRange), DATA.d1FiltLed3(thisRange), ...
    'Color', red, 'LineWidth', 1.2);
ylim([ -rms(DATA.d1FiltLed3(thisRange)) * 3 ...
    rms(DATA.d1FiltLed3(thisRange)) * 3 ]);
hold on;
if (exist('thisCC','var'))
    yPos = s2.YLim(2) - ((s2.YLim(2) - s2.YLim(1)) / 10);
    ftPlotCues(ccTimes, ccCues, yPos);    
    xline(aoTimes, 'Color', goldenrod, 'LineWidth', 1.2);    
end
hold off;
xlabel('Time, local');
ylabel('A.U.');
title('Optical cardiography (\lambda = 1200 nm)');
grid;
legend('1200 nm','cardiac cycle', 'Location', 'southeast');

s3 = subplot(413);
plot(DATA.Time(thisRange), DATA.ax(thisRange), 'Color', maroon);
hold on;
if (exist('thisCC','var'))
    yPos = s3.YLim(2) - ((s3.YLim(2) - s3.YLim(1)) / 10);
    ftPlotCues(ccTimes, ccCues, yPos);    
    xline(aoTimes, 'Color', goldenrod, 'LineWidth', 1.2);    
end
hold off;
xlabel('Time, local');
ylabel('m ?? s^-^2');
title('Seismocardiography (band-passed x-axis accelerometer');
grid;
legend('bpAx','cardiac cycle', 'Location', 'southeast');

s4 = subplot(414);
plot(DATA.Time(thisRange), DATA.depth(thisRange), 'Color', cyan);
% hold on;
% if (exist('thisCC','var'))
%     yPos = s4.YLim(2) - ((s4.YLim(2) - s4.YLim(1)) / 10);
%     ftPlotCues(ccTimes, ccCues, yPos);    
%     xline(aoTimes, 'Color', goldenrod, 'LineWidth', 1.2);    
% end
% hold off;
xlabel('Time, local');
ylabel('Tag depth, m');
title('Tag depth (positive means the tag is above water level)');
grid;
%legend('depth','cardiac cycle');

linkaxes([s1 s2 s3 s4],'x');


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
    fprintf('\t[1] - 105nm only            | [2] - 1200nm only\n');
    fprintf('\t[3] - 1050 & 1200 nm both   | [z] - neither channel\n');
    fprintf('\t[a] - aortic valve open     | [x] - exclusion zone\n');
    fprintf('\t[v] - ventilation           | [p] - probable optic vent detect\n');
    fprintf('\t[!] - all AO == 3           | [0] - all AO == z\n');
    fprintf('\t[u] - undo last cue         | [r] - remove cue(s) within L-R range\n');
    fprintf('\t[b] - back five seconds     | [f] - forward five seconds\n');
    fprintf('\t[,] - back by one second    | [.] - forward by one second\n');
    fprintf('\t[<] - go to epoch start     | [>] - go to epoch end\n');
    fprintf('\t[w] - adjust window size    | [y] - adjust ?? y-axis scale\n');
    fprintf('\t[j] - jump to (hh:mm:ss)    | [d] - clean & dump thisCC to console\n');
    fprintf('\t[s] - save & end annotate   | [q] - quit w/o save\n');
    
    
    [x,~,button] = ginput(1);
    
    xAxis = gca;
    
    switch(button)
        case '1'                                    % 1 - 1050 nm pulse detection
            Time(i) = num2ruler(x, xAxis.XAxis);
            cueType(i) = {'1'};
            
            fprintf('1 -> time: %s | sample: %d\n', Time(i), dT_s(i));
            fprintf('Remaining: %s\n', duration(thisEpoch.Time(end) - Time(end)))
            thisCC = timetable(Time',cueType','VariableNames',{'cueType'});
            thisCC = sortrows(thisCC);
            i = i + 1;
            adjustTime = false;
            updatePlot = true;
            annotationComplete = false;
        case '2'                                    % 1 - 1200 nm pulse detection
            Time(i) = num2ruler(x, xAxis.XAxis);
            cueType(i) = {'2'};
            
            fprintf('2 -> time: %s | sample: %d\n', Time(i), dT_s(i));
            fprintf('Remaining: %s\n', duration(thisEpoch.Time(end) - Time(end)))
            thisCC = timetable(Time',cueType','VariableNames',{'cueType'});
            thisCC = sortrows(thisCC);
            i = i + 1;
            adjustTime = false;
            updatePlot = true;
            annotationComplete = false;
        case '3'                                    % 3 - 1050 nm AND 1200 nm pulse detection
            Time(i) = num2ruler(x, xAxis.XAxis);
            cueType(i) = {'3'};
            
            fprintf('3 -> time: %s | sample: %d\n', Time(i), dT_s(i));
            fprintf('Remaining: %s\n', duration(thisEpoch.Time(end) - Time(end)))
            thisCC = timetable(Time',cueType','VariableNames',{'cueType'});
            thisCC = sortrows(thisCC);
            i = i + 1;
            adjustTime = false;
            updatePlot = true;
            annotationComplete = false;            
        case 'p'                                    % P - probable optic vent detection
            Time(i) = num2ruler(x, xAxis.XAxis);
            cueType(i) = {'P'};
            
            fprintf('P -> time: %s | sample: %d\n', Time(i), dT_s(i));
            fprintf('Remaining: %s\n', duration(thisEpoch.Time(end) - Time(end)))
            thisCC = timetable(Time',cueType','VariableNames',{'cueType'});
            thisCC = sortrows(thisCC);
            i = i + 1;
            adjustTime = false;
            updatePlot = true;
            annotationComplete = false;
        case 'z'                                    % z - neither 1050 nm & 1200 nm have pulsatile signal
            Time(i) = num2ruler(x, xAxis.XAxis);
            cueType(i) = {'z'};
            
            fprintf('z -> time: %s | sample: %d\n', Time(i), dT_s(i));
            fprintf('Remaining: %s\n', duration(thisEpoch.Time(end) - Time(end)))
            thisCC = timetable(Time',cueType','VariableNames',{'cueType'});
            thisCC = sortrows(thisCC);
            i = i + 1;
            adjustTime = false;
            updatePlot = true;
            annotationComplete = false;            
       case 'a'                                    % AO - aortic valve open
            Time(i) = num2ruler(x, xAxis.XAxis);
            cueType(i) = {'ao'};
            fprintf('AO: %s\n', Time(i));
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
            fprintf('V -> time: %s | sample: %d\n', Time(i), dT_s(i));
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
                round(seconds(datetime(jumpTo) - thisCC.Time(1)),2) * dFs;
            startTime = sampleJump;
            adjustTime = true;
            updatePlot = true;
            annotationComplete = false;
            
        case 1
            
            thisJumpTime = num2ruler(x, xAxis.XAxis);
            sampleJump = round(seconds(datetime(thisJumpTime) - ...
                thisCC.Time(1)),2) * dFs;
            startTime = sampleJump;
            adjustTime = true;
            updatePlot = true;
            annotationComplete = false;            
            
        case 'w'
            newWinSizeTxt = 'Enter new window size (default is 5s): ';
            winSize = double(input(newWinSizeTxt));
            fprintf('New window size set to %d\n', winSize);
            timeWindow   = winSize;                   % five seconds
            sampleWindow = timeWindow * dFs;    % number of samples per window
            startTime = startTime;
            endTime = startTime + sampleWindow;
            adjustTime = true;
            updatePlot = true;
            annotationComplete = false;   
        case 'y'
            if (overRideYscale == true)
               fprintf('Resetting auto y-scaling...\n');
               overRideYscale = false;   
            else
                newYtxt = 'Enter the new ?? y-scale (single number): ';
                newY = double(input(newYtxt));
                overRideYscale = true;
            end
            adjustTime = true;
            updatePlot = true;
            annotationComplete = false;
            
        case 'o'
            fprintf('Choose one of the following variables to display...\n');
            fprintf('\t1 = led* (unfilter but decimated)\n');
            fprintf('\t2 = filtLed* (band-pass filtered)\n');
            fprintf('\t3 = d1FiltLed* (1?? band-pass filtered)\n');
            fprintf('\t4 = cleanLed* (outlier & cmddenoised)\n');
            fprintf('\t5 = sgFilt* (Savitzky-Golay filtered)\n');
            whichOptVar = input('Enter the number of the var type: ');
            switch( double(whichOptVar) )
                case 1
                    
                case 2
                    
                case 3
                    
                case 4
                    
                case 5
                    
                otherwise
                    fprintf('Bad entry type. Sticking with what ya got already.\n');
            end
            
            
            
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
            startTime = startTime + dFs;
            adjustTime = true;
            updatePlot = true;
            annotationComplete = false;
        case ','                                    % move -1 second window
            fprintf('[,] rewinding the window by one second...\n');
            timeJump = seconds(-1);
            startTime = startTime - dFs;
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

        case '!'
            
            insertTimes = thisCC.Time( string(thisCC.cueType)=='ao' | ...
                string(thisCC.cueType) == 'x' );
            allOpticsTimes = insertTimes + seconds(0.1);
            numObjs = numel(allOpticsTimes);
            allOpticsCues(1:numObjs,1) = {'3'};
            newCCTimes = [thisCC.Time ; allOpticsTimes ];
            newCCCues  = [thisCC.cueType; allOpticsCues ];

            if (exist('Time','var'))
                tempTime = Time;
                clear Time;
                Time = newCCTimes;    
            else
                Time = newCCTimes;    
            end

            ccTT = timetable(Time, newCCCues, 'VariableNames', {'cueType'} );
            ccTT = sortrows(ccTT);

            if (exist('ccTT','var'))
                cardiacCyclesDefined = true;
                clear natPositive;
                clear natIndex;
                clear numNats;

                fprintf('Checking for NATs in thisCC before proceeding...\n');
                q = 1;
                for n = 1:height(ccTT)
                   if ( isnat(ccTT.Time(n)) )
                        fprintf('Found a NAT in ccTT row %d. Adding to natIndex...\n', ...
                            n);
                        natIndex(q) = n;
                        q = q + 1;
                   end
                end
                if (~exist('natIndex','var'))
                    fprintf('\tNo NATs detected in ccTT. Proceeding!\n')
                    natIndex = 0;
                else
                    numNats = numel(natIndex);
                    if ( numNats > 1 )
                        fprintf('NATs are consecutive. Dropping last %d rows.\n', ...
                            numel(natIndex) );
                        newEnd = natIndex(1) - 1;
                        Time = ccTT.Time(1:newEnd);
                        cueType = ccTT.cueType(1:newEnd);
                        ccTT = timetable(Time,cueType,'VariableNames',{'cueType'});
                        ccTT = sortrows(ccTT);
                    elseif (numNats == 1)
                        fprintf('NAT last entry is NAT. Fixing...\n');
                        newEnd = natIndex(1) - 1;
                        Time = ccTT.Time(1:newEnd);
                        cueType = ccTT.cueType(1:newEnd);
                        ccTT= timetable(Time,cueType,'VariableNames',{'cueType'});
                        ccTT = sortrows(ccTT);        
                    elseif (numNats == 0)
                        fprintf('\tNo NATs in ccTT! Proceeding...\n');
                    else
                        fprintf('NATs are in nonconsecutive order. Handle manually.\n');
                        return;
                    end
                end
            end
            disp(ccTT)
            txt = 'Proceeding with replacement of thisCC with contents of ccTT (y/n) [n]: ';
            replaceStr = input(txt, 's');
            switch(replaceStr)
                case 'y'
                    fprintf('Okay... replacing thisCC with ccTT... \n');
                    thisCC = ccTT;
                otherwise
                    fprintf('Okay... doing nothing... \n');
            end            
            adjustTime = false;
            updatePlot = true;
            annotationComplete = false;            
            
        case '0'
            
            insertTimes = thisCC.Time( string(thisCC.cueType)=='ao' | ...
                string(thisCC.cueType) == 'x' );
            noOpticsTimes = insertTimes + seconds(0.1);
            numObjs = numel(noOpticsTimes);
            noOpticsCues(1:numObjs,1) = {'z'};
            newCCTimes = [thisCC.Time ;noOpticsTimes ];
            newCCCues  = [thisCC.cueType; noOpticsCues ];

            if (exist('Time','var'))
                tempTime = Time;
                clear Time;
                Time = newCCTimes;    
            else
                Time = newCCTimes;    
            end

            ccTT = timetable(Time, newCCCues, 'VariableNames', {'cueType'} );
            ccTT = sortrows(ccTT);

            if (exist('ccTT','var'))
                cardiacCyclesDefined = true;
                clear natPositive;
                clear natIndex;
                clear numNats;

                fprintf('Checking for NATs in thisCC before proceeding...\n');
                q = 1;
                for n = 1:height(ccTT)
                   if ( isnat(ccTT.Time(n)) )
                        fprintf('Found a NAT in ccTT row %d. Adding to natIndex...\n', ...
                            n);
                        natIndex(q) = n;
                        q = q + 1;
                   end
                end
                if (~exist('natIndex','var'))
                    fprintf('\tNo NATs detected in ccTT. Proceeding!\n')
                    natIndex = 0;
                else
                    numNats = numel(natIndex);
                    if ( numNats > 1 )
                        fprintf('NATs are consecutive. Dropping last %d rows.\n', ...
                            numel(natIndex) );
                        newEnd = natIndex(1) - 1;
                        Time = ccTT.Time(1:newEnd);
                        cueType = ccTT.cueType(1:newEnd);
                        ccTT = timetable(Time,cueType,'VariableNames',{'cueType'});
                        ccTT = sortrows(ccTT);
                    elseif (numNats == 1)
                        fprintf('NAT last entry is NAT. Fixing...\n');
                        newEnd = natIndex(1) - 1;
                        Time = ccTT.Time(1:newEnd);
                        cueType = ccTT.cueType(1:newEnd);
                        ccTT= timetable(Time,cueType,'VariableNames',{'cueType'});
                        ccTT = sortrows(ccTT);        
                    elseif (numNats == 0)
                        fprintf('\tNo NATs in ccTT! Proceeding...\n');
                    else
                        fprintf('NATs are in nonconsecutive order. Handle manually.\n');
                        return;
                    end
                end
            end
            disp(ccTT)
            txt = 'Proceeding with replacement of thisCC with contents of ccTT (y/n) [n]: ';
            replaceStr = input(txt, 's');
            switch(replaceStr)
                case 'y'
                    fprintf('Okay... replacing thisCC with ccTT... \n');
                    thisCC = ccTT;
                otherwise
                    fprintf('Okay... doing nothing... \n');
            end            
            adjustTime = false;
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
                                E1_DECIMATED = DATA;
                                E1_KE = thisKE;
                                t
                                fprintf('E1_CC added to RAW.\n');
                            case 2
                                E2_DECIMATED = DATA;
                                E2_KE = thisKE;
                                fprintf('\nAppending E2_CC timetable to RAW file...\n');
                                save(fullRawFileName, 'E2_CC', 'E2_KE', 'E2_DECIMATED', '-append');               
                                fprintf('E2_CC added to RAW.\n');
                            case 3
                                E3_DECIMATED = DATA;
                                E3_KE = thisKE;
                                fprintf('\nAppending E3_CC timetable to RAW file...\n');
                                save(fullRawFileName, 'E3_CC', 'E3_KE', 'E3_DECIMATED', '-append');               
                                fprintf('E3_CC added to RAW.\n');
                            case 4
                                E4_DECIMATED = DATA;
                                E4_KE = thisKE;
                                fprintf('\nAppending E4_CC timetable to RAW file...\n');
                                save(fullRawFileName, 'E4_CC', 'E4_KE', 'E4_DECIMATED', '-append');
                                fprintf('E4_CC added to RAW.\n');
                            case 5
                                E5_DECIMATED = DATA;
                                E5_KE = thisKE;
                                fprintf('\nAppending E5_CC timetable to RAW file...\n');
                                save(fullRawFileName, 'E5_CC', 'E5_KE', 'E5_DECIMATED', '-append');
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
        if ( endTime > numel(dT) )
            endTime = numel(dT); 
            startTime = endTime - (dFs * 5);
        end
        updatePlot = true;      
        
    end
    
    if (updatePlot)
        
        thisRange = timerange(DATA.Time(uint32(startTime)), ...
            DATA.Time(uint32(endTime)) );

        useCC = true;
        ccTimes = thisCC.Time(thisRange);
        ccCues = thisCC.cueType(thisRange);
            aoTimes =  ccTimes(string(ccCues)=='ao');
   
        figure(figAnalysisRegion);    
        s1 = subplot(411);
        plot(DATA.Time(thisRange), DATA.d1FiltLed2(thisRange), ...
            'Color', blue, 'LineWidth', 1.2);
        if (overRideYscale)
            yScale = newY;
            ylim([ (-1 * yScale)  yScale ]); 
        else
            ylim([ -rms(DATA.d1FiltLed2(thisRange)) * 3 ...
                rms(DATA.d1FiltLed2(thisRange)) * 3 ]);                    
        end     
        hold on;
        if (exist('thisCC','var'))
            yPos = s1.YLim(2) - ((s1.YLim(2) - s1.YLim(1)) / 10);
            ftPlotCues(ccTimes, ccCues, yPos);  
            if (~isempty(aoTimes))
                xline(aoTimes, 'Color', goldenrod, 'LineWidth', 1.2);    
            end
            
        end
        hold off;
        xlabel('Time, local');
        ylabel('A.U.');
        title('Optical cardiography (\lambda = 1050 nm)');
        grid;
        legend('1050 nm','cardiac cycle', 'Location', 'southeast');

        s2 = subplot(412);
        plot(DATA.Time(thisRange), DATA.d1FiltLed3(thisRange), ...
            'Color', red, 'LineWidth', 1.2);
        if (overRideYscale)
            yScale = newY;
            ylim([ (-1 * yScale)  yScale ]); 
        else
            ylim([ -rms(DATA.d1FiltLed3(thisRange)) * 3 ...
                rms(DATA.d1FiltLed3(thisRange)) * 3 ]);                    
        end
        

        hold on;
        if (exist('thisCC','var'))
            yPos = s2.YLim(2) - ((s2.YLim(2) - s2.YLim(1)) / 10);
            ftPlotCues(ccTimes, ccCues, yPos);            
            if (~isempty(aoTimes))
                xline(aoTimes, 'Color', goldenrod, 'LineWidth', 1.2);    
            end 
        end
        hold off;
        xlabel('Time, local');
        ylabel('A.U.');
        title('Optical cardiography (\lambda = 1200 nm)');
        grid;
        legend('1200 nm','cardiac cycle', 'Location', 'southeast');

        s3 = subplot(413);
        plot(DATA.Time(thisRange), DATA.ax(thisRange), 'Color', red);
        hold on;
        if (exist('thisCC','var'))
            yPos = s3.YLim(2) - ((s3.YLim(2) - s3.YLim(1)) / 10);
            ftPlotCues(ccTimes, ccCues, yPos);            
            if (~isempty(aoTimes))
                xline(aoTimes, 'Color', goldenrod, 'LineWidth', 1.2);    
            end
        end
        hold off;
        xlabel('Time, local');
        ylabel('m ?? s^-^2');
        title('Seismocardiography (band-passed x-axis accelerometer');
        grid;
        legend('bpAx','cardiac cycle', 'Location', 'southeast');

        s4 = subplot(414);
        plot(DATA.Time(thisRange), DATA.depth(thisRange), 'Color', cyan);
%         hold on;
%         if (exist('thisCC','var'))
%             yPos = s4.YLim(2) - ((s4.YLim(2) - s4.YLim(1)) / 10);
%             ftPlotCues(ccTimes, ccCues, yPos);            
%             if (~isempty(aoTimes))
%                 xline(aoTimes, 'Color', goldenrod, 'LineWidth', 1.2);    
%             end   
%         end
%         hold off;
        xlabel('Time, local');
        ylabel('Depth, meters');
        title('Tag depth (positive values mean tag is above water level)');
        grid;
        %legend('bpGy','cardiac cycle');
        

        linkaxes([s1 s2 s3 s4],'x');
    
        adjustTime = false;
    
    end
    
end             % end annotationComplete test


