%%  ftHeartRate_Kinematics.m

%   Written by Dave Haas between 9 July abd 7 August 2021
%
%   Last updated by Dave Haas on 7 August 2021, to make wsstridges from
%   findpeaks data.

clc;
clear;

%%

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

%% Make introductions...

fprintf('Welcome to the FaunaTag kinematic heart rate analysis app.\n');


%% define makePlots

makePlot = lower(input('Want to see filter plots? (y/n) [n]: ','s'));
switch(makePlot)
    case 'y'
        makePlot = true;
    otherwise
        makePlot = false;
end


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
        tag = TRIAL.tag;
        
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
                confText = 'Analysis results already exist. Redo? (y/n): ';
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
                confText = 'Analysis results already exist. Redo? (y/n): ';
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
                confText = 'Analysis results already exist. Redo? (y/n): ';
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
                confText = 'Analysis results already exist. Redo? (y/n): ';
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


%% Do some regional selections with wsst and wsstridge finding

screenSize = get(0,'screensize');

figAnalysisRegion = figure;
figSize = figAnalysisRegion.Position;
figSize(3) = screenSize(3);
figSize(4) = screenSize(4) * 0.65;
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
    goodRegionTxt = 'Is this region good for your analysis? (y/n): ';
    thatsGood = lower(input(goodRegionTxt,'s'));        

    switch(thatsGood)
        case 'y'
            
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
            
        otherwise
            analysisRegionDefined = false;
            
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

[denoiseLed1,~,thrParams,~,bestNbOfInt] = cmddenoise(oLed1,'db3',3);
[denoiseLed2,~,thrParams,~,bestNbOfInt] = cmddenoise(oLed2,'db3',3);
[denoiseLed3,~,thrParams,~,bestNbOfInt] = cmddenoise(oLed3,'db3',3);
[denoiseLed4,~,thrParams,~,bestNbOfInt] = cmddenoise(oLed4,'db3',3);

[denoise_rLed1,~,thrParams,~,bestNbOfInt] = cmddenoise(rLed1,'db3',3);
[denoise_rLed2,~,thrParams,~,bestNbOfInt] = cmddenoise(rLed2,'db3',3);
[denoise_rLed3,~,thrParams,~,bestNbOfInt] = cmddenoise(rLed3,'db3',3);
[denoise_rLed4,~,thrParams,~,bestNbOfInt] = cmddenoise(rLed4,'db3',3);

oHampelWin = oFs / 25;      % 10 + 1 samples for optical window
oSigma = 2;                 % +/- 2 standard deviations

[cleanedLed1, o1OutlierIndex, o1Median, o1Std] = hampel(denoiseLed1, oHampelWin, oSigma);
[cleanedLed2, o2OutlierIndex, o2Median, o2Std] = hampel(denoiseLed2, oHampelWin, oSigma);
[cleanedLed3, o3OutlierIndex, o3Median, o3Std] = hampel(denoiseLed3, oHampelWin, oSigma);
[cleanedLed4, o4OutlierIndex, o4Median, o4Std] = hampel(denoiseLed4, oHampelWin, oSigma);



%% Start with a Gaussian window filter 

% set a 100 ms windowSize to use for gaussian filtering
%optGaussianWindowSize = oFs / 10;      % = 0.040 s
%kinGaussianWindowSize = kFs / 4;       % = 0.040 s
optGaussianWindowSize = oFs / 25;       % = 0.100 s
kinGaussianWindowSize = kFs / 10;       % = 0.100 s

gfBo = ones(1, optGaussianWindowSize) / optGaussianWindowSize;
gfBk = ones(1, kinGaussianWindowSize) / kinGaussianWindowSize;

gfA = [1];

gfLed1 = filtfilt(gfBo, gfA, led1);
gfLed2 = filtfilt(gfBo, gfA, led2);
gfLed3 = filtfilt(gfBo, gfA, led3);
gfLed4 = filtfilt(gfBo, gfA, led4);

gf_denoisedLed1 = filtfilt(gfBo, gfA, denoiseLed1);
gf_denoisedLed2 = filtfilt(gfBo, gfA, denoiseLed2);
gf_denoisedLed3 = filtfilt(gfBo, gfA, denoiseLed3);
gf_denoisedLed4 = filtfilt(gfBo, gfA, denoiseLed4);

gf_ODLed1 = filtfilt(gfBo, gfA,(log(max(led1)./led1)));
gf_ODLed2 = filtfilt(gfBo, gfA,(log(max(led2)./led2)));
gf_ODLed3 = filtfilt(gfBo, gfA,(log(max(led3)./led3)));
gf_ODLed4 = filtfilt(gfBo, gfA,(log(max(led4)./led4)));

gfAx = filtfilt(gfBk, gfA, kAx);
gfAy = filtfilt(gfBk, gfA, kAy);
gfAz = filtfilt(gfBk, gfA, kAz);
gfGx = filtfilt(gfBk, gfA, kGx);
gfGy = filtfilt(gfBk, gfA, kGy);
gfGz = filtfilt(gfBk, gfA, kGz);

gfOdba = filtfilt(gfBk, gfA, kOdba);
gfOdav = filtfilt(gfBk, gfA, kOdav);

fprintf('Sliding Gaussian 100 ms window filter data is prefixed as gf*\n');

%% Plot to have a look at gaussian filtered signals

if (makePlot)
    
figGaussOpt = figure;
figGaussOpt.Name = 'Gaussian sliding window filter with 100 ms window size';

g1 = subplot(211);
plot(oT_s, gfLed1, 'Color', led1Color);
hold on;
plot(oT_s, gfLed2, 'Color', led2Color);
plot(oT_s, gfLed3, 'Color', led3Color);
plot(oT_s, gfLed4, 'Color', led4Color);
hold off;
xlabel('Time, seconds)');
ylabel('Intensity');
title('Raw optical intensities smoothed with 100 ms Gaussian filtering');
legend('ambient1','1050 nm','1200 nm','ambient2');
grid;

g2 = subplot(212);
plot(oT_s, gf_ODLed1, 'Color', led1Color);
hold on;
plot(oT_s, gf_ODLed2, 'Color', led2Color);
plot(oT_s, gf_ODLed3, 'Color', led3Color);
plot(oT_s, gf_ODLed4, 'Color', led4Color);
hold off;
xlabel('Time, seconds');
ylabel('OD');
title('Optical Density (OD) smoothed with 100 ms Gaussian filtering');
legend('ambient1','1050 nm','1200 nm','ambient2');
grid;

linkaxes([g1,g2],'x');


% Look at rolled-up kinematics 

figGaussKin = figure;

gfk1 = subplot(211);
plot(kT_s, gfOdba, 'Color', blue, 'LineWidth', 1.2);
xlabel('Time, seconds');
ylabel('m/s^2');
title('Gaussing filtered (100 ms) ODBA');
grid;

gfk2 = subplot(212);
plot(kT_s, gfOdav, 'Color', red, 'LineWidth', 1.2);
xlabel('Time, seconds');
ylabel('°/s');
title('Gaussing filtered (100 ms) ODBA');
grid;

linkaxes([gfk1 gfk2],'x');

end         % end if (plotIt == true)

%% Next up, a high-pass filter for DC offset removal

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


[AkDCHP, BkDCHP, CkDCHP, DkDCHP] = ...
    butter(dcoFiltOrder, ( dcoPassbandFreqHz / (2) ), 'high' );
    % butter(dcoFiltOrder, ( dcoPassbandFreqHz / (kNyquistFs) ), 'high' );

[kSOS,kG] = ss2sos(AkDCHP, BkDCHP, CkDCHP, DkDCHP);

% design the high-pass filter. NOTE: using samplingRate = 2, which is the
% designfilt default, results in normalized frequencies

dcoOpticalFilter = designfilt('highpassiir','FilterOrder',dcoFiltOrder, ...
    'PassbandFrequency', dcoPassbandFreqHz, ... 
    'PassbandRipple', 0.2, ...
    'SampleRate', 2);

dcoKinematicFilter = designfilt('highpassiir','FilterOrder',dcoFiltOrder, ...
    'PassbandFrequency', dcoPassbandFreqHz, ... 
    'PassbandRipple', 0.2, ...
    'SampleRate', 2);

if (makePlot)
    
    fvt = fvtool(oSOS, dcoOpticalFilter, 'Fs', oFs);
    legend(fvt,'butter','designfilt');
    
    
end

dcoLed1 = filtfilt(dcoOpticalFilter', led1);
dcoLed2 = filtfilt(dcoOpticalFilter', led2);
dcoLed3 = filtfilt(dcoOpticalFilter', led3);
dcoLed4 = filtfilt(dcoOpticalFilter', led4);

dcodLed1 = filtfilt(dcoOpticalFilter', denoiseLed1);
dcodLed2 = filtfilt(dcoOpticalFilter', denoiseLed2);
dcodLed3 = filtfilt(dcoOpticalFilter', denoiseLed3);
dcodLed4 = filtfilt(dcoOpticalFilter', denoiseLed4);
 
dco_rLed1 = filtfilt(dcoOpticalFilter', rLed1);
dco_rLed2 = filtfilt(dcoOpticalFilter', rLed2);
dco_rLed3 = filtfilt(dcoOpticalFilter', rLed3);
dco_rLed4 = filtfilt(dcoOpticalFilter', rLed4);

dcoAx = filtfilt(dcoKinematicFilter', kAx);
dcoAy = filtfilt(dcoKinematicFilter', kAy);
dcoAz = filtfilt(dcoKinematicFilter', kAz);
dcoGx = filtfilt(dcoKinematicFilter', kGx);
dcoGy = filtfilt(dcoKinematicFilter', kGy);
dcoGz = filtfilt(dcoKinematicFilter', kGz);

dcoOdba = sqrt(dcoAx .^ 2 + dcoAy .^ 2 + dcoAz .^ 2);
dcoOdav = sqrt(dcoGx .^ 2 + dcoGy .^ 2 + dcoGz .^ 2);

dcoName = sprintf('HPF/DC-remove - order: %d | passbandFreq: %.2f Hz',...
    dcoFiltOrder, dcoPassbandFreqHz);
    
fprintf('DC-removed (high-pass above 0.01 Hz) data is prefixed as dco*\n');

%% plot the differences after DC removal

if (makePlot)
    
figDco = figure;

figDco.Name = dcoName;

s0 = subplot(5,2,[1,2]);
plot(oT_s, dcoLed1, 'Color', led1Color);
hold on;
plot(oT_s, dcoLed2, 'Color', led2Color);
plot(oT_s, dcoLed3, 'Color', led3Color);
plot(oT_s, dcoLed4, 'Color', led4Color);
xlabel('Time, seconds');
ylabel('Optical intensity');
title('Optics');
grid;
legend('Ambient1','1050 nm','1200 nm','Ambient2');


s1 = subplot(523);
plot(kT_s, dcoAx, 'Color', kxColor);
xlabel('Time, seconds');
ylabel('m/s^2');
title('Accel X-axis');
grid;
legend('dcoAx');


s2 = subplot(524);
plot(kT_s, dcoGx, 'Color', kxColor);
xlabel('Time, seconds');
ylabel('°/s');
title('Gyro X-axis');
grid;
legend('dcoGx');


s3 = subplot(525);
plot(kT_s, dcoAy, 'Color', kyColor);
xlabel('Time, seconds');
ylabel('m/s^2');
title('Accel Y-axis');
grid;
legend('dcoAy');

s4 = subplot(526);
plot(kT_s, dcoGy, 'Color', kyColor);
xlabel('Time, seconds');
ylabel('°/s');
title('Gyro Y-axis');
grid;
legend('dcoGy');

s5 = subplot(527);
plot(kT_s, dcoAz, 'Color', kzColor);
xlabel('Time, seconds');
ylabel('m/s^2');
title('Accel Z-axis');
grid;
legend('dcoAz');

s6 = subplot(528);
plot(kT_s, dcoGz, 'Color', kzColor);
xlabel('Time, seconds');
ylabel('°/s');
title('Gyro Z-axis');
grid;
legend('dcoGz');

s7 = subplot(529);
plot(kT_s, dcoOdba, 'k');
xlabel('Time, seconds');
ylabel('m/s^2');
title('dcoODBA');
grid;
legend('dcoOdba');

s8 = subplot(5,2,10);
plot(kT_s, dcoOdav, 'k');
xlabel('Time, seconds');
ylabel('°/s');
title('dcoODAV');
grid;
legend('dcoOdav');

linkaxes([s0 s1 s2 s3 s4 s5 s6 s7 s8],'x');

end         % end if (plotIt == true)

%% next up, a low-pass filter for removal of high frequency signals

lpFiltOrder       = 6;        
oPassbandFreqHz = oFs / 25;
kPassbandFreqHz = kFs / 4;

fprintf('Low-pass filter for removal of high frequency noise signals.\n');
fprintf('\tFilter order: %d\n', lpFiltOrder);
fprintf('\tOptics passband frequency: %d Hz\n', oPassbandFreqHz);
fprintf('\tKinematics passband frequency: %d Hz\n', kPassbandFreqHz);

% filter coefficients for high-pass

[Aolp, Bolp, Colp, Dolp] = butter(lpFiltOrder, ( oPassbandFreqHz / (oNyquistFs) ) );
[Aklp, Bklp, Cklp, Dklp] = butter(lpFiltOrder, ( kPassbandFreqHz / (kNyquistFs) ) );

% design the high-pass filter for optics
lowPassOpticalFilter = designfilt('lowpassiir','FilterOrder', lpFiltOrder, ...
    'PassbandFrequency', oPassbandFreqHz, ... 
    'PassbandRipple', 0.2, ...
    'SampleRate', oNyquistFs);

% design the high-pass filter for kinematics
lowPassKinematicFilter = designfilt('lowpassiir','FilterOrder', lpFiltOrder, ...
    'PassbandFrequency', kPassbandFreqHz, ... 
    'PassbandRipple', 0.2, ...
    'SampleRate', kNyquistFs);

if (makePlot)
    
    sos = ss2sos(Aolp,Bolp,Colp,Dolp);
    fvt = fvtool(sos, lowPassOpticalFilter, 'Fs', oNyquistFs);
    legend(fvt,'butter','designfilt');

    sos = ss2sos(Aklp, Bklp, Cklp, Dklp);
    fvt = fvtool(sos, lowPassKinematicFilter, 'Fs', kNyquistFs);
    legend(fvt,'butter','designfilt');
    
end


lpLed1 = filtfilt(lowPassOpticalFilter', dcoLed1);
lpLed2 = filtfilt(lowPassOpticalFilter', dcoLed2);
lpLed3 = filtfilt(lowPassOpticalFilter', dcoLed3);
lpLed4 = filtfilt(lowPassOpticalFilter', dcoLed4);

lpdLed1 = filtfilt(lowPassOpticalFilter', dcodLed1);
lpdLed2 = filtfilt(lowPassOpticalFilter', dcodLed2);
lpdLed3 = filtfilt(lowPassOpticalFilter', dcodLed3);
lpdLed4 = filtfilt(lowPassOpticalFilter', dcodLed4);

lpAx = filtfilt(lowPassKinematicFilter', dcoAx);
lpAy = filtfilt(lowPassKinematicFilter', dcoAy);
lpAz = filtfilt(lowPassKinematicFilter', dcoAz);
lpGx = filtfilt(lowPassKinematicFilter', dcoGx);
lpGy = filtfilt(lowPassKinematicFilter', dcoGy);
lpGz = filtfilt(lowPassKinematicFilter', dcoGz);

lpOdba = sqrt( lpAx .^ 2 + lpAy .^ 2 + lpAz .^ 2);
lpOdav = sqrt( lpGx .^ 2 + lpGy .^ 2 + lpGz .^ 2);
        


%% plot the differences after low-pass (high-frequency) removal

if (makePlot)
    
figLowPass = figure;

lpFigName = sprintf('Low-pass filter. oPassband: %d Hz, kPassband: %d Hz', ...
    oPassbandFreqHz, kPassbandFreqHz);

figLowPass.Name = lpFigName;

lp1 = subplot(5,2,[1,2]);
plot(oT, lpLed1,'k');
hold on;
plot(oT, lpLed2,'b');
plot(oT, lpLed3,'r');
plot(oT, lpLed4,'g');
xlabel('Time, local');
ylabel('Optical intensity');
title('Low-pass (high-frequency removal) optical signals');
grid;
legend('Ambient1','1050 nm','1200 nm','Ambient2');

lp2 = subplot(523);
plot(kT, dcoAx, 'b');
hold on;
plot(kT, lpAx, 'Color', kxColor);
hold off;
xlabel('Time, seconds');
ylabel('m/s^2');
title('Accel X-axis');
grid;
legend('dcoAx','lpAx');

lp3 = subplot(524);
plot(kT, dcoGx, 'b');
hold on;
plot(kT, lpGx, 'Color', kxColor);
hold off;
xlabel('Time, seconds');
ylabel('°/s');
title('Gyro X-axis');
grid;
legend('dcoGx','lpGx');

lp4 = subplot(525);
plot(kT, dcoAy, 'b');
hold on;
plot(kT, lpAy, 'Color', kyColor);
hold off;
xlabel('Time, seconds');
ylabel('m/s^2');
title('Accel Y-axis');
grid;
legend('dcoAy','lpAy');

lp5 = subplot(526);
plot(kT, dcoGy, 'b');
hold on;
plot(kT, lpGy, 'Color', kyColor);
hold off;
xlabel('Time, seconds');
ylabel('°/s');
title('Gyro Y-axis');
grid;
legend('dcoGy','lpGy');

lp6 = subplot(527);
plot(kT, dcoAz, 'b');
hold on;
plot(kT, lpAz, 'Color', kzColor);
hold off;
xlabel('Time, seconds');
ylabel('m/s^2');
title('Accel Z-axis');
grid;
legend('dcoAz','lpAz');

lp7 = subplot(528);
plot(kT, dcoGz, 'b');
hold on;
plot(kT, lpGz, 'Color', kzColor);
hold off;
xlabel('Time, seconds');
ylabel('°/s');
title('Gyro Z-axis');
grid;
legend('dcoGz','lpGz');

lp8 = subplot(529);
plot(kT, dcoOdba, 'k');
xlabel('Time, seconds');
ylabel('m/s^2');
title('Low-pass ODBA');
grid;
legend('lpOdba');

lp9 = subplot(5,2,10);
plot(kT, dcoOdav, 'k');
xlabel('Time, seconds');
ylabel('°/s');
title('Low-pass ODAV');
grid;
legend('lpOdav');

linkaxes([lp1 lp2 lp3 lp4 lp5 lp6 lp7 lp8 lp9],'x');

end         % end if (plotIt == true)

%% Try using a 30th order, frameSize frame Savitzky-Golay filter

%       Supposedly helpful for separating slow-varying motion wander


% add a logic test to shrink SG frame size for small number of samples

if (numel(kAx) <= 401)
    frameSize = 51;
else 
   frameSize = 401; 
end


fprintf('Savitzky-Golay (SG) filtering, 30th order, 401-frame.\n');
fprintf('\tsgoLed* = SG filter of raw OPTICS.led*\n');
fprintf('\tsgoA*/G* = SG filter of raw KINEMATICS.*\n');
fprintf('\tsgfLed* = SG filter of filtered Led* (-DC +LP optics)\n');
fprintf('\tsgfA*/G* = SG filter of filtered A*/G* (-DC +LP kinematics)\n');
fprintf('\tsgr* = SG filter subtraction from raw signals\n');

sgoLed1 = sgolayfilt(led1, 30, frameSize);
sgoLed2 = sgolayfilt(led2, 30, frameSize);
sgoLed3 = sgolayfilt(led3, 30, frameSize);
sgoLed4 = sgolayfilt(led4, 30, frameSize);

sgodLed1 = sgolayfilt(denoiseLed1, 30, frameSize);
sgodLed2 = sgolayfilt(denoiseLed2, 30, frameSize);
sgodLed3 = sgolayfilt(denoiseLed3, 30, frameSize);
sgodLed4 = sgolayfilt(denoiseLed4, 30, frameSize);

% high-pass (DC-removed) + low-pass optical signals with an SG filter

sgfLed1 = sgolayfilt(lpLed1, 30, frameSize);
sgfLed2 = sgolayfilt(lpLed2, 30, frameSize);
sgfLed3 = sgolayfilt(lpLed3, 30, frameSize);
sgfLed4 = sgolayfilt(lpLed4, 30, frameSize);

sgfdLed1 = sgolayfilt(lpdLed1, 30, frameSize);
sgfdLed2 = sgolayfilt(lpdLed2, 30, frameSize);
sgfdLed3 = sgolayfilt(lpdLed3, 30, frameSize);
sgfdLed4 = sgolayfilt(lpdLed4, 30, frameSize);

% SG-removed optics

sgrLed1 = led1 - sgoLed1;
sgrLed2 = led2 - sgoLed2;
sgrLed3 = led3 - sgoLed3;
sgrLed4 = led4 - sgoLed4;

sgfrLed1 = lpLed1 - sgfLed1;
sgfrLed2 = lpLed2 - sgfLed2;
sgfrLed3 = lpLed3 - sgfLed3;
sgfrLed4 = lpLed4 - sgfLed4;

sgrdLed1 = denoiseLed1 - sgodLed1;
sgrdLed2 = denoiseLed2 - sgodLed2;
sgrdLed3 = denoiseLed3 - sgodLed3;
sgrdLed4 = denoiseLed4 - sgodLed4;

% now do the kinematics

sgoAx = sgolayfilt(kAx, 30, frameSize);
sgoAy = sgolayfilt(kAy, 30, frameSize);
sgoAz = sgolayfilt(kAz, 30, frameSize);
sgoGx = sgolayfilt(kGx, 30, frameSize);
sgoGy = sgolayfilt(kGy, 30, frameSize);
sgoGz = sgolayfilt(kGz, 30, frameSize);

sgfAx = sgolayfilt(lpAx, 30, frameSize);
sgfAy = sgolayfilt(lpAy, 30, frameSize);
sgfAz = sgolayfilt(lpAz, 30, frameSize);
sgfGx = sgolayfilt(lpGx, 30, frameSize);
sgfGy = sgolayfilt(lpGy, 30, frameSize);
sgfGz = sgolayfilt(lpGz, 30, frameSize);

sgrAx = kAx - sgoAx;
sgrAy = kAy - sgoAy;
sgrAz = kAz - sgoAz;
sgrGx = kGx - sgoGx;
sgrGy = kGy - sgoGy;
sgrGz = kGz - sgoGz;

dcoOdba = sqrt(dcoAx.^2 + dcoAy.^2 + dcoAz.^2);
lpOdba  = sqrt(lpAx.^2  + lpAy.^2  + lpAz.^2 );
sgrOdba = sqrt(sgrAx.^2 + sgrAy.^2 + sgrAz.^2);

odav = sqrt( (kGx).^2 + (kGy).^2 + (kGz).^2 );
dcoOdav = sqrt(dcoGx.^2 + dcoGy.^2 + dcoGz.^2);
lpOdav  = sqrt(lpGx.^2  + lpGy.^2  + lpGz.^2 );
sgrOdav = sqrt(sgrGx.^2 + sgrGy.^2 + sgrGz.^2);


%% look at raw signals with subtracted Savitzky-Golay 

if (makePlot)
figSgolay = figure;
figSgolay.Name = 'Sovitzy-Golay filter removal from RAW signals';

plotSgrLed24 = subplot(421);
plot(oT_s, sgrLed2, 'Color', led2Color, 'LineWidth', 1.5);
hold on;
plot(oT_s, sgrLed4, 'Color', led4Color, 'LineWidth', 1.5);
hold off;
xlabel('Time, seconds');
ylabel('Optical intensity');
title('Savitzy-Golay removal from raw optical signals');
grid;
legend('1050 nm','Ambient2');

plotSgrLed31 = subplot(422);
plot(oT_s, sgrLed3, 'Color', led3Color, 'LineWidth', 1.5);
hold on;
plot(oT_s, sgrLed1, 'Color', led1Color, 'LineWidth', 1.5);
hold off;
xlabel('Time, seconds');
ylabel('Optical intensity');
title('Savitzy-Golay removal from raw optical signals');
grid;
legend('1200 nm','Ambient1');

plotSgrAx = subplot(423);
plot(kT_s, sgrAx, 'Color', kxColor, 'LineWidth', 1.2);
xlabel('Time, seconds');
ylabel('m/s^2');
title('Accel X-axis');
grid;
legend('sgrAx');

plotSgrGx = subplot(424);
plot(kT_s, sgrGx, 'Color', kxColor, 'LineWidth', 1.2);
xlabel('Time, seconds');
ylabel('°/s');
title('Gyro X-axis');
grid;
legend('sgrGx');

plotSgrAy = subplot(425);
plot(kT_s, sgrAy, 'Color', kyColor, 'LineWidth', 1.2);
xlabel('Time, seconds');
ylabel('m/s^2');
title('Accel Z-axis');
legend('sgrAz');
grid;

plotSgrGy = subplot(426);
plot(kT_s, sgrGy, 'Color', kyColor, 'LineWidth', 1.2);
xlabel('Time, seconds');
ylabel('°/s');
title('Gyro Y-axis');
legend('sgrGy');
grid;

plotSgrAz = subplot(427);
plot(kT_s, sgrAz, 'Color', kzColor, 'LineWidth', 1.2);
xlabel('Time, seconds');
ylabel('m/s^2');
title('Accel Z-axis');
legend('sgrAz');
grid;

plotSgrGz = subplot(428);
plot(kT_s, sgrGz, 'Color', kzColor, 'LineWidth', 1.2);
xlabel('Time, seconds');
ylabel('°/s');
title('Gyro Z-axis');
legend('sgrGz');
grid;

linkaxes([plotSgrLed24 plotSgrLed31 plotSgrAx plotSgrGx...
    plotSgrAy plotSgrGy plotSgrAz plotSgrGz],'x');

end         % end if (plotIt == true)

%% plot the differences after low-pass (high-frequency) removal

if (makePlot)
    
figDifferences = figure;
figDifferences.Name = 'Filter effects: visualizing the differences';

plot(oT_s, sgfLed1,'k');
hold on;
plot(oT_s, sgfLed2,'b');
plot(oT_s, sgfLed3,'r');
plot(oT_s, sgfLed4,'g');
xlabel('Time, seconds');
ylabel('Optical intensity');
title('-DC +LP +SG filtered (smoothed) optical signals');
legend('Ambient1','1050 nm','1200 nm','Ambient2');
grid;

figure;

s3a = subplot(321);
plot(kT_s, dcoAx, 'b:');
hold on;
plot(kT_s, sgrAx, 'r');
hold off;
xlabel('Time, seconds');
ylabel('m/s^2');
title('Accel X-axis');
legend('dcoAx','sgrAx');
grid;

s3b = subplot(322);
plot(kT_s, dcoGx, 'b:');
hold on;
plot(kT_s, sgrGx, 'r');
hold off;
xlabel('Time, seconds');
ylabel('°/s');
title('Gyro X-axis');
legend('dcoGx','sgrGx');
grid;

s3c = subplot(323);
plot(kT_s, dcoAy, 'b:');
hold on;
plot(kT_s, sgrAy, 'r');
hold off;
xlabel('Time, seconds');
ylabel('m/s^2');
title('Accel Y-axis');
legend('dcoAy','sgrAy');
grid;

s3d = subplot(324);
plot(kT_s, dcoGy, 'b:');
hold on;
plot(kT_s, sgrGy, 'r');
hold off;
xlabel('Time, seconds');
ylabel('°/s');
title('Gyro Y-axis');
legend('dcoGy','sgrGy');
grid;

s3e = subplot(325);
plot(kT_s, dcoAz, 'b:');
hold on;
plot(kT_s, sgrAz, 'r');
hold off;
xlabel('Time, seconds');
ylabel('m/s^2');
title('Accel Z-axis');
legend('dcoAz','sgrAz');
grid;

s3f = subplot(326);
plot(kT_s, dcoGz, 'b:');
hold on;
plot(kT_s, sgrGz, 'r');
hold off;
xlabel('Time, seconds');
ylabel('°/s');
title('Gyro Z-axis');
legend('dcoGz','sgrGz');
grid;

linkaxes([s3a s3b s3c s3d s3e s3f],'x');

end         % end if (plotIt == true)

%% Band-pass optical and kinematic signals

bandpassFiltOrder = 6;          % 6th order butterworth filter

oLowCut            = 0.3333;    % 0.3333 Hz = ~20 beats per minute
oHighCut           = 8.5;      % 3.5 Hz    = 210 beats per minute
opticalTestFs      = 25;        % 8 produces big pulses around odba HR in bpLed2                      


kLowCut             = 5.625;    % lower end for maximizing QRS ECG energy
kHighCut            = 22.5;     % higher end for maximizing QRS ECG energy

fprintf('Band-pass butterworth filter for optical cardio energies:\n');
fprintf('\tFilter order: %d\n', bandpassFiltOrder);
fprintf('\tLow cut frequency: %2.2f Hz\n', oLowCut);
fprintf('\tHigh cut frequency: %2.2f Hz\n', oHighCut);
fprintf('\tOptical BPF sample rate (2Hz = normalized frequencies): %d Hz\n', ...
    opticalTestFs);

[Aobp, Bobp, Cobp, Dobp] = butter(bandpassFiltOrder/2, ...
    [oLowCut oHighCut] / (opticalTestFs) );

[sosBpo,gBpo] = ss2sos(Aobp,Bobp,Cobp,Dobp);

bandPassHrOpticalFilter = designfilt('bandpassiir','FilterOrder',bandpassFiltOrder, ...
    'HalfPowerFrequency1',oLowCut,'HalfPowerFrequency2', oHighCut, ...
	'SampleRate', opticalTestFs);

[Akbp, Bkbp, Ckbp, Dkbp] = butter(bandpassFiltOrder/2, ...
    [kLowCut kHighCut] / (kNyquistFs) );

[sosBpk,gBpk] = ss2sos(Akbp,Bkbp,Ckbp,Dkbp);
  
bandPassHrKinematicFilter = designfilt('bandpassiir','FilterOrder',bandpassFiltOrder, ...
    'HalfPowerFrequency1', kLowCut,'HalfPowerFrequency2', kHighCut, ...
	'SampleRate', kNyquistFs);

fprintf('Band-pass butterworth filter for seismocardiography energies:\n');
fprintf('\tFilter order: %d\n', bandpassFiltOrder);
fprintf('\tLow cut frequency: %d Hz\n', kLowCut);
fprintf('\tHigh cut frequency: %d Hz\n', kHighCut);


if (makePlot)
    
    fvtO = fvtool(sosBpo, bandPassHrOpticalFilter, 'Fs', opticalTestFs);
    legend(fvtO,'butter','designfilt');


    fvtK = fvtool(sosBpk, bandPassHrKinematicFilter, 'Fs', kNyquistFs);
    legend(fvtK,'butter','designfilt');

end

% band-pass raw optical and kinematic signals

bpLed1 = filtfilt(bandPassHrOpticalFilter', led1);
bpLed2 = filtfilt(bandPassHrOpticalFilter', led2);
bpLed3 = filtfilt(bandPassHrOpticalFilter', led3);
bpLed4 = filtfilt(bandPassHrOpticalFilter', led4);

bprLed1 = filtfilt(bandPassHrOpticalFilter', rLed1);
bprLed2 = filtfilt(bandPassHrOpticalFilter', rLed2);
bprLed3 = filtfilt(bandPassHrOpticalFilter', rLed3);
bprLed4 = filtfilt(bandPassHrOpticalFilter', rLed4);

bpDcoLed1 = filtfilt(bandPassHrOpticalFilter', dcoLed1);
bpDcoLed2 = filtfilt(bandPassHrOpticalFilter', dcoLed2);
bpDcoLed3 = filtfilt(bandPassHrOpticalFilter', dcoLed3);
bpDcoLed4 = filtfilt(bandPassHrOpticalFilter', dcoLed4);

bpLpLed1 = filtfilt(bandPassHrOpticalFilter', lpLed1);
bpLpLed2 = filtfilt(bandPassHrOpticalFilter', lpLed2);
bpLpLed3 = filtfilt(bandPassHrOpticalFilter', lpLed3);
bpLpLed4 = filtfilt(bandPassHrOpticalFilter', lpLed4);

bpSgrLed1 = filtfilt(bandPassHrOpticalFilter', sgrLed1);
bpSgrLed2 = filtfilt(bandPassHrOpticalFilter', sgrLed2);
bpSgrLed3 = filtfilt(bandPassHrOpticalFilter', sgrLed3);
bpSgrLed4 = filtfilt(bandPassHrOpticalFilter', sgrLed4);



bpAx   = filtfilt(bandPassHrKinematicFilter',kAx);
bpAy   = filtfilt(bandPassHrKinematicFilter',kAy);
bpAz   = filtfilt(bandPassHrKinematicFilter',kAz);
bpGx   = filtfilt(bandPassHrKinematicFilter',kGx);
bpGy   = filtfilt(bandPassHrKinematicFilter',kGy);
bpGz   = filtfilt(bandPassHrKinematicFilter',kGz);

bpOdba = sqrt( bpAx.^2 + bpAy.^2 + bpAz.^2);
bpOdav = sqrt( bpGx.^2 + bpGy.^2 + bpGz.^2);


%% Now all the heart rate investigation variables should be in place...

if (makePlot) 
    
figBpOpt = figure;
figBpTitle = sprintf('Band-passed optics | Order: %d | %2.1f Hz to %2.1f Hz, (no other filters)', ...
    bandpassFiltOrder, kLowCut, kHighCut);
figBpOpt.Name = figBpTitle;

pBpLed1 = subplot(221);
plot(oT, bpLed1,'Color', black);
xlabel('Time, seconds');
ylabel('Optical intensity');
titleTxt = sprintf('%s - sgr optics', tag);
title('Band-passed optics (no HP+LP+SG subtraction)');
grid;
legend('ambient1');

pBpLed2 = subplot(222);
plot(oT, bpLed2, 'Color', blue);
xlabel('Time, seconds');
ylabel('Optical intensity');
titleTxt = sprintf('%s - sgr optics', tag);
title('Band-passed optics (no HP+LP+SG subtraction)');
grid;
legend('1050 nm');

pBpLed3 = subplot(223);
plot(oT, bpLed3, 'Color', red);
xlabel('Time, local');
ylabel('Optical intensity');
titleTxt = sprintf('%s - sgr optics', tag);
title('Band-passed optics (no HP+LP+SG subtraction)');
grid;
legend('1200 nm');

pBpLed4 = subplot(224);
plot(oT, bpLed4,'Color', goldenrod);
xlabel('Time, local');
ylabel('Optical intensity');
titleTxt = sprintf('%s - sgr optics', tag);
title('Band-passed optics (no HP+LP+SG subtraction)');
grid;
legend('ambient2');

linkaxes([pBpLed1 pBpLed2 pBpLed3 pBpLed4],'x');


figBpSgrOpt = figure;
figBpSgrOpt.Name = 'Band-passed Sovitky-Golay removed optics';

pBpSgrLed1 = subplot(221);
plot(oT, bpSgrLed1,'Color', black);
xlabel('Time, seconds');
ylabel('Optical intensity');
titleTxt = sprintf('%s - sgr optics', tag);
title('Optics after HP+LP+SG subtraction & 6°BPBW Filter');
grid;
legend('ambient1');

pBpSgrLed2 = subplot(222);
plot(oT, bpSgrLed2, 'Color', blue);
xlabel('Time, seconds');
ylabel('Optical intensity');
titleTxt = sprintf('%s - sgr optics', tag);
title('Optics after HP+LP+SG subtraction & 6°BPBW Filter');
grid;
legend('1050 nm');

pBpSgrLed3 = subplot(223);
plot(oT, bpSgrLed3, 'Color', red);
xlabel('Time, local');
ylabel('Optical intensity');
titleTxt = sprintf('%s - sgr optics', tag);
title('Optics after HP+LP+SG subtraction & 6°BPBW Filter');
grid;
legend('1200 nm');

pBpSgrLed4 = subplot(224);
plot(oT, bpSgrLed4,'Color', goldenrod);
xlabel('Time, local');
ylabel('Optical intensity');
titleTxt = sprintf('%s - sgr optics', tag);
title('Optics after HP+LP+SG subtraction & 6°BPBW Filter');
grid;
legend('ambient2');

linkaxes([pBpSgrLed1 pBpSgrLed2 pBpSgrLed3 pBpSgrLed4],'x');

% Plot some band-passed kinematic data


figBpKin = figure('Position',[100 100 800 600]); 

s5a = subplot(211);
plot(kT_s, bpAx, kT_s, bpAy, kT_s, bpAz);
xlabel('Time, seconds');
ylabel('bpHrA*, m/s^2');
title('Accelerometry HR estimator');
grid;
legend('ax', 'ay', 'az');


s5b = subplot(212);

plot(kT_s, bpGx, kT_s, bpGy, kT_s, bpGz); grid;
xlabel('Time, seconds');
ylabel('bpHrG*, °/s');
title('Gyroscopy HR estimator');
grid;
legend('gx', 'gy', 'gz');

linkaxes([s5a s5b],'x');

% plot individual A* and G* windows to separate what may be inferrable

figBpKin = figure;
figBpKin.Name = 'Band-passed aggregated kinematics';

s6a = subplot(421);
plot(kT, bpAx, 'Color', blue);
xlabel('Time, local');
ylabel('m/s^2');
title('Accel X-axis');
grid;
legend('bpAx');

s6b = subplot(422);
plot(kT, bpGx, 'Color', blue);
xlabel('Time, local');
ylabel('°/s');
title('Gyro X-axis');
grid;
legend('bpGx');


s6c = subplot(423);
plot(kT, bpAy, 'Color', red);
xlabel('Time, local');
ylabel('m/s^2');
title('Accel Y-axis');
grid;
legend('bpAy');


s6d = subplot(424);
plot(kT, bpGy, 'Color', red);
xlabel('Time, local');
ylabel('°/s');
title('Gyro Y-axis');
grid;
legend('bpGy');


s6e = subplot(425);
plot(kT, bpAz, 'Color', goldenrod);
xlabel('Time, local');
ylabel('m/s^2');
title('Accel Z-axis');
grid;
legend('bpAz');

s6f = subplot(426);
plot(kT, bpGz, 'Color', goldenrod);
xlabel('Time, local');
ylabel('°/s');
title('Gyro Z-axis');
grid;
legend('bpGz');

s6g = subplot(427);
plot(kT, bpOdba, 'Color', green);
xlabel('Time, local');
ylabel('m/s^2');
title('ODBA');
grid;
legend('bpOdba');


s6h = subplot(428);
plot(kT, bpOdav, 'Color', green);
xlabel('Time, local');
ylabel('°/s');
title('ODAV');
grid;
legend('bpOdav');

linkaxes([s6a s6b s6c s6d s6e s6f s6g s6h],'x');

end         % end if (plotIt == true) for band-pass plots

%% Prompt to harshly close all open figures before proceeding...
% 
% closeFigsStr = lower(input('Close all open figures? (y/n): ','s'));
% switch(closeFigsStr)
%     case 'y'
%         fprintf('Closing open figures...\n');
%         figHandles = get(groot, 'Children');
%         if (~isempty(figHandles))
%             close(figHandles);
%         end
%     otherwise
%         fprintf('Leaving figures in place...\n');
% end

%% Okay, everything is in place for heart rate analysis...

%  Pause here and let the user select which analysis to run

fprintf('\n\n');
fprintf('Signal processing is complete on this epoch.\n');
fprintf('Press [enter] to continue analysis, or [q] to go manual.\n');
waitForInput = lower(input('Take your pick... (enter/q) [enter]: ','s'));
switch(waitForInput)
    case 'q'
        clc;
        fprintf('Bailing ftHeartRate_Kinematics.m around line 1390...\n');
        fprintf('Pick it up from there!\n');
        return;
	otherwise
        fprintf('Continuing with analysis...\n');
end
        

%% plot band-passed kinematics, so user knows what they have to work with

fprintf('Here are band-passed kinematics for reference...\n');

figKinematics = figure('Position',[50 50 1800 1200]);
figKinematics.Name = 'Band-passed aggregated kinematics';

s6a = subplot(421);
plot(kT, bpAx, 'Color', blue);
xlabel('Time, local');
ylabel('m/s^2');
title('Accel X-axis');
grid;
legend('bpAx');

s6b = subplot(422);
plot(kT, bpGx, 'Color', blue);
xlabel('Time, local');
ylabel('°/s');
title('Gyro X-axis');
grid;
legend('bpGx');

s6c = subplot(423);
plot(kT, bpAy, 'Color', red);
xlabel('Time, local');
ylabel('m/s^2');
title('Accel Y-axis');
grid;
legend('bpAy');

s6d = subplot(424);
plot(kT, bpGy, 'Color', red);
xlabel('Time, local');
ylabel('°/s');
title('Gyro Y-axis');
grid;
legend('bpGy');

s6e = subplot(425);
plot(kT, bpAz, 'Color', goldenrod);
xlabel('Time, local');
ylabel('m/s^2');
title('Accel Z-axis');
grid;
legend('bpAz');

s6f = subplot(426);
plot(kT, bpGz, 'Color', goldenrod);
xlabel('Time, local');
ylabel('°/s');
title('Gyro Z-axis');
grid;
legend('bpGz');

s6g = subplot(427);
plot(kT, bpOdba, 'Color', red);
xlabel('Time, local');
ylabel('m/s^2');
title('ODBA');
grid;
legend('bpOdba');

s6h = subplot(428);
plot(kT, bpOdav, 'Color', green);
xlabel('Time, local');
ylabel('°/s');
title('ODAV');
grid;
legend('bpOdav');

linkaxes([s6a s6b s6c s6d s6e s6f s6g s6h],'x');

%% Perform variational mode decomposition on bpOdba and bpOdav

%   Output from this section will produce candidate hrBpOdba & hrBpOdav

fprintf('Computing variational mode decompositions of bpOdba and bpOdav...\n');

tic

fprintf('Computing band-passed ODBA IMFs...\n');
[imfBpOdba, residualBpOdba, infoBpOdba] = vmd(bpOdba, 'NumIMF', 10);

fprintf('Computing band-passed ODAV IMFs...\n');
[imfBpOdav, residualBpOdav, infoBpOdav] = vmd(bpOdav, 'NumIMF', 10);

% fprintf('Computing HP-DC ODBA IMFs...\n');
% [imfDcoOdba, residualDcoOdba, infoDcoOdba] = vmd(dcoOdba, 'NumIMF', 10);
% 
% fprintf('Computing HP-DC ODAV IMFs...\n');
% [imfDcoOdav, residualDcoOdav, infoDcoOdav] = vmd(dcoOdav, 'NumIMF', 10);

fprintf('Completed variational mode decompositions for bpOdba and bpOdav...\n');

toc

%% Figure out how many IMFs were computed for bpOdba and bpOdav

numImfsBpOdba = width(imfBpOdba);
numImfsBpOdav = width(imfBpOdav);

% compute the sum of all but the first and last IMFs

sumImfBpOdba = sum(imfBpOdba(:,2:numImfsBpOdba-1),2);
sumImfBpOdav = sum(imfBpOdav(:,2:numImfsBpOdav-1),2);

% bpOdba VMD IMFs

figImfsOdba = figure('Position', [200 100 1000 1000]);

for i = 1:numImfsBpOdba
   
    pOdba(i) = subplot((round(numImfsBpOdba/2)),2,i);
    plot(kT, imfBpOdba(:,i), 'Color', red);
    xlabel('Time, local');
    ylabel('IMF');
    titleTxt = sprintf('BpOdba IMF %d', i);
    title(titleTxt);
    grid;    
    
end

linkaxes(pOdba(1:numImfsBpOdba),'x');

% bpOdav VMD IMFs
figImfsOdav = figure('Position',[1300 100 1000 1000]);

for i = 1:numImfsBpOdav
   
    pOdav(i) = subplot((round(numImfsBpOdav/2)),2,i);
    plot(kT, imfBpOdav(:,i), 'Color', green);
    xlabel('Time, local');
    ylabel('IMF');
    titleTxt = sprintf('BpOdav IMF %d', i);
    title(titleTxt);
    grid;    
    
end

linkaxes(pOdav(1:numImfsBpOdav),'x');


%% Set candidate for hrBpOdba and hrBpOdav

% allow the user to manually select which IMF works best 

bestImfOdba = 0;
bestImfOdav = 0;

while(bestImfOdba == 0)
    
    imfOdbaTxt = 'Select bpOdba IMF number that best characterizes HR signal: ';
    proposedImfOdba = str2double(input(imfOdbaTxt,'s'));
    if (proposedImfOdba >= 2 & proposedImfOdba <= numImfsBpOdba)
        fprintf('You chose IMF number %d. Proceeding...\n', proposedImfOdba);
        bestImfOdba = proposedImfOdba;
    else
        fprintf('Enter a valid IMF number.\n');
    end
    
end

while(bestImfOdav == 0)
    
    imfOdavTxt = 'Select bpOdav IMF number that best characterizes HR signal: ';
    proposedImfOdav = str2double(input(imfOdavTxt,'s'));
    if (proposedImfOdav >= 2 & proposedImfOdav <= numImfsBpOdav)
        fprintf('You chose IMF number %d. Proceeding...\n', proposedImfOdav);
        bestImfOdav = proposedImfOdav;
    else
        fprintf('Enter a valid IMF number.\n');
    end
    
end

% 1. If breath-hold, maybe IMF 9 always works well (as in tt21_142d BH) ?
% 2. If free-breathe and recover, maybe IMF 9 always works best as well?
% 3. After much review, IMF 9 always works best

bestImfBpOdba = imfBpOdba(:,bestImfOdba);
bestImfBpOdav = imfBpOdav(:,bestImfOdav);


%% hrBpOdba and hrBpOdav - the best VMD IMFs for HR sqwave & peak-finding

figure;

p1Odba = subplot(211);
plot(kT, bpOdba, 'Color', red);
hold on;
plot(kT, bestImfBpOdba, 'Color', blue);
hold off;
xlabel('Time, local');
ylabel('SCG');
title('bpOdba and bestImfBpOdba (bestImf for HR estimation)');
grid;
legend('bpOdba','bestImfBpOdba');

p2Odav = subplot(212);
plot(kT, bpOdav, 'Color', green);
hold on;
plot(kT, bestImfBpOdav, 'Color', blue);
hold off;
xlabel('Time, local');
ylabel('GCG');
title('bpOdav and bestImfBpOdav (bestImf for HR estimation)');
grid;
legend('bpOdav','bestImfBpOdav');

linkaxes([p1Odba p2Odav],'x');
    
%% bpOdav (for reference) + bpOdav IMFs - first and last + IMF 3 + IMF 4

figBestImfs = figure('Position',[100 100 1400 900]);

p1ImfOdba = subplot(211);
plot(kT, bpOdba, 'Color', red);
hold on;
plot(kT, smooth(sumImfBpOdba), 'Color', blue);
plot(kT, bestImfBpOdba, 'Color', maroon);
hold off;
xlabel('Time, local');
ylabel('bpOdba and sumImfBpOdba');
grid;
legend('bpOdba','sumImfBpOdba','hrBpOdba');

p2ImfOdav = subplot(212);
plot(kT, bpOdav, 'Color', green);
hold on;
plot(kT, smooth(sumImfBpOdav), 'Color', blue);
plot(kT, bestImfBpOdav, 'Color', purple);
hold off;
xlabel('Time, local');
ylabel('bpOdav');
grid;
legend('bpOdav','sumImfBpOdav','hrBpOdav');

linkaxes([p1ImfOdba p2ImfOdav],'x');

%% Use the best IMFs for bpOdba and bpOdav as candidate square wave targets

meanBestImfBpOdba = mean(bestImfBpOdba);
sdBestImfBpOdba = std(bestImfBpOdba);
meanBestImfBpOdav = mean(bestImfBpOdav);
sdBestImfBpOdav = std(bestImfBpOdav);

fprintf('Mean bestImfBpOdba: %d\n', meanBestImfBpOdba);
fprintf('SD bestImfBpOdba: %d\n', sdBestImfBpOdba);
fprintf('Mean bestImfBpOdbav): %d\n', meanBestImfBpOdav);
fprintf('SD bestImfBpOdav: %d\n', sdBestImfBpOdav);

thresholdBestImfOdba = meanBestImfBpOdba + (sdBestImfBpOdba * 0.5);
thresholdBestImfOdav = meanBestImfBpOdav + (sdBestImfBpOdav * 0.5);

targetBestImfBpOdba = (bestImfBpOdba <= thresholdBestImfOdba);
targetBestImfBpOdav = (bestImfBpOdav <= thresholdBestImfOdav);

%% compute proposedHr and proposedAntiHr ODBA and ODAV signals

fprintf('Creating best IMF bpOdba-based AutC square wave...\n');
threshSqwaveBestImfBpOdba = double(~bwareaopen(targetBestImfBpOdba, (kFs/5)));
threshSqwaveAntiBestImfBpOdba = double(bwareaopen(targetBestImfBpOdba, (kFs/5)));

fprintf('Creating best IMF bpOdav-based AutC square wave...\n');
threshSqwaveBestImfBpOdav = double(~bwareaopen(targetBestImfBpOdav, (kFs/5)));
threshSqwaveAntiBestImfBpOdav = double(bwareaopen(targetBestImfBpOdav, (kFs/5)));


%% plot hrBpOdba and hrBpOdav with their proposed Hr and AntiHr square waves

figSqwaves = figure('Position',[100 300 1400 800]);

p1SqwaveOdba = subplot(211);
plot(kT, bestImfBpOdba, 'Color', red);
hold on;
plot(kT, threshSqwaveBestImfBpOdba, '-', 'Color', blue);
% plot(kT, threshSqwaveAntiBestImfBpOdba, '-', 'Color', goldenrod);
yline(thresholdBestImfOdba, '--', 'Color', red);
hold off;
xlabel('Time, local');
ylabel('BestIMF BpOdba');
grid;
legend('bestImfBpOdba','sqwaveBestImf','threshold');

p1SqwaveOdav = subplot(212);
plot(kT, bestImfBpOdav, 'Color', green);
hold on;
plot(kT, threshSqwaveBestImfBpOdav, '-', 'Color', blue);
% plot(kT, threshSqwaveAntiBestImfBpOdav, '-', 'Color', goldenrod);
yline(thresholdBestImfOdav,'--','Color', green);
hold off;
xlabel('Time, local');
ylabel('BestIMF BpOdav');
grid;
legend('bestImfBpOdba','sqwaveBestImf','threshold');

linkaxes([p1SqwaveOdba p1SqwaveOdav],'x');

%% Bring in the peak finders now...

%   Do findpeaks on hrBpOdba and hrBpOdav.
%   Use these to try calculating HR frequency -and-
%   respiratory frequency

[pkBestImfBpOdba, locBestImfBpOdba, widBestImfBpOdba, promBestImfBpOdba] = ...
    findpeaks(bestImfBpOdba, 'MinPeakHeight', thresholdBestImfOdba, ...
    'MinPeakDistance', 50);
    
[pkBestImfBpOdav, locBestImfBpOdav, widBestImfBpOdav, promBestImfBpOdav] = ...
    findpeaks(bestImfBpOdav, 'MinPeakHeight', thresholdBestImfOdav, ...
    'MinPeakDistance', 50);

% create timetables from these findpeak results

HR_ODBA = timetable(kT(locBestImfBpOdba), pkBestImfBpOdba, widBestImfBpOdba, ... 
    promBestImfBpOdba, 'VariableNames', {'peaks','width','prominence'});

HR_ODAV = timetable(kT(locBestImfBpOdav), pkBestImfBpOdav, widBestImfBpOdav, ...
    promBestImfBpOdav, 'VariableNames', {'peaks','width','prominence'});

% now spline interpolate in context of kT

bestImfBpOdbaSpline = interp1(HR_ODBA.Time, HR_ODBA.peaks, kT, 'spline');
bestImfBpOdbaRespCycle = csaps(kT_s(locBestImfBpOdba), pkBestImfBpOdba, 0.5, kT_s);

bestImfBpOdavSpline = interp1(HR_ODAV.Time, HR_ODAV.peaks, kT, 'spline');
bestImfBpOdavRespCycle = csaps(kT_s(locBestImfBpOdba), pkBestImfBpOdba, 0.5, kT_s);

% generate CUEs within the range of kT

cueTimes = CUE.Time(CUE.Time < kT(end) & CUE.Time > kT(1));
cueTypes = CUE.type(unique(cueTimes));

f1 = figure;
f1.Name = 'Respiration test based on IMF4 BpOdba';
figSize = f1.Position;
figSize(3) = screenSize(3);
figSize(4) = screenSize(4) * 0.5;
f1.Position = figSize;  

pPeaksOdba = subplot(211);
plot(kT, bestImfBpOdba, 'Color', red, 'LineWidth', 1.1);
hold on;
yline(thresholdBestImfOdba, '--', 'Color', red);
plot(kT(locBestImfBpOdba), pkBestImfBpOdba, 'd', 'Color', blue, 'MarkerFace', blue);
plot(kT, bestImfBpOdbaSpline, 'Color', blue);
plot(kT, bestImfBpOdbaRespCycle, 'Color', maroon, 'LineWidth', 1.1);
yPos = pPeaksOdba.YLim(2) - (pPeaksOdba.YLim(2) * 0.1);
numCuePlots = ftPlotCues(cueTimes, cueTypes, yPos);
hold off;
xlabel('Time, local');
ylabel('hrBpOdba');
grid;
legend('bestImfBpOdba','threshold','Peaks','peakSpline','respCycle');

pPeaksOdav = subplot(212);
plot(kT, bestImfBpOdav, 'Color', green, 'LineWidth', 1.1);
hold on;
yline(thresholdBestImfOdav, '--', 'Color', green);
plot(kT(locBestImfBpOdav), pkBestImfBpOdav, 'd', 'Color', blue, 'MarkerFace', blue);
plot(kT, bestImfBpOdavSpline, 'Color', green);
plot(kT, bestImfBpOdavRespCycle, 'Color', maroon, 'LineWidth', 1.1);
yPos = pPeaksOdav.YLim(2) - (pPeaksOdav.YLim(2) * 0.1);
numCuePlots = ftPlotCues(cueTimes, cueTypes, yPos);
hold off;
xlabel('Time, local');
ylabel('HrBpOdba');
grid;
legend('bestImfBpOdav','threshold','Peaks','peakSpline','respCycle');

linkaxes([pPeaksOdba pPeaksOdav],'x');

%% Try to make some frequency conversions of these peaks

bestImfBpOdbaLength = numel(bestImfBpOdba);
bestImfBpOdavLength = numel(bestImfBpOdav);

peakBestImfBpOdbaElements = numel(locBestImfBpOdba);
peakBestImfBpOdavElements = numel(locBestImfBpOdav);

peaksBestImfBpOdba = zeros(bestImfBpOdbaLength,1);
peaksBestImfBpOdav = zeros(bestImfBpOdavLength,1);

for i = 1:peakBestImfBpOdbaElements

    thisStart = locBestImfBpOdba(i);
    thisEnd = locBestImfBpOdba(i) + round( widBestImfBpOdba(i));
    if (thisEnd > bestImfBpOdbaLength)
        peaksBestImfBpOdba(thisStart:bestImfBpOdbaLength,1) = 1;
    else
        peaksBestImfBpOdba(thisStart:thisEnd,1) = 1;
    end
    

end

for i = 1:peakBestImfBpOdavElements
   
    thisStart = locBestImfBpOdav(i);
    thisEnd = locBestImfBpOdav(i) + round( widBestImfBpOdav(i));
    if (thisEnd > bestImfBpOdavLength)
        peaksBestImfBpOdav(thisStart:bestImfBpOdavLength,1) = 1;    
    else
        peaksBestImfBpOdav(thisStart:thisEnd,1) = 1;    
    end
    
    
end


figPeakFreq = figure('Position',[100 500 1400 1000]);

pPeak1 = subplot(211);
plot(kT, peaksBestImfBpOdba, 'Color', red);
hold off;
xlabel('Time, seconds');
ylabel('HR Pulses');
title('BestImfBpOdba - binary peak heart rate estimate');

pPeak2 = subplot(212);
plot(kT, peaksBestImfBpOdav, 'Color', green);
hold off;
xlabel('Time, seconds');
ylabel('HR pulses');
title('BestImfBpOdav - binary peak heart rate estimate');

linkaxes([pPeak1 pPeak2],'x');


%% compute WSSTs

fprintf('Starting WSST computations...\n');

tic

fprintf('Creating WSSTs for band-passed ODBA (bpOdba) and ODAV (bpOdav)...\n');
[wsstBpOdba, fBpOdba] = wsst(bpOdba, kFs);
[wsstBpOdav, fBpOdav] = wsst(bpOdav, kFs);

fprintf('Creating WSSTs for sumImfBpOd** signals...\n');
[wsstSumImfBpOdba, fSumImfBpOdba] = wsst(sumImfBpOdba, kFs);
[wsstSumImfBpOdav, fSumImfBpOdav] = wsst(sumImfBpOdav, kFs);

fprintf('Creating WSSTs for bestImf (%d) BpOdba & bestImf (%d) bpOdav...\n', ...
    bestImfOdba, bestImfOdav);
[wsstBestImfBpOdba, fBestImfBpOdba] = wsst(bestImfBpOdba, kFs);
[wsstBestImfBpOdav, fBestImfBpOdav] = wsst(bestImfBpOdav, kFs);

fprintf('Creating WSSTs for threshold-based sqwaves of bestImfs bpOdba & bpOdav)...\n');
[wsstThreshSqwaveBestImfBpOdba, fThreshSqwaveBestImfBpOdba] = ...
    wsst(threshSqwaveBestImfBpOdba, kFs);
[wsstThreshSqwaveBestImfBpOdav, fThreshSqwaveBestImfBpOdav] = ...
    wsst(threshSqwaveBestImfBpOdav, kFs);

fprintf('Creating WSSTs for peak-based sqwaves of bestImfs of bpOdba & bpOdav...\n');
[wsstPeakSqwaveBestImfBpOdba, fPeakSqwaveBestImfBpOdba] = ...
    wsst(peaksBestImfBpOdba, kFs);
[wsstPeakSqwaveBestImfBpOdav, fPeakSqwaveBestImfBpOdav] = ...
    wsst(peaksBestImfBpOdav, kFs);

fprintf('Completed WSST computations for all signals of interest...\n');
toc;

%% Plot all the bpOdba WSSTs

tic

fprintf('Plotting all the WSSTs...BE PATIENT! This can take a minute!\n');

figWsstPlots = figure('Position',[50 50 800 1200 ]);

p1wsst1 = subplot(521);
pcolor(kT, fBpOdba, abs(wsstBpOdba));
shading interp;
ylim([0 4]);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('WSST for bpOdba');

p1wsst3 = subplot(523);
pcolor(kT, fSumImfBpOdba, abs(wsstSumImfBpOdba));
shading interp;
ylim([0 4]);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('WSST for SumImfBpOdba');

p1wsst5 = subplot(525);
pcolor(kT, fBestImfBpOdba, abs(wsstBestImfBpOdba));
shading interp;
ylim([0 4]);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('WSST for bestImfBpOdba');

p1wsst7 = subplot(527);
pcolor(kT, fThreshSqwaveBestImfBpOdba, abs(wsstThreshSqwaveBestImfBpOdba));
shading interp;
ylim([0 4]);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('WSST for threshold-based square-wave of bestImfBpOdba');

p1wsst9 = subplot(529);
pcolor(kT, fPeakSqwaveBestImfBpOdba, abs(wsstPeakSqwaveBestImfBpOdba));
shading interp;
ylim([0 4]);
xlabel('Time, seconds');
ylabel('Frequency, Hz');
title('WSST for peak-based square-wave of bestImfBpOdba');

p1wsst2 = subplot(522);
pcolor(kT, fBpOdav, abs(wsstBpOdav));
shading interp;
ylim([0 4]);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('WSST for bpOdav');

p1wsst4 = subplot(524);
pcolor(kT, fSumImfBpOdav, abs(wsstSumImfBpOdav));
shading interp;
ylim([0 4]);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('WSST for SumImfBpOdav');

p1wsst6 = subplot(526);
pcolor(kT, fBestImfBpOdav, abs(wsstBestImfBpOdav));
shading interp;
ylim([0 4]);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('WSST for bestImfBpOdav');

p1wsst8 = subplot(528);
pcolor(kT, fThreshSqwaveBestImfBpOdav, abs(wsstThreshSqwaveBestImfBpOdav));
shading interp;
ylim([0 4]);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('WSST for threshold-based square-wave of bestImfBpOdba');

p1wsst10 = subplot(5,2,10);
pcolor(kT, fPeakSqwaveBestImfBpOdav, abs(wsstPeakSqwaveBestImfBpOdav));
shading interp;
ylim([0 4]);
xlabel('Time, seconds');
ylabel('Frequency, Hz');
title('WSST for peak-based square-wave of bestImfBpOdav');

fprintf('Completed plots of WSSTs for hrOdba and hrOdav.\n');

toc


%% before computing wsstridges, allow user to exclude low frequencies

useEntireSpectrumTxt = 'Use entire frequency spectrum? (y/n) [y]: ';

useEntireSpectrum = lower(input(useEntireSpectrumTxt,'s'));

switch(useEntireSpectrum)
    
    case 'n'
        
        % Sort to find the min and max frequencies buckets representing
        % >= ~20 BPM (0.3 to 0.35 Hz) and <= ~240 BPM (3.9 to 4.1 Hz);
        % use the respective min and max, to preserve max frequency space.
        % This will prevent super-low frequency hopping from confusion 
        % with presumed Mayer waves / respiratory signals at ~0.1 to 0.2 Hz
        
        lowF1 = 0.3;    % 18 BPM
        lowF2 = 0.35;   % 21 BPM
        highF1 = 3.9;   % 234 BPM
        highF2 = 4.1;   % 246 BPM
        
        % do fBpOdba
        minF_BpOdba = min(find(fBpOdba > lowF1 & fBpOdba < lowF2));
        maxF_BpOdba = max(find(fBpOdba > highF1 & fBpOdba < highF2));
        
        % do fBpOdav
        minF_BpOdav = min(find(fBpOdav > lowF1 & fBpOdav < lowF2));
        maxF_BpOdav = max(find(fBpOdav > highF1 & fBpOdav < highF2));        
        
        % do fHrBpOdba
        minF_HrBpOdba = min(find(fBestImfBpOdba > lowF1 & fBestImfBpOdba < lowF2));
        maxF_HrBpOdba = max(find(fBestImfBpOdba > highF1 & fBestImfBpOdba < highF2));
        
        % do fHrBpOdav
        minF_HrBpOdav = min(find(fBestImfBpOdav > lowF1 & fBestImfBpOdav < lowF2));
        maxF_HrBpOdav = max(find(fBestImfBpOdav > highF1 & fBestImfBpOdav < highF2));        
        
        % do fSqwaveOdba
        minF_SqwaveOdba = min(find(fThreshSqwaveBestImfBpOdba > lowF1 & fThreshSqwaveBestImfBpOdba < lowF2));
        maxF_SqwaveOdba = max(find(fThreshSqwaveBestImfBpOdba > highF1 & fThreshSqwaveBestImfBpOdba < highF2));
        
        % do fSqwaveOdav
        minF_SqwaveOdav = min(find(fThreshSqwaveBestImfBpOdav > lowF1 & fThreshSqwaveBestImfBpOdav < lowF2));
        maxF_SqwaveOdav = max(find(fThreshSqwaveBestImfBpOdav > highF1 & fThreshSqwaveBestImfBpOdav < highF2));
        
        % ensure that all newF values have the same min and max
        minTestSum = minF_BpOdba + minF_BpOdav + minF_HrBpOdba + ...
            minF_HrBpOdav + minF_SqwaveOdba + minF_SqwaveOdav;
        
        maxTestSum = maxF_BpOdba + maxF_BpOdav + maxF_HrBpOdba + ...
            maxF_HrBpOdav + maxF_SqwaveOdba + maxF_SqwaveOdav;
        
        if ( ((minTestSum / 6) == minF_BpOdba) & ((maxTestSum /6) == maxF_BpOdba) )
            
            fprintf('Sizes of new Fs match. Reducing before ridge finding.\n');

            % recast the Fs from minF to maxF
            
            fBpOdba_        = fBpOdba(1,minF_BpOdba:maxF_BpOdba);
            fBpOdav_        = fBpOdav(1,minF_BpOdav:maxF_BpOdav);
            fHrBpOdba_      = fBestImfBpOdba(1,minF_HrBpOdba:maxF_HrBpOdba);
            fHrBpOdav_      = fBestImfBpOdav(1,minF_HrBpOdav:maxF_HrBpOdav);
            fSqwaveOdba_    = fThreshSqwaveBestImfBpOdba(1,minF_SqwaveOdba:maxF_SqwaveOdba);
            fSqwaveOdav_    = fThreshSqwaveBestImfBpOdav(1,minF_SqwaveOdav:maxF_SqwaveOdav);
            
            % recast the WSSTs which are 1:416 x 1:lots

            wsstBpOdba_     = wsstBpOdba(minF_BpOdba:maxF_BpOdba,:);
            wsstBpOdav_     = wsstBpOdav(minF_BpOdav:maxF_BpOdav,:);
            wsstHrBpOdba_   = wsstBpOdba(minF_HrBpOdba:maxF_HrBpOdba,:);
            wsstHrBpOdav_   = wsstBpOdba(minF_HrBpOdav:maxF_HrBpOdav,:);
            wsstSqwaveOdba_ = wsstBpOdba(minF_SqwaveOdba:maxF_SqwaveOdba,:);
            wsstSqwaveOdav_ = wsstBpOdba(minF_SqwaveOdav:maxF_SqwaveOdav,:);
            
            useReducedRidgeSpace = true;
            
        else
            fprintf('There is a mismatch in size of various Fs. Figure this out.\n');
            return;
        end
        
    otherwise
        
        fprintf('Cool, proceeding with entire frequency space...\n');
        useReducedRidgeSpace = false;
        
end

% compute wsstridges for all signals of interest...

if (useReducedRidgeSpace)
    
    fprintf('Starting *reduced* WSST ridge computations for ridge 1 and 2...\n');
    fprintf('\t*** BE PATIENT! This might takes 30-60 seconds!\n');

    tic

    fprintf('Computing WSST ridges for bpOdba...\n');
    [ridgeBpOdba, iRidgeBpOdba] = wsstridge(wsstBpOdba_, 5, fBpOdba_, ...
        'NumRidges', 1);

    fprintf('Computing WSST ridges for hrBpOdba (IMF %d of bpOdba)...\n', ...
        bestImfOdba);
    [ridgeBestImfBpOdba, iRidgeBestImfBpOdba] = wsstridge(wsstHrBpOdba_, 5, ...
        fHrBpOdba_, 'NumRidges', 2);

    fprintf('Computing WSST ridges for sqwave of hrBpOdba...\n');
    [ridgeThreshSqwaveBestImfBpOdba, iRidgeThreshSqwaveBestImfBpOdba] = wsstridge(wsstSqwaveOdba_, ...
        20, fSqwaveOdba_, 'NumRidges', 2);

    fprintf('Computing WSST ridges for bpOdav...\n');
    [ridgeBpOdav, iRidgeBpOdav] = wsstridge(wsstBpOdav_, 30, fBpOdav_, ...
        'NumRidges', 1);

    fprintf('Computing WSST ridges for hrBpOdav (IMF %d of bpOdav)...\n', ...
        bestImfOdav);
    [ridgeBestImfBpOdav, iRidgeBestImfBpOdav] = wsstridge(wsstHrBpOdav_, 60, ...
        fHrBpOdav_, 'NumRidges', 2);

    fprintf('Computing WSST ridges for square wave of hrBpOdav...\n');
    [ridgeThreshSqwaveBestImfBpOdav, iRidgeThreshSqwaveBestImfBpOdav] = wsstridge(wsstSqwaveOdav_, ...
        30, fSqwaveOdav_, 'NumRidges', 2);

    fprintf('Completed WSST ridge computations for all signals of interest...\n');
    toc;
    
else

    fprintf('Starting WSST ridge computations for ridge 1 and 2...\n');
    fprintf('\t*** BE VERY PATIENT! This always takes 60+ seconds!\n');

    tic

    % first, do the bpOdba and bpOdav pair
    
    fprintf('Computing WSST ridges for bpOdba and bpOdav...\n');
    [ridgeBpOdba, iRidgeBpOdba] = wsstridge(wsstBpOdba, 30, fBpOdba, ...
        'NumRidges', 2);    
    [ridgeBpOdav, iRidgeBpOdav] = wsstridge(wsstBpOdav, 30, fBpOdav, ...
        'NumRidges', 2);    

    % second, do the sumimf bpOdba and bpOdav pair
    
    fprintf('Computing WSST ridges for sumImf of bpOdba & bpOdav)...\n');
    [ridgeSumImfBpOdba, iRidgeSumImfBpOdba] = wsstridge(wsstSumImfBpOdba, ...
        30, fSumImfBpOdba, 'NumRidges', 2);       
    [ridgeSumImfBpOdav, iRidgeSumImfBpOdav] = wsstridge(wsstSumImfBpOdav, ...
        30, fSumImfBpOdav, 'NumRidges', 2);  
    
    % third, do the bestImf bpOdba and bpOdav pair
    
    fprintf('Computing WSST ridges for bestImfs of bpOdba (%d) & bpOdav (%d)...\n', ...
        bestImfOdba, bestImfOdav);
    [ridgeBestImfBpOdba, iRidgeBestImfBpOdba] = ...
        wsstridge(wsstBestImfBpOdba, 60, fBestImfBpOdba, 'NumRidges', 2);        
    [ridgeBestImfBpOdav, iRidgeBestImfBpOdav] = ...
        wsstridge(wsstBestImfBpOdav, 60, fBestImfBpOdav, 'NumRidges', 2);
    
    % fourth, do threshold-based square-wave of bestImfs of bpOdba & bpOdav
    
    fprintf('Computing WSST ridges for threshold-based bestImfBpOd** square waves...\n');
    [ridgeThreshSqwaveBestImfBpOdba, iRidgeThreshSqwaveBestImfBpOdba] = ...
        wsstridge(wsstThreshSqwaveBestImfBpOdba, 30, ...
        fThreshSqwaveBestImfBpOdba, 'NumRidges', 2);
    [ridgeThreshSqwaveBestImfBpOdav, iRidgeThreshSqwaveBestImfBpOdav] = ...
        wsstridge(wsstThreshSqwaveBestImfBpOdav, 30, ...
        fThreshSqwaveBestImfBpOdav, 'NumRidges', 2);    
    
    % fifth, and finally, do peak-based square-wave bestImf bpOdba & bpOdav
    
    fprintf('Computing WSST ridges for peak-based bestImfBpOd** square waves...\n');
    [ridgePeakSqwaveBestImfBpOdba, iRidgePeakSqwaveBestImfBpOdba] = ...
        wsstridge(wsstPeakSqwaveBestImfBpOdba, 5, ...
        fPeakSqwaveBestImfBpOdba, 'NumRidges', 2);
    [ridgePeakSqwaveBestImfBpOdav, iRidgePeakSqwaveBestImfBpOdav] = ...
        wsstridge(wsstPeakSqwaveBestImfBpOdav, 5, ...
        fPeakSqwaveBestImfBpOdav, 'NumRidges', 2);
    
    fprintf('Completed WSST ridge computations for all signals of interest...\n');
    
    toc;    
    
end


%% examine first and second ridges for these respective predictors

figRidges = figure('Position', [ 0 0 screenSize(3)/2 screenSize(4)-100 ]);

% wsstridges for bpOdba

pRidge1 = subplot(521);
plot(kT, ridgeBpOdba(:,1), 'Color', blue);
hold on;
plot(kT, ridgeBpOdba(:,2), 'Color', red);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
ylim([0 4]);
title('WSST ridges for bpOdba');
grid;
legend('r1', 'r2');

% wsstridges for sumImfBpOdba

pRidge3 = subplot(523);
plot(kT, ridgeSumImfBpOdba(:,1), 'Color', blue);
hold on;
plot(kT, ridgeSumImfBpOdba(:,2), 'Color', red);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
ylim([0 4]);
title('WSST ridges for SumImfBpOdba');
grid;
legend('r1', 'r2');

% wsstridges for bestImf bpOdba

pRidge5 = subplot(525);
plot(kT, ridgeBestImfBpOdba(:,1), 'Color', blue);
hold on;
plot(kT, ridgeBestImfBpOdba(:,2), 'Color', red);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
ylim([0 4]);
title('WSST ridges for BestImfBpOdba');
grid;
legend('r1', 'r2');

% wsstridges for threshold-based square wave of bestImfBpOdba

pRidge7 = subplot(527);
plot(kT, ridgeThreshSqwaveBestImfBpOdba(:,1), 'Color', blue);
hold on;
plot(kT, ridgeThreshSqwaveBestImfBpOdba(:,2), 'Color', red);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
ylim([0 4]);
title('WSST ridges for threshold-based square waves of bestImfBpOdba');
grid;
legend('r1', 'r2');

% wsstridges for peak-based square wave of bestImfBpOdba
pRidge9 = subplot(529);
plot(kT, ridgePeakSqwaveBestImfBpOdba(:,1), 'Color', blue);
hold on;
plot(kT, ridgePeakSqwaveBestImfBpOdba(:,2), 'Color', red);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
ylim([0 4]);
title('WSST ridges for peak-based square waves of bestImfBpOdba');
grid;
legend('r1', 'r2');

% now the bpOdav-based ones...

% wsstridges for bpOdav

pRidge2 = subplot(522);
plot(kT, ridgeBpOdav(:,1), 'Color', green);
hold on;
plot(kT, ridgeBpOdav(:,2), 'Color', maroon);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
ylim([0 4]);
title('WSST ridges for bpOdav');
grid;
legend('r1', 'r2');

% wsstridges for sumImfBpOdav

pRidge4 = subplot(524);
plot(kT, ridgeSumImfBpOdav(:,1), 'Color', green);
hold on;
plot(kT, ridgeSumImfBpOdav(:,2), 'Color', maroon);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
ylim([0 4]);
title('WSST ridges for SumImfBpOdav');
grid;
legend('r1', 'r2');

% wsstridges for bestImfBpOdav

pRidge6 = subplot(526);
plot(kT, ridgeBestImfBpOdav(:,1), 'Color', green);
hold on;
plot(kT, ridgeBestImfBpOdav(:,2), 'Color', maroon);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
ylim([0 4]);
title('WSST ridges for BestImfBpOdav');
grid;
legend('r1', 'r2');

% wsstridges for threshold-based square wave of bestImfBpOdav

pRidge8 = subplot(528);
plot(kT, ridgeThreshSqwaveBestImfBpOdav(:,1), 'Color', green);
hold on;
plot(kT, ridgeThreshSqwaveBestImfBpOdav(:,2), 'Color', maroon);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
ylim([0 4]);
title('WSST ridges for threshold-based square waves of bestImfBpOdav');
grid;
legend('r1', 'r2');

% wsstridges for peak-based square wave of bestImfBpOdav
pRidge9 = subplot(5,2,10);
plot(kT, ridgePeakSqwaveBestImfBpOdav(:,1), 'Color', green);
hold on;
plot(kT, ridgePeakSqwaveBestImfBpOdav(:,2), 'Color', maroon);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
ylim([0 4]);
title('WSST ridges for peak-based square waves of bestImfBpOdba');
grid;
legend('r1', 'r2');


%% Have the user confirm which ridge corresponds with HR ridges

%   This step will ideally be automated eventually, but for now it's tricky
%   for the ridge detection associated with heart rate to be automatically
%   detected in a correct way. So prompt user to specify which ridge is the
%   HR ridge, if any, to use for consensus HR estimation


confirmed = false;

while(~confirmed)

bestHrRidgeBpOdba                       = 100;
bestHrRidgeBpOdav                       = 100;
bestHrRidgeSumImfBpOdba                 = 100;
bestHrRidgeSumImfBpOdav                 = 100;
bestHrRidgeBestImfBpOdba                = 100;
bestHrRidgeBestImfBpOdav                = 100;
bestHrRidgeThreshSqwaveBestImfBpOdba    = 100;
bestHrRidgeThreshSqwaveBestImfBpOdav    = 100;
bestHrRidgePeakSqwaveBestImfBpOdba      = 100;
bestHrRidgePeakSqwaveBestImfBpOdav      = 100;


fprintf('For each of the following, choose the ridge that best fits HR...\n');
fprintf('Valid ridge values are:\n');
fprintf('\t0 = neither ridge\n');
fprintf('\t1 = ridge 1\n');
fprintf('\t2 = ridge 2\n');

% bpOdba

while (bestHrRidgeBpOdba == 100)
    bestHrRidgeBpOdba = str2double(input('Best HR ridge for bpOdba: ','s'));
    switch(bestHrRidgeBpOdba)
        case 0
            fprintf('\tYou selected neither ridge for bpOdba.\n');
        case 1
            fprintf('\tYou selected ridge 1 for bpOdba.\n')
        case 2
            fprintf('\tYou selected ridge 2 for bpOdba.\n')
        otherwise
            fprintf('\tYou selected an invalid option. Try again.\n');
            bestHrRidgeBpOdba = 100;
    end
end


% bpOdav

while (bestHrRidgeBpOdav == 100)
    bestHrRidgeBpOdav = str2double(input('Best HR ridge for bpOdav: ','s'));
    switch(bestHrRidgeBpOdav)
        case 0
            fprintf('\tYou selected neither ridge for bpOdav.\n');
        case 1
            fprintf('\tYou selected ridge 1 for bpOdav.\n')
        case 2
            fprintf('\tYou selected ridge 2 for bpOdav.\n')
        otherwise
            fprintf('\tYou selected an invalid option. Try again.\n');
            bestHrRidgeBpOdav = 100;
    end
end

% sumImfBpOdba

while (bestHrRidgeSumImfBpOdba == 100)
    bestHrRidgeSumImfBpOdba = str2double(input('Best HR ridge for sumImfBpOdba: ','s'));
    switch(bestHrRidgeSumImfBpOdba)
        case 0
            fprintf('\tYou selected neither ridge for sumImfBpOdba.\n');
        case 1
            fprintf('\tYou selected ridge 1 for sumImfBpOdba.\n')
        case 2
            fprintf('\tYou selected ridge 2 for sumImfBpOdba.\n')
        otherwise
            fprintf('\tYou selected an invalid option. Try again.\n');
            bestHrRidgeSumImfBpOdba = 100;
    end
end


% sumImfBpOdav

while (bestHrRidgeSumImfBpOdav == 100)
    bestHrRidgeSumImfBpOdav = str2double(input('Best HR ridge for sumImfBpOdav: ','s'));
    switch(bestHrRidgeSumImfBpOdav)
        case 0
            fprintf('\tYou selected neither ridge for sumImfBpOdav.\n');
        case 1
            fprintf('\tYou selected ridge 1 for sumImfBpOdav.\n')
        case 2
            fprintf('\tYou selected ridge 2 for sumImfBpOdav.\n')
        otherwise
            fprintf('\tYou selected an invalid option. Try again.\n');
            bestHrRidgeSumImfBpOdav = 100;
    end
end

% bestImfBpOdba

while (bestHrRidgeBestImfBpOdba == 100)
    bestHrRidgeBestImfBpOdba = str2double(input('Best HR ridge for bestImfBpOdba: ','s'));
    switch(bestHrRidgeBestImfBpOdba)
        case 0
            fprintf('\tYou selected neither ridge for bestImfBpOdba.\n');
        case 1
            fprintf('\tYou selected ridge 1 for bestImfBpOdba.\n')
        case 2
            fprintf('\tYou selected ridge 2 for bestImfBpOdba.\n')
        otherwise
            fprintf('\tYou selected an invalid option. Try again.\n');
            bestHrRidgeBestImfBpOdba = 100;
    end
end


% bestImfBpOdav

while (bestHrRidgeBestImfBpOdav == 100)
    bestHrRidgeBestImfBpOdav = str2double(input('Best HR ridge for bestImfBpOdav: ','s'));
    switch(bestHrRidgeBestImfBpOdav)
        case 0
            fprintf('\tYou selected neither ridge for bestImfBpOdav.\n');
        case 1
            fprintf('\tYou selected ridge 1 for bestImfBpOdav.\n')
        case 2
            fprintf('\tYou selected ridge 2 for bestImfBpOdav.\n')
        otherwise
            fprintf('\tYou selected an invalid option. Try again.\n');
            bestHrRidgeBestImfBpOdav = 100;
    end
end

% threshSqwaveBestImfBpOdba

while (bestHrRidgeThreshSqwaveBestImfBpOdba == 100)
    bestHrRidgeThreshSqwaveBestImfBpOdba = ...
        str2double(input('Best HR ridge for threshSqwaveBestImfBpOdba: ','s'));
    switch(bestHrRidgeThreshSqwaveBestImfBpOdba)
        case 0
            fprintf('\tYou selected neither ridge for threshSqwaveBestImfBpOdba.\n');
        case 1
            fprintf('\tYou selected ridge 1 for threshSqwaveBestImfBpOdba.\n')
        case 2
            fprintf('\tYou selected ridge 2 for threshSqwaveBestImfBpOdba.\n')
        otherwise
            fprintf('\tYou selected an invalid option. Try again.\n');
            bestHrRidgeThreshSqwaveBestImfBpOdba = 100;
    end
end


% threshSqwaveBestImfBpOdav

while (bestHrRidgeThreshSqwaveBestImfBpOdav == 100)
    bestHrRidgeThreshSqwaveBestImfBpOdav = ...
        str2double(input('Best HR ridge for threshSqwaveBestImfBpOdav: ','s'));
    switch(bestHrRidgeThreshSqwaveBestImfBpOdav)
        case 0
            fprintf('\tYou selected neither ridge for threshSqwaveBestImfBpOdav.\n');
        case 1
            fprintf('\tYou selected ridge 1 for threshSqwaveBestImfBpOdav.\n')
        case 2
            fprintf('\tYou selected ridge 2 for threshSqwaveBestImfBpOdav.\n')
        otherwise
            fprintf('\tYou selected an invalid option. Try again.\n');
            bestHrRidgeThreshSqwaveBestImfBpOdav = 100;
    end
end

% peakSqwaveBestImfBpOdba

while (bestHrRidgePeakSqwaveBestImfBpOdba == 100)
    bestHrRidgePeakSqwaveBestImfBpOdba = ...
        str2double(input('Best HR ridge for peakSqwaveBestImfBpOdba: ','s'));
    switch(bestHrRidgePeakSqwaveBestImfBpOdba)
        case 0
            fprintf('\tYou selected neither ridge for peakSqwaveBestImfBpOdba.\n');
        case 1
            fprintf('\tYou selected ridge 1 for peakSqwaveBestImfBpOdba.\n')
        case 2
            fprintf('\tYou selected ridge 2 for peakSqwaveBestImfBpOdba.\n')
        otherwise
            fprintf('\tYou selected an invalid option. Try again.\n');
            bestHrRidgePeakSqwaveBestImfBpOdba = 100;
    end
end

% peakSqwaveBestImfBpOdav

while (bestHrRidgePeakSqwaveBestImfBpOdav == 100)
    bestHrRidgePeakSqwaveBestImfBpOdav = ...
        str2double(input('Best HR ridge for peakSqwaveBestImfBpOdav: ','s'));
    switch(bestHrRidgePeakSqwaveBestImfBpOdav)
        case 0
            fprintf('\tYou selected neither ridge for peakSqwaveBestImfBpOdav.\n');
        case 1
            fprintf('\tYou selected ridge 1 for peakSqwaveBestImfBpOdav.\n')
        case 2
            fprintf('\tYou selected ridge 2 for peakSqwaveBestImfBpOdav.\n')
        otherwise
            fprintf('\tYou selected an invalid option. Try again.\n');
            bestHrRidgePeakSqwaveBestImfBpOdav = 100;
    end
end

fprintf('Best HR ridge summary:\n');
fprintf('\tbpOdba: ridge %d\n', bestHrRidgeBpOdba);
fprintf('\tbpOdav: ridge %d\n', bestHrRidgeBpOdav);
fprintf('\tsumImfBpOdba: ridge %d\n', bestHrRidgeSumImfBpOdba);
fprintf('\tsumImfBpOdav: ridge %d\n', bestHrRidgeSumImfBpOdav);
fprintf('\tbestImfBpOdba: ridge %d\n', bestHrRidgeBestImfBpOdba);
fprintf('\tbestImfBpOdav: ridge %d\n', bestHrRidgeBestImfBpOdav);
fprintf('\tthreshSqwaveBestImfBpOdba: ridge %d\n', ...
    bestHrRidgeThreshSqwaveBestImfBpOdba);
fprintf('\tthreshSqwaveBestImfBpOdav: ridge %d\n', ...
    bestHrRidgeThreshSqwaveBestImfBpOdav);
fprintf('\tpeakSqwaveBestImfBpOdba: ridge %d\n', ...
    bestHrRidgePeakSqwaveBestImfBpOdba);
fprintf('\tpeakSqwaveBestImfBpOdav: ridge %d\n', ...
    bestHrRidgePeakSqwaveBestImfBpOdav);


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





%% Bland-Altman plots of key measures

figSizeBA = [ 100 100 600 600 ];

if (bestHrRidgeBpOdba ~= 0 & bestHrRidgeBpOdav ~= 0)

fig1 = figure('Position', figSizeBA, 'Color', white, 'NumberTitle', 'off');
[ rpc1, fig1, stats1 ] = BlandAltman(fig1, ...
    ridgeBpOdba(:,bestHrRidgeBpOdba),ridgeBpOdav(:,bestHrRidgeBpOdav), ...
    {'bpOdba','bpOdav'}, ...
    'Bland-Altman test for bpOdba:bpOdav - WSST 1st Ridge Frequency', ...
    {'bpOdba:bpOdav','LoBF'}, ...
    'corrInfo', {'eq','r','r2','p','SSE','RMSE','n'}, ...
    'baStatsMode','Non-parametric');

else
    fprintf('Skipping Bland-Altman for bpOd** pair.\n');
end


if (bestHrRidgeSumImfBpOdba ~= 0 & bestHrRidgeSumImfBpOdav ~= 0)

fig3 = figure('Position', figSizeBA, 'Color', white, 'NumberTitle', 'off');
[ rpc3, fig3, stats3 ]  = BlandAltman(fig3, ...
    ridgeSumImfBpOdba(:,bestHrRidgeSumImfBpOdba),...
    ridgeSumImfBpOdav(:,bestHrRidgeSumImfBpOdav), ...
    {'sumImfBpOdba','sumImfBpOdav'}, ...
    'Bland-Altman test for sumImfBpOd** - WSST 1st Ridge', ...
    {'sumImfBpOd**','LoBF'}, ...
    'corrInfo', {'eq','r','r2','p','SSE','RMSE','n'}, ...
    'baStatsMode','Non-parametric');

else
    fprintf('Skipping Bland-Altman for sumImfBpOd** pair.\n');
end


if (bestHrRidgeBestImfBpOdba ~= 0 & bestHrRidgeBestImfBpOdav ~= 0)
    
fig5 = figure('Position', figSizeBA, 'Color', white, 'NumberTitle', 'off');
[ rpc5, fig5, stats5 ] = BlandAltman(fig5, ...
    ridgeBestImfBpOdba(:,bestHrRidgeBestImfBpOdba),...
    ridgeBestImfBpOdav(:,bestHrRidgeBestImfBpOdav), ...
    {'bestImfBpOdba','bestImfBpOdav'}, ...
    'Bland-Altman test for BestImfBpOd** - WSST 1st Ridge', ...
    {'bestImfBpOd**','LoBF'}, ...
    'corrInfo', {'eq','r','r2','p','SSE','RMSE','n'}, ...
    'baStatsMode','Non-parametric');

else
    fprintf('Skipping Bland-Altman for bestImfBpOd** pair.\n');
end

if (bestHrRidgeThreshSqwaveBestImfBpOdba ~= 0 & ...
        bestHrRidgeThreshSqwaveBestImfBpOdav ~= 0)
    
fig7 = figure('Position', figSizeBA, 'Color', white, 'NumberTitle', 'off');
[ rpc7, fig7, stats7 ]  = BlandAltman(fig7, ...
    ridgeThreshSqwaveBestImfBpOdba(:,bestHrRidgeThreshSqwaveBestImfBpOdba), ...
    ridgeThreshSqwaveBestImfBpOdav(:,bestHrRidgeThreshSqwaveBestImfBpOdav), ...
    {'threshSqwaveOdba','threshSqwaveOdav'}, ...
    'Bland-Altman test for threshSqwaveBestImfBpOd** - WSST 1st Ridge', ...
    {'threshSqwaveBestImfBpOd**','LoBF'}, ...
    'corrInfo', {'eq','r','r2','p','SSE','RMSE','n'}, ...
    'baStatsMode','Non-parametric');

else
    fprintf('Skipping Bland-Altman for threshSqwaveBestImfBpOd** pair.\n');
end

if (bestHrRidgePeakSqwaveBestImfBpOdba ~= 0 & ...
        bestHrRidgePeakSqwaveBestImfBpOdav ~= 0)

fig9 = figure('Position', figSizeBA, 'Color', white, 'NumberTitle', 'off');
[ rpc9, fig9, stats9 ]  = BlandAltman(fig9, ...
    ridgePeakSqwaveBestImfBpOdba(:,bestHrRidgePeakSqwaveBestImfBpOdba), ...
    ridgePeakSqwaveBestImfBpOdav(:,bestHrRidgePeakSqwaveBestImfBpOdav), ...
    {'peakSqwaveOdba','peakSqwaveOdav'}, ...
    'Bland-Altman test for peakSqwaveBestImfBpOd** - WSST 1st Ridge', ...
    {'peakSqwaveBestImfBpOd**','LoBF'}, ...
    'corrInfo', {'eq','r','r2','p','SSE','RMSE','n'}, ...
    'baStatsMode','Non-parametric');

else
    fprintf('Skipping Bland-Altman for peakSqwaveBestImfBpOd** pair.\n');
end

%% Calculate mean and standard deviations for bpOdba-based ones first...

meanBpOdbaRidges = [ mean(ridgeBpOdba(:,1)), mean(ridgeBpOdba(:,2)) ];
stdBpOdbaRidges = [ std(ridgeBpOdba(:,1)), std(ridgeBpOdba(:,2)) ];

meanSumImfBpOdbaRidges = [ mean(ridgeSumImfBpOdba(:,1)), ...
    mean(ridgeSumImfBpOdba(:,2)) ];
stdSumImfBpOdbaRidges = [ std(ridgeSumImfBpOdba(:,1)), ...
    std(ridgeSumImfBpOdba(:,2)) ];

meanBestImfBpOdbaRidges = [ mean(ridgeBestImfBpOdba(:,1)), ...
    mean(ridgeBestImfBpOdba(:,2)) ];
stdBestImfBpOdbaRidges = [ std(ridgeBestImfBpOdba(:,1)), ...
    std(ridgeBestImfBpOdba(:,2)) ];

meanThreshSqwaveBestImfBpOdbaRidges = [ mean(ridgeThreshSqwaveBestImfBpOdba(:,1)), ...
    mean(ridgeThreshSqwaveBestImfBpOdba(:,2)) ];
stdThreshSqwaveBestImfBpOdbaRidges = [ std(ridgeThreshSqwaveBestImfBpOdba(:,1)), ...
    std(ridgeThreshSqwaveBestImfBpOdba(:,2)) ];

meanPeakSqwaveBestImfBpOdbaRidges = ...
    [ mean(ridgePeakSqwaveBestImfBpOdba(:,1)), ...
    mean(ridgePeakSqwaveBestImfBpOdba(:,2)) ];
stdPeakSqwaveBestImfBpOdbaRidges = ...
    [ std(ridgePeakSqwaveBestImfBpOdba(:,1)), ...
    std(ridgePeakSqwaveBestImfBpOdba(:,2)) ];

% do the bpOdav-based ones next...

meanBpOdavRidges = [ mean(ridgeBpOdav(:,1)), mean(ridgeBpOdav(:,2)) ];
stdBpOdavRidges = [ std(ridgeBpOdav(:,1)), std(ridgeBpOdav(:,2)) ];

meanSumImfBpOdavRidges = [ mean(ridgeSumImfBpOdav(:,1)), ...
    mean(ridgeSumImfBpOdav(:,2)) ];
stdSumImfBpOdavRidges = [ std(ridgeSumImfBpOdav(:,1)), ...
    std(ridgeSumImfBpOdav(:,2)) ];

meanBestImfBpOdavRidges = [ mean(ridgeBestImfBpOdav(:,1)), ...
    mean(ridgeBestImfBpOdav(:,2)) ];
stdBestImfBpOdavRidges = [ std(ridgeBestImfBpOdav(:,1)), ...
    std(ridgeBestImfBpOdav(:,2)) ];

meanThreshSqwaveBestImfBpOdavRidges = ...
    [ mean(ridgeThreshSqwaveBestImfBpOdav(:,1)), ...
    mean(ridgeThreshSqwaveBestImfBpOdav(:,2)) ];
stdThreshSqwaveBestImfBpOdavRidges = ...
    [ std(ridgeThreshSqwaveBestImfBpOdav(:,1)), ...
    std(ridgeThreshSqwaveBestImfBpOdav(:,2)) ];

meanPeakSqwaveBestImfBpOdavRidges = ...
    [ mean(ridgePeakSqwaveBestImfBpOdav(:,1)), ...
    mean(ridgePeakSqwaveBestImfBpOdav(:,2)) ];
stdPeakSqwaveBestImfBpOdavRidges = ...
    [ std(ridgePeakSqwaveBestImfBpOdav(:,1)), ...
    std(ridgePeakSqwaveBestImfBpOdav(:,2)) ];

% fprintf('Some stats on the frequencies of the ridges...\n');
% fprintf('\tbpOdba mean [sd]: %2.2f Hz [%2.2f Hz]\n', ...
%     meanBpOdbaRidges(bestHrRidgeBpOdba), ...
%     stdBpOdbaRidges(bestHrRidgeBpOdba) );

%% Make histograms

figHrFreq = figure('Position',[0 0 screenSize(3)/3 screenSize(4) ], ...
    'Color', white, 'NumberTitle', 'off');

sHrF1 = subplot(521);
if (bestHrRidgeBpOdba ~=0)
    histogram( ridgeBpOdba(:,bestHrRidgeBpOdba), 'Normalization', 'pdf', ...
        'FaceColor', red, 'EdgeColor', black, 'Orientation', 'horizontal');
    hold on;
    yline(meanBpOdbaRidges(bestHrRidgeBpOdba), 'Color', blue, ...
        'LineWidth', 1.2);
    yline(meanBpOdbaRidges(bestHrRidgeBpOdba) + ...
        1.96 * stdBpOdbaRidges(bestHrRidgeBpOdba), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    yline(meanBpOdbaRidges(bestHrRidgeBpOdba) + ...
        -1.96 * stdBpOdbaRidges(bestHrRidgeBpOdba), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    meanStr = sprintf('%1.2f Hz', meanBpOdbaRidges(bestHrRidgeBpOdba));
    text(3.5,1.35, ['{\mu}: ' meanStr],'FontSize',12);
    stdStr = sprintf('%1.2f Hz', stdBpOdbaRidges(bestHrRidgeBpOdba));
    text(3.5,1.25, ['{\sigma}: ' stdStr],'FontSize',12);
    % text(4, 1.2, meanStr, 'FontSize', 14);
    hold off;
    ylabel('Frequency, Hz');
    grid;
end

sHrF2 = subplot(522);
if (bestHrRidgeBpOdav ~=0)
    histogram(ridgeBpOdav(:,bestHrRidgeBpOdav), 'Normalization', 'pdf', ...
        'FaceColor', green, 'EdgeColor', black, 'Orientation', 'horizontal');
    hold on;
    yline(meanBpOdavRidges(bestHrRidgeBpOdav), 'Color', blue, ...
        'LineWidth', 1.2);
    yline(meanBpOdavRidges(bestHrRidgeBpOdav) + ...
        1.96 * stdBpOdavRidges(bestHrRidgeBpOdav), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    yline(meanBpOdavRidges(bestHrRidgeBpOdav) + ...
        -1.96 * stdBpOdavRidges(bestHrRidgeBpOdav), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    meanStr = sprintf('%1.2f Hz', meanBpOdavRidges(bestHrRidgeBpOdav));
    text(3.5,1.35, ['{\mu}: ' meanStr],'FontSize',12);
    stdStr = sprintf('%1.2f Hz', stdBpOdavRidges(bestHrRidgeBpOdav));
    text(3.5,1.25, ['{\sigma}: ' stdStr],'FontSize',12);
    hold off;
    ylabel('Frequency, Hz');
    grid;
end

sHrF3 = subplot(523);
if (bestHrRidgeSumImfBpOdba ~= 0)
    histogram( ridgeSumImfBpOdba(:,bestHrRidgeSumImfBpOdba), ...
        'Normalization', 'pdf', 'FaceColor', red, 'EdgeColor', black, ...
        'Orientation', 'horizontal');
    hold on;
    yline(meanSumImfBpOdbaRidges(bestHrRidgeSumImfBpOdba), 'Color', blue, ...
        'LineWidth', 1.2);
    yline(meanSumImfBpOdbaRidges(bestHrRidgeSumImfBpOdba) + ...
        1.96 * stdSumImfBpOdbaRidges(bestHrRidgeSumImfBpOdba), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    yline(meanSumImfBpOdbaRidges(bestHrRidgeSumImfBpOdba) + ...
        -1.96 * stdSumImfBpOdbaRidges(bestHrRidgeSumImfBpOdba), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    meanStr = sprintf('%1.2f Hz', ...
        meanSumImfBpOdbaRidges(bestHrRidgeSumImfBpOdba));
    text(3.5,1.35, ['{\mu}: ' meanStr],'FontSize',12);
    stdStr = sprintf('%1.2f Hz', ...
        stdSumImfBpOdbaRidges(bestHrRidgeSumImfBpOdba));
    text(3.5,1.25, ['{\sigma}: ' stdStr],'FontSize',12);
    hold off;
    ylabel('Frequency, Hz');
    grid;
end

sHrF4 = subplot(524);
if (bestHrRidgeSumImfBpOdav ~= 0)
    histogram( ridgeSumImfBpOdav(:,bestHrRidgeSumImfBpOdav), ...
        'Normalization', 'pdf', 'FaceColor', green, 'EdgeColor', black, ...
        'Orientation', 'horizontal');
    hold on;
    yline(meanSumImfBpOdavRidges(bestHrRidgeSumImfBpOdav), 'Color', blue,...
        'LineWidth', 1.2);
    yline(meanSumImfBpOdavRidges(bestHrRidgeSumImfBpOdav) + ...
        1.96 * stdSumImfBpOdavRidges(bestHrRidgeSumImfBpOdav), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    yline(meanSumImfBpOdavRidges(bestHrRidgeSumImfBpOdav) + ...
        -1.96 * stdSumImfBpOdavRidges(bestHrRidgeSumImfBpOdav), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    meanStr = sprintf('%1.2f Hz', ...
        meanSumImfBpOdavRidges(bestHrRidgeSumImfBpOdav));
    text(3.5,1.35, ['{\mu}: ' meanStr],'FontSize',12);
    stdStr = sprintf('%1.2f Hz', ...
        stdSumImfBpOdavRidges(bestHrRidgeSumImfBpOdav));
    text(3.5,1.25, ['{\sigma}: ' stdStr],'FontSize',12);
    hold off;
    ylabel('Frequency, Hz');
    grid;
end

sHrF5 = subplot(525);
if (bestHrRidgeBestImfBpOdba ~= 0)
    histogram( ridgeBestImfBpOdba(:,bestHrRidgeBestImfBpOdba), ...
        'Normalization', 'pdf', 'FaceColor', red, 'EdgeColor', black, ...
        'Orientation', 'horizontal');
    hold on;
    yline(meanBestImfBpOdbaRidges(bestHrRidgeBestImfBpOdba), 'Color', blue, 'LineWidth', 1.2);
    yline(meanBestImfBpOdbaRidges(bestHrRidgeBestImfBpOdba) + ...
        1.96 * stdBestImfBpOdbaRidges(bestHrRidgeBestImfBpOdba), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    yline(meanBestImfBpOdbaRidges(bestHrRidgeBestImfBpOdba) + ...
        -1.96 * stdBestImfBpOdbaRidges(bestHrRidgeBestImfBpOdba), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    meanStr = sprintf('%1.2f Hz', meanBestImfBpOdbaRidges(bestHrRidgeBestImfBpOdba));
    text(3.5,1.35, ['{\mu}: ' meanStr],'FontSize',12);
    stdStr = sprintf('%1.2f Hz', stdBestImfBpOdbaRidges(bestHrRidgeBestImfBpOdba));
    text(3.5,1.25, ['{\sigma}: ' stdStr],'FontSize',12);
    hold off;
    ylabel('Frequency, Hz');
    grid;
end

sHrF6 = subplot(526);
if (bestHrRidgeBestImfBpOdav ~= 0)
    histogram( ridgeBestImfBpOdav(:,bestHrRidgeBestImfBpOdav), ...
        'Normalization', 'pdf', 'FaceColor', green, 'EdgeColor', black, ...
        'Orientation', 'horizontal');
    hold on;
    yline(meanBestImfBpOdavRidges(bestHrRidgeBestImfBpOdav), ...
        'Color', blue, 'LineWidth', 1.2);
    yline(meanBestImfBpOdavRidges(bestHrRidgeBestImfBpOdav) + ...
        1.96 * stdBestImfBpOdavRidges(bestHrRidgeBestImfBpOdav), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    yline(meanBestImfBpOdavRidges(bestHrRidgeBestImfBpOdav) + ...
        -1.96 * stdBestImfBpOdavRidges(bestHrRidgeBestImfBpOdav), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    meanStr = sprintf('%1.2f Hz', ...
        meanBestImfBpOdavRidges(bestHrRidgeBestImfBpOdav));
    text(3.5,1.35, ['{\mu}: ' meanStr],'FontSize',12);
    stdStr = sprintf('%1.2f Hz', ...
        stdBestImfBpOdavRidges(bestHrRidgeBestImfBpOdav));
    text(3.5,1.25, ['{\sigma}: ' stdStr],'FontSize',12);
    hold off;
    ylabel('Frequency, Hz');
    grid;
end

sHrF7 = subplot(527);

if (bestHrRidgeThreshSqwaveBestImfBpOdba ~= 0)
    histogram(...
        ridgeThreshSqwaveBestImfBpOdba(:,bestHrRidgeThreshSqwaveBestImfBpOdba ), ...
        'Normalization', 'pdf', 'FaceColor', red, 'EdgeColor', black, ...
        'Orientation', 'horizontal');
    hold on;
    yline(...
        meanThreshSqwaveBestImfBpOdbaRidges(bestHrRidgeThreshSqwaveBestImfBpOdba), ...
        'Color', blue, 'LineWidth', 1.2);
    yline(...
        meanThreshSqwaveBestImfBpOdbaRidges(bestHrRidgeThreshSqwaveBestImfBpOdba) + ...
        1.96 * stdThreshSqwaveBestImfBpOdbaRidges(bestHrRidgeThreshSqwaveBestImfBpOdba), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    yline(...
        meanThreshSqwaveBestImfBpOdbaRidges(bestHrRidgeThreshSqwaveBestImfBpOdba) + ...
        -1.96 * ...
        stdThreshSqwaveBestImfBpOdbaRidges(bestHrRidgeThreshSqwaveBestImfBpOdba), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    meanStr = sprintf('%1.2f Hz', ...
        meanThreshSqwaveBestImfBpOdbaRidges(bestHrRidgeThreshSqwaveBestImfBpOdba));
    text(3.5,1.35, ['{\mu}: ' meanStr],'FontSize',12);
    stdStr = sprintf('%1.2f Hz', ...
        stdThreshSqwaveBestImfBpOdbaRidges(bestHrRidgeThreshSqwaveBestImfBpOdba));
    text(3.5,1.25, ['{\sigma}: ' stdStr],'FontSize',12);
    hold off;
    ylabel('Frequency, Hz');
    grid;
end

sHrF8 = subplot(528);
if (bestHrRidgeThreshSqwaveBestImfBpOdav ~= 0)
    histogram(ridgeThreshSqwaveBestImfBpOdav(:,bestHrRidgeThreshSqwaveBestImfBpOdav ), ...
        'Normalization', 'pdf', 'FaceColor', green, 'EdgeColor', black, ...
        'Orientation', 'horizontal');
    hold on;
    yline(meanThreshSqwaveBestImfBpOdavRidges(bestHrRidgeThreshSqwaveBestImfBpOdav), ...
        'Color', blue, 'LineWidth', 1.2);
    yline(meanThreshSqwaveBestImfBpOdavRidges(bestHrRidgeThreshSqwaveBestImfBpOdav) + ...
        1.96 * stdThreshSqwaveBestImfBpOdavRidges(bestHrRidgeThreshSqwaveBestImfBpOdav), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    yline(meanThreshSqwaveBestImfBpOdavRidges(bestHrRidgeThreshSqwaveBestImfBpOdav) + ...
        -1.96 * stdThreshSqwaveBestImfBpOdavRidges(bestHrRidgeThreshSqwaveBestImfBpOdav), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    meanStr = sprintf('%1.2f Hz', ...
        meanThreshSqwaveBestImfBpOdavRidges(bestHrRidgeThreshSqwaveBestImfBpOdav));
    text(3.5,1.35, ['{\mu}: ' meanStr],'FontSize',12);
    stdStr = sprintf('%1.2f Hz', ...
        stdThreshSqwaveBestImfBpOdavRidges(bestHrRidgeThreshSqwaveBestImfBpOdav));
    text(3.5,1.25, ['{\sigma}: ' stdStr],'FontSize',12);
    hold off;
    ylabel('Frequency, Hz');
    grid;
end

sHrF9 = subplot(529);
if (bestHrRidgePeakSqwaveBestImfBpOdba ~= 0)
    histogram(ridgePeakSqwaveBestImfBpOdba(:,bestHrRidgePeakSqwaveBestImfBpOdba ), ...
        'Normalization', 'pdf', 'FaceColor', red, 'EdgeColor', black, ...
        'Orientation', 'horizontal');
    hold on;
    yline(meanPeakSqwaveBestImfBpOdbaRidges(bestHrRidgePeakSqwaveBestImfBpOdba), ...
        'Color', blue, 'LineWidth', 1.2);
    yline(meanPeakSqwaveBestImfBpOdbaRidges(bestHrRidgePeakSqwaveBestImfBpOdba) + ...
        1.96 * stdPeakSqwaveBestImfBpOdbaRidges(bestHrRidgePeakSqwaveBestImfBpOdba), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    yline(meanPeakSqwaveBestImfBpOdbaRidges(bestHrRidgePeakSqwaveBestImfBpOdba) + ...
        -1.96 * stdPeakSqwaveBestImfBpOdbaRidges(bestHrRidgePeakSqwaveBestImfBpOdba), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    meanStr = sprintf('%1.2f Hz', ...
        meanPeakSqwaveBestImfBpOdbaRidges(bestHrRidgePeakSqwaveBestImfBpOdba));
    text(3.5,1.35, ['{\mu}: ' meanStr],'FontSize',12);
    stdStr = sprintf('%1.2f Hz', ...
        stdPeakSqwaveBestImfBpOdbaRidges(bestHrRidgePeakSqwaveBestImfBpOdba));
    text(3.5,1.25, ['{\sigma}: ' stdStr],'FontSize',12);
    hold off;
    ylabel('Frequency, Hz');
    grid;
end

sHrF10 = subplot(5,2,10);
if (bestHrRidgePeakSqwaveBestImfBpOdav ~= 0)
    histogram(ridgePeakSqwaveBestImfBpOdav(:,bestHrRidgePeakSqwaveBestImfBpOdav ), ...
        'Normalization', 'pdf', 'FaceColor', green, 'EdgeColor', black, ...
        'Orientation', 'horizontal');
    hold on;
    yline(meanPeakSqwaveBestImfBpOdavRidges(bestHrRidgePeakSqwaveBestImfBpOdav), ...
        'Color', blue, 'LineWidth', 1.2);
    yline(meanPeakSqwaveBestImfBpOdavRidges(bestHrRidgePeakSqwaveBestImfBpOdav) + ...
        1.96 * stdPeakSqwaveBestImfBpOdavRidges(bestHrRidgePeakSqwaveBestImfBpOdav), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    yline(meanPeakSqwaveBestImfBpOdavRidges(bestHrRidgePeakSqwaveBestImfBpOdav) + ...
        -1.96 * stdPeakSqwaveBestImfBpOdavRidges(bestHrRidgePeakSqwaveBestImfBpOdav), ...
        '--', 'Color', blue, 'LineWidth', 1.1);
    meanStr = sprintf('%1.2f Hz', ...
        meanPeakSqwaveBestImfBpOdavRidges(bestHrRidgePeakSqwaveBestImfBpOdav));
    text(3.5,1.35, ['{\mu}: ' meanStr],'FontSize',12);
    stdStr = sprintf('%1.2f Hz', ...
        stdPeakSqwaveBestImfBpOdavRidges(bestHrRidgePeakSqwaveBestImfBpOdav));
    text(3.5,1.25, ['{\sigma}: ' stdStr],'FontSize',12);
    hold off;
    ylabel('Frequency, Hz');
    grid;
end

linkaxes([sHrF1 sHrF2 sHrF3 sHrF4 sHrF5 sHrF6 sHrF7 sHrF8 sHrF9 sHrF10 ],'x');

%% construct a single set of candidate HR ODBA ridges 

numOdbaRidges = ( bestHrRidgeBpOdba > 0 ) + ...
    ( bestHrRidgeSumImfBpOdba > 0 ) + ...
    ( bestHrRidgeBestImfBpOdba > 0 ) + ...
    ( bestHrRidgeThreshSqwaveBestImfBpOdba > 0 ) + ...
    ( bestHrRidgePeakSqwaveBestImfBpOdba > 0 ) ;

hrOdbaRidges = zeros(length(bpOdba), numOdbaRidges);

% keep track of ridgesAdded, to make sure right numRidges are added...

nextRidge = 0;

% build the hrOdbaRidges 

switch(bestHrRidgeBpOdba)
    case 1
        nextRidge = nextRidge + 1;
        hrOdbaRidges(:,nextRidge) = ridgeBpOdba(:,1);
        fprintf('Adding bpOdba ridge 1 to hrOdbaRidges.\n');
    case 2
        nextRidge = nextRidge + 1;
        hrOdbaRidges(:,nextRidge) = ridgeBpOdba(:,2);
        fprintf('Adding bpOdba ridge 2 to hrOdbaRidges.\n');
    otherwise
        fprintf('No bpOdba ridge added to hrOdbaRidges.\n');
end

switch(bestHrRidgeSumImfBpOdba)
    case 1
        nextRidge = nextRidge + 1;
        hrOdbaRidges(:,nextRidge) = ridgeSumImfBpOdba(:,1);
        fprintf('Adding sumImfBpOdba ridge 1 to hrOdbaRidges.\n');
    case 2
        nextRidge = nextRidge + 1;
        hrOdbaRidges(:,nextRidge) = ridgeSumImfBpOdba(:,2);
        fprintf('Adding sumImfBpOdba ridge 2 to hrOdbaRidges.\n');
    otherwise
        fprintf('No sumImfBpOdba ridge added to hrOdbaRidges.\n');
end

switch(bestHrRidgeBestImfBpOdba)
    case 1
        nextRidge = nextRidge + 1;
        hrOdbaRidges(:,nextRidge) = ridgeBestImfBpOdba(:,1);
        fprintf('Adding bestImfBpOdba ridge 1 to hrOdbaRidges.\n');
    case 2
        nextRidge = nextRidge + 1;
        hrOdbaRidges(:,nextRidge) = ridgeBestImfBpOdba(:,2);
        fprintf('Adding bestImfBpOdba ridge 2 to hrOdbaRidges.\n');
    otherwise
        fprintf('No bestImfBpOdba ridge added to hrOdbaRidges.\n');
end

switch(bestHrRidgeThreshSqwaveBestImfBpOdba)
    case 1
        nextRidge = nextRidge + 1;
        hrOdbaRidges(:,nextRidge) = ridgeThreshSqwaveBestImfBpOdba(:,1);
        fprintf('Adding threshSqwaveBestImfBpOdba ridge 1 to hrOdbaRidges.\n');
    case 2
        nextRidge = nextRidge + 1;
        hrOdbaRidges(:,nextRidge) = ridgeThreshSqwaveBestImfBpOdba(:,2);
        fprintf('Adding threshSqwaveBestImfBpOdba ridge 2to hrOdbaRidges.\n');
    otherwise
        fprintf('No threshSqwaveBestImfBpOdba ridge added to hrOdbaRidges.\n');
end

switch(bestHrRidgePeakSqwaveBestImfBpOdba)
    case 1
        nextRidge = nextRidge + 1;
        hrOdbaRidges(:,nextRidge) = ridgePeakSqwaveBestImfBpOdba(:,1);
        fprintf('Adding peakSqwaveBestImfBpOdba ridge 1 to hrOdbaRidges.\n');
    case 2
        nextRidge = nextRidge + 1;
        hrOdbaRidges(:,nextRidge) = ridgePeakSqwaveBestImfBpOdba(:,2);
        fprintf('Adding peakSqwaveBestImfBpOdba ridge 2 to hrOdbaRidges.\n');
    otherwise
        fprintf('No peakSqwaveBestImfBpOdba ridge added to hrOdbaRidges.\n');
end

%% now do the hrOdavRidges

numOdavRidges = ( bestHrRidgeBpOdav > 0 ) + ...
    ( bestHrRidgeSumImfBpOdav > 0 ) + ...
    ( bestHrRidgeBestImfBpOdav > 0 ) + ...
    ( bestHrRidgeThreshSqwaveBestImfBpOdav > 0 ) + ...
    ( bestHrRidgePeakSqwaveBestImfBpOdav > 0 ) ;

hrOdavRidges = zeros(length(bpOdav), numOdavRidges);

nextRidge = 0;

switch(bestHrRidgeBpOdav)
    case 1
        nextRidge = nextRidge + 1;
        hrOdavRidges(:,nextRidge) = ridgeBpOdav(:,1);
        fprintf('Adding BpOdav ridge 1 to hrOdavRidges.\n');
    case 2
        nextRidge = nextRidge + 1;
        hrOdavRidges(:,nextRidge) = ridgeBpOdav(:,2);
        fprintf('Adding BpOdav ridge 2 to hrOdavRidges.\n');
    otherwise
        fprintf('No BpOdav ridge added to hrOdavRidges.\n');
end

switch(bestHrRidgeSumImfBpOdav)
    case 1
        nextRidge = nextRidge + 1;
        hrOdavRidges(:,nextRidge) = ridgeSumImfBpOdav(:,1);
        fprintf('Adding sumImfBpOdav ridge 1 to hrOdavRidges.\n');
    case 2
        nextRidge = nextRidge + 1;
        hrOdavRidges(:,nextRidge) = ridgeSumImfBpOdav(:,2);
        fprintf('Adding sumImfBpOdav ridge 2 to hrOdavRidges.\n');
    otherwise
        fprintf('No sumImfBpOdav ridge added to hrOdavRidges.\n');
end

switch(bestHrRidgeBestImfBpOdav)
    case 1
        nextRidge = nextRidge + 1;
        hrOdavRidges(:,nextRidge) = ridgeBestImfBpOdav(:,1);
        fprintf('Adding bestImfBpOdav ridge 1 to hrOdavRidges.\n');
    case 2
        nextRidge = nextRidge + 1;
        hrOdavRidges(:,nextRidge) = ridgeBestImfBpOdav(:,2);
        fprintf('Adding bestImfBpOdav ridge 2 to hrOdavRidges.\n');
    otherwise
        fprintf('No bestImfBpOdav ridge added to hrOdavRidges.\n');
end

switch(bestHrRidgeThreshSqwaveBestImfBpOdav)
    case 1
        nextRidge = nextRidge + 1;
        hrOdavRidges(:,nextRidge) = ridgeThreshSqwaveBestImfBpOdav(:,1);
        fprintf('Adding threshSqwaveBestImfBpOdav ridge 1 to hrOdavRidges.\n');
    case 2
        nextRidge = nextRidge + 1;
        hrOdavRidges(:,nextRidge) = ridgeThreshSqwaveBestImfBpOdav(:,2);
        fprintf('Adding threshSqwaveBestImfBpOdav ridge 2to hrOdavRidges.\n');
    otherwise
        fprintf('No threshSqwaveBestImfBpOdav ridge added to hrOdavRidges.\n');
end

switch(bestHrRidgePeakSqwaveBestImfBpOdav)
    case 1
        nextRidge = nextRidge + 1;
        hrOdavRidges(:,nextRidge) = ridgePeakSqwaveBestImfBpOdav(:,1);
        fprintf('Adding peakSqwaveBestImfBpOdav ridge 1 to hrOdavRidges.\n');
    case 2
        nextRidge = nextRidge + 1;
        hrOdavRidges(:,nextRidge) = ridgePeakSqwaveBestImfBpOdav(:,2);
        fprintf('Adding peakSqwaveBestImfBpOdav ridge 2 to hrOdavRidges.\n');
    otherwise
        fprintf('No peakSqwaveBestImfBpOdav ridge added to hrOdavRidges.\n');
end

%% Calculate mu & sigma for consensus ODBA, ODAV, and ODBA + ODAV

% first do the accel stuff...

meanHrOdba = mean(hrOdbaRidges, 'All');
stdHrOdba = std(hrOdbaRidges, 1);

consensusHrOdba = mean(hrOdbaRidges, 2);
consensusHrOdbaStd = mean(stdHrOdba);



% ... then do the gyro stuff...

meanHrOdav = mean(hrOdavRidges, 'All');
stdHrOdav = std(hrOdavRidges, 1);

consensusHrOdav = mean(hrOdavRidges, 2);
consensusHrOdavStd = mean(stdHrOdav);

% ... then put it all together ...

meanConsensusHr = mean( [consensusHrOdba, consensusHrOdav], 'All');
stdConsensusHr = mean( [consensusHrOdbaStd, consensusHrOdavStd] );

consensusHr = mean( [consensusHrOdba, consensusHrOdav], 2);


%% 

figure('Position', [100 400 800 800]);

p1a = subplot(311);
plot(kT, consensusHr .* 60, 'Color', maroon, 'LineWidth', 1.2);
xlabel('Time, local');
ylabel('Frequency, BPM');
title('Consensus Heart Rate Estimate (in BPM) from Kinematics');
grid;

p1b = subplot(312);
plot(kT, consensusHr, 'Color', blue, 'LineWidth', 1.2);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Consensus Heart Rate Estimate (in Hz) from Kinematics');
grid;

p1c = subplot(313);
plot(kT, consensusHrOdba, 'Color', red, 'LineWidth', 1.2);
hold on;
plot(kT, consensusHrOdav, 'Color', green, 'LineWidth', 1.2);
hold off;
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Heart Rate Estimate (in Hz) from Accel & Gyro');
grid;

linkaxes([p1a p1b p1c],'x');


%% Bland-Altman of final ODBA:ODAV synthesis

figSizeBA = [ 100 100 600 600 ];

figHrConsensus = figure('Position', figSizeBA, 'Color', white, 'NumberTitle', 'off');
[ rpcHrConsensus, figHrConsensus, stats1 ] = BlandAltman(figHrConsensus, ...
    consensusHrOdba, consensusHrOdav, ...
    {'hrOdba','hrOdav'}, ...
    'Bland-Altman test for ODBA vs ODAV consensus', ...
    {'ODBA','ODAV'}, ...
    'corrInfo', {'eq','r','r2','p','SSE','RMSE','n'}, ...
    'baStatsMode','Non-parametric');





%% 
% 
% figure('Position',[100 400 800 600]);
% 
% p2a = subplot(311);
% plot(kT, meanCandidateHrRidge, 'Color', blue, 'LineWidth', 1.2); 
% hold on; 
% plot(kT, meanCandidateHrRidge + (1.96 * sdCandidateHrRidge), ...
%     'Color', green); 
% plot(kT, meanCandidateHrRidge + (-1.96 * sdCandidateHrRidge), ...
%     'Color', green); 
% hold off;
% xlabel('Time, local');
% ylabel('Frequency, Hz');
% title('Consensus heart rate frequency from kinematics');
% grid;
% 
% p2b = subplot(312);
% plot(kT, meanCandidateOdbaRidge, 'Color', red, 'LineWidth', 1.2);
% hold on;
% plot(kT, meanCandidateOdbaRidge + (1.96 * sdCandidateOdbaRidge), ...
%     '-.', 'Color', maroon);
% plot(kT, meanCandidateOdbaRidge - (1.96 * sdCandidateOdbaRidge), ...
%     '-.', 'Color', maroon);
% hold off;
% xlabel('Time, local');
% ylabel('Frequency, Hz');
% title('Candidate HR Frequency from accelerometer data');
% grid;
% 
% p2c = subplot(313);
% plot(kT, meanCandidateOdavRidge, 'Color', green, 'LineWidth', 1.2);
% hold on;
% plot(kT, meanCandidateOdavRidge + (1.96 * sdCandidateOdavRidge), ...
%     '-.', 'Color', cyan);
% plot(kT, meanCandidateOdavRidge - (1.96 * sdCandidateOdavRidge), ...
%     '-.', 'Color', cyan);
% hold off;
% xlabel('Time, local');
% ylabel('Frequency, Hz');
% title('Candidate HR Frequency from gyroscope data');
% grid;
% 
% linkaxes([p2a p2b p2c],'x');
% 


%%
% 
% figure;
% plot(kT, meanCandidateOdbaRidge, 'Color', green);
% hold on;
% plot(kT, meanCandidateOdbaRidge + (1.96 * sdCandidateOdbaRidge), ...
%     '--', 'Color', blue);
% plot(kT, meanCandidateOdbaRidge - (1.96 * sdCandidateOdbaRidge), ...
%     '--', 'Color', blue);
% hold off;
% xlabel('Time, local');
% ylabel('Frequency, Hz');
% title('Candidate HR Frequency from gyroscopy');
% grid;


%% polyfit ridge data
% 
% widthCandidateHrRidge = width(candidateHrRidge);
% 
% [p,err,mu] = polyfit(kT_s, candidateHrRidge(:,1:widthCandidateHrRidge)', 9);
% 
% [p,err,mu] = polyfit(kT_s, meanCandidateHrRidge', 9);
% 
% [ridgeFit, deltaFit] = polyval(p,kT_s,err,mu);
% 
% figure; 
% plot(kT_s, meanCandidateHrRidge, 'Color', black);
% hold on;
% plot(kT_s, ridgeFit, 'b--');
% plot(kT_s, ridgeFit + 2 * deltaFit, 'r--');
% plot(kT_s, ridgeFit + -2 * deltaFit, 'r--');
% hold off;
% grid;
% 
% opt = fitoptions('Method','NonlinearLeastSquares', ...
%     'Startpoint',[1,0.2]);
% 
% f = fittype('p1*cos(p2*x)+p2*sin(p1*x)','options', opt);
% 
% fitobject = fit(kT_s, ridges1(:,1), f);


%% make proposals for freqHr and timeHr from first and second ridges

approveForSave = lower(input('Approve these candidates? (y/n): ','s'));

switch(approveForSave)
    
    case 'y'
        
        % hrOdbaRidges <- has all accel-based consensus candidates
        % hrOdavRidges <- has all gyro-based consensus candidates
        % consensusHrOdba <- has the consensus estimate of hrOdbaRidges
        % consensusHrOdav <- has the consensus estimate of hrOdavRidges
        % consensusHr <- has consensus estimate using ODBA + ODAV

        hrTimeOdbaRidges    = hrOdbaRidges .* 60;
        hrTimeOdavRidges    = hrOdavRidges .* 60;
        hrTimeConsensusOdba = consensusHrOdba .* 60;
        hrTimeConsensusOdav = consensusHrOdav .* 60;
        hrTimeConsensus     = consensusHr .* 60;

    otherwise
        fprintf('Halting. Go back and fix whatever is needed.\n');
        return;
        
end
        




%%

% [curve, goodness, output] = fit(kT_s, ridgeBpOdba(:,1), 'smoothingspline');

[curve, goodness, output] = fit(kT_s, consensusHr, 'smoothingspline');

figure; 
plot(curve, kT_s, consensusHr, 'residuals' );
xlabel('Time, seconds');
ylabel('Frequency');

figure; 
plot( consensusHr, output.residuals, 'r.' );
xlabel('ConsensusHr');
ylabel('Residuals');


%% Ask user if they are satisfied with analysis results for save to RAW

prepToSave = lower(input('Make HR timetables and peak structs? (y/n): ','s'));

switch(prepToSave)
    
    case 'y'

        fprintf('Building HR timetable and peak structs...\n');
        
        Time = kT;

        % Build a generic timetable with all necessary constituent parts

        TT = timetable(Time, bpOdba, bpOdav, imfBpOdba, imfBpOdav, ...
            sumImfBpOdba, sumImfBpOdav, bestImfBpOdba, bestImfBpOdav, ...
            threshSqwaveBestImfBpOdba, threshSqwaveAntiBestImfBpOdba, ...
            threshSqwaveBestImfBpOdav, threshSqwaveAntiBestImfBpOdav, ...
            peaksBestImfBpOdba, peaksBestImfBpOdav, ...
            hrOdbaRidges, hrOdavRidges, consensusHrOdba, ...
            consensusHrOdav, consensusHr, ...
            hrTimeOdbaRidges, hrTimeOdavRidges, hrTimeConsensusOdba, ...
            hrTimeConsensusOdav, hrTimeConsensus);

        % Build the correct timetable and structs based on epochChoice

        switch(epochChoice)

            case 1

                fprintf('Creating E1 timetable and peak structures...\n');
                E1_KINEMATICS = TT;
                E1_PEAKS_ODBA = HR_ODBA;
                E1_PEAKS_ODAV = HR_ODAV;


            case 2

                fprintf('Creating E2 timetable and peak structures...\n');
                E2_KINEMATICS = TT;
                E2_PEAKS_ODBA = HR_ODBA;
                E2_PEAKS_ODAV = HR_ODAV;


            case 3

                fprintf('Creating E3 timetable and peak structures...\n');
                E3_KINEMATICS = TT;
                E3_PEAKS_ODBA = HR_ODBA;
                E3_PEAKS_ODAV = HR_ODAV;

            case 4

                fprintf('Creating E4 timetable and peak structures...\n');
                E4_KINEMATICS = TT;
                E4_PEAKS_ODBA = HR_ODBA;
                E4_PEAKS_ODAV = HR_ODAV;


            case 5

                fprintf('Creating E5 timetable and peak structures...\n');
                E5_KINEMATICS = TT;
                E5_PEAKS_ODBA = HR_ODBA;
                E5_PEAKS_ODAV = HR_ODAV;

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




