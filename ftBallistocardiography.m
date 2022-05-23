%%  ftBallistocardiography.m

%   Written by Dave Haas on 2 October

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

makePlot = false;


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

%% Do accelerometer accelerometry ballistocardiography estimate

% set dT

dT = downsample(kT, 2);

% take 100 Hz band-passed A* and decimate to 50 Hz
dBpAx = decimate(bpAx, 2);
dBpAy = decimate(bpAy, 2);
dBpAz = decimate(bpAz, 2);

% calculate the square root of the square of dBpA*
rssDBpAx = sqrt(dBpAx.^2);
rssDBpAy = sqrt(dBpAy.^2);
rssDBpAz = sqrt(dBpAz.^2);

% set the threshold factor for s standard deviations for clipping

stdThr = 2.5;


% cap rssDBpAx at mean + 5·SD, so vents don't overwhelm signal

rssDBpAx(rssDBpAx >= (mean(rssDBpAx) + ...
    stdThr * (std(rssDBpAx)))) = ...
    (mean(rssDBpAx)+ stdThr * (std(rssDBpAx)));

rssDBpAy(rssDBpAy >= (mean(rssDBpAy) + ...
    stdThr * (std(rssDBpAy)))) = ...
    (mean(rssDBpAy) + stdThr * (std(rssDBpAy)));

rssDBpAz(rssDBpAz >= (mean(rssDBpAz) + ...
    stdThr * (std(rssDBpAz)))) = ...
    (mean(rssDBpAz) + stdThr * (std(rssDBpAz)));

figure;

r1 = subplot(311);
plot(dT, rssDBpAx, 'Color', blue); grid;

r2 = subplot(312);
plot(dT, rssDBpAy, 'Color', red); grid;

r3 = subplot(313);
plot(dT, rssDBpAz, 'Color', goldenrod); grid;

sumRss = rssDBpAx + rssDBpAy + rssDBpAz;

figure; 
plot(dT, sumRss);
grid;

[imf] = vmd(sumRss, 'NumIMFs', 10);

bestImf = imf(:,9);

figure; 
plot(dT, bestImf);
grid;

figure;
[wsstThis, fThis] = wsst(bestImf, 50);
pcolor(dT, fThis, abs(wsstThis));
shading interp;
ylim([0 4]);
title('bestImf WSST');

[ridge, iridge] = wsstridge(wsstThis, 20, fThis, 'NumRidges', 2);

figure;
plot(dT, ridge);
grid;
title('bestImf ridge estimates');