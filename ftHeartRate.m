%%  ftHeartRate.m

%   Written by Dave Haas between 9-19 July 2021
%   Updated by Dave Haas on 26 July 2021 to accept new raw data format

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
if ( exist(fullRawFileName, 'file') )
    fprintf('Raw file located for this trial... loading it!\n');
    load(fullRawFileName);
    
    % confirm presence of OPTICS variable
    if ( exist('OPTICS','var') )
        
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
hold off;
xlabel('Time, local');
ylabel('Intensity (A.U.)');
grid;
legend('1050 nm','1200 nm');

p1b = subplot(312);
plot(KINEMATICS.Time, KINEMATICS.odba, 'Color', green);
xlabel('Time, local');
ylabel('m/s^2');
title('ODBA');
grid;

p1c = subplot(313);
plot(PRESSURE.Time, PRESSURE.depth, 'Color', blue');
xlabel('Time, local');
ylabel('Depth, meters');
title('Tag Depth');
set(gca,'YDir','reverse');
grid;

linkaxes([p1a p1b p1c],'x');

analysisRegionDefined = false;

while(~analysisRegionDefined)

    fprintf('Click once to set the left-edge of the analysis region.\n');
    [aStart, ~, ~] = ginput(1);

    fprintf('Click once to set the right-edge of the analysis region.\n');
    [aEnd, ~, ~] = ginput(1);

    xAxis = gca;
    xStart = num2ruler(aStart, xAxis.XAxis);
    xEnd = num2ruler(aEnd, xAxis.XAxis);

    analysisRange = timerange(xStart, xEnd);

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
        otherwise
            analysisRegionDefined = false;
    end
    
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

%% Apply some denoise tests...

oLength = numel(oLed2);
kLength = numel(kAx);
okRatio = round( (oLength / kLength), 2);

if ( okRatio ~= 2.5000 )
   fprintf('Time series mismatch. oLength / kLength = %1.2f.\n', ...
       okRatio);
   fprintf('Halt and investigate before proceeding.\n');
   return;
else
    fprintf('Time series matches. O:K ratio = %1.2f.\n', okRatio);
end

% led1 = zeros(1,oLength);
% led2 = zeros(1,oLength);
% led3 = zeros(1,oLength);
% led4 = zeros(1,oLength);
% 
% led1(1,:) = OPTICS.led1;
% led2(1,:) = OPTICS.led2;
% led3(1,:) = OPTICS.led3;
% led4(1,:) = OPTICS.led4;

[dnLed1,~,thrParams,~,bestNbOfInt] = cmddenoise(oLed1,'db3',3);
[dnLed2,~,thrParams,~,bestNbOfInt] = cmddenoise(oLed2,'db3',3);
[dnLed3,~,thrParams,~,bestNbOfInt] = cmddenoise(oLed3,'db3',3);
[dnLed4,~,thrParams,~,bestNbOfInt] = cmddenoise(oLed4,'db3',3);

[dnRLed1,~,thrParams,~,bestNbOfInt] = cmddenoise(rLed1,'db3',3);
[dnRLed2,~,thrParams,~,bestNbOfInt] = cmddenoise(rLed2,'db3',3);
[dnRLed3,~,thrParams,~,bestNbOfInt] = cmddenoise(rLed3,'db3',3);
[dnRLed4,~,thrParams,~,bestNbOfInt] = cmddenoise(rLed4,'db3',3);

oHampelWin = oFs / 25;      % 10 + 1 samples for optical window
oSigma = 2;                 % +/- 2 standard deviations

[hdnLed1, o1OutlierIndex, o1Median, o1Std] = hampel(dnLed1, oHampelWin, oSigma);
[hdnLed2, o2OutlierIndex, o2Median, o2Std] = hampel(dnLed2, oHampelWin, oSigma);
[hdnLed3, o3OutlierIndex, o3Median, o3Std] = hampel(dnLed3, oHampelWin, oSigma);
[hdnLed4, o4OutlierIndex, o4Median, o4Std] = hampel(dnLed4, oHampelWin, oSigma);

% bpHdnLed1 = filtfilt(bandPassHrOpticalFilter', hdnLed1);
% bpHdnLed2 = filtfilt(bandPassHrOpticalFilter', hdnLed2);
% bpHdnLed3 = filtfilt(bandPassHrOpticalFilter', hdnLed3);
% bpHdnLed4 = filtfilt(bandPassHrOpticalFilter', hdnLed4);
% 
% fprintf('Doing a quick WSST for sgrHdnLed2+3...\n');
% tic
% [wsstBpHdnLed2, fBpHdnLed2] = wsst(bpHdnLed2, oFs);
% [wsstBpHdnLed3, fBpHdnLed3] = wsst(bpHdnLed3, oFs);
% toc;
% 
% figHdnLed23 = figure;
% 
% p2a = subplot(211);
% pcolor(oT, fBpHdnLed2, abs(wsstBpHdnLed2));
% % ylim([0 4]);
% shading interp;
% xlabel('Time, local');
% ylabel('Frequency, Hz');
% title('1050 nm');
% 
% p2b = subplot(212);
% pcolor(oT, fBpHdnLed3, abs(wsstBpHdnLed3));
% % ylim([0 4]);
% shading interp;
% xlabel('Time, local');
% ylabel('Frequency, Hz');
% title('1050 nm');
% 
% % The following three blocks don't produce usable WSSTs
% % sqrDnLed1 = dnLed1 .^ 2;
% % sqrDnLed2 = dnLed2 .^ 2;
% % sqrDnLed3 = dnLed3 .^ 2;
% % sqrDnLed4 = dnLed4 .^ 2;
% % 
% % [hdnSqrLed1, o1OutlierIndex, o1Median, o1Std] = hampel(sgrDnLed1, oHampelWin, oSigma);
% % [hdnSqrLed2, o2OutlierIndex, o2Median, o2Std] = hampel(sqrDnLed2, oHampelWin, oSigma);
% % [hdnSqrLed3, o3OutlierIndex, o3Median, o3Std] = hampel(sqrDnLed3, oHampelWin, oSigma);
% % [hdnSqrLed4, o4OutlierIndex, o4Median, o4Std] = hampel(sqrDnLed4, oHampelWin, oSigma);
% % 
% % bpHSqLed1 = filtfilt(bandPassHrOpticalFilter', hdnSqrLed1);
% % bpHSqLed2 = filtfilt(bandPassHrOpticalFilter', hdnSqrLed2);
% % bpHSqLed3 = filtfilt(bandPassHrOpticalFilter', hdnSqrLed3);
% % bpHSqLed4 = filtfilt(bandPassHrOpticalFilter', hdnSqrLed4);
% % 
% % fprintf('Doing a quick WSST for sgrHdnLed2+3...\n');
% % tic
% % [wsstBpHSqLed2, fBpHSqLed2] = wsst(bpHSqLed2, oFs);
% % [wsstBpHSqLed3, fBpHSqLed3] = wsst(bpHSqLed3, oFs);
% % toc;
% % 
% % figHdnLed23 = figure;
% % 
% % p2a = subplot(211);
% % pcolor(oT, fBpHSqLed2, abs(wsstBpHSqLed2));
% % ylim([0 4]);
% % shading interp;
% % xlabel('Time, local');
% % ylabel('Frequency, Hz');
% % title('1050 nm');
% % 
% % p2b = subplot(212);
% % pcolor(oT, fBpHSqLed3, abs(wsstBpHSqLed3));
% % ylim([0 4]);
% % shading interp;
% % xlabel('Time, local');
% % ylabel('Frequency, Hz');
% % title('1050 nm');
% 
% figure;
% subplot(221);
% plot(oT,hdnSqrLed1);grid;
% subplot(222);
% plot(oT,hdnSqrLed2);grid;
% subplot(223);
% plot(oT,hdnSqrLed3);grid;
% subplot(224);
% plot(oT,hdnSqrLed4);grid;
% 
% 
% sgolayLed1 = sgolayfilt(hdnLed1, 20, 401);
% sgolayLed2 = sgolayfilt(hdnLed2, 20, 401);
% sgolayLed3 = sgolayfilt(hdnLed3, 20, 401);
% sgolayLed4 = sgolayfilt(hdnLed4, 20, 401);
% 
% sgrHdnLed1 = hdnLed1 - sgolayLed1;
% sgrHdnLed2 = hdnLed2 - sgolayLed2;
% sgrHdnLed3 = hdnLed3 - sgolayLed3;
% sgrHdnLed4 = hdnLed4 - sgolayLed4;
% 
% figure;
% subplot(221);
% plot(oT,hdnLed1,oT,sgrHdnLed1);grid;
% subplot(222);
% plot(oT,hdnLed2,oT,sgrHdnLed2);grid;
% subplot(223);
% plot(oT,hdnLed3,oT,sgrHdnLed3);grid;
% subplot(224);
% plot(oT,hdnLed4,oT,sgrHdnLed4);grid;
% 
% sqrSgrHdnLed1 = (sgrHdnLed1 .^ 2);
% sqrSgrHdnLed2 = (sgrHdnLed2 .^ 2);
% sqrSgrHdnLed3 = (sgrHdnLed3 .^ 2);
% sqrSgrHdnLed4 = (sgrHdnLed4 .^ 2);


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

% gfdLed1 = filtfilt(gfBo, gfA, denoiseLed1);
% gfdLed2 = filtfilt(gfBo, gfA, denoiseLed2);
% gfdLed3 = filtfilt(gfBo, gfA, denoiseLed3);
% gfdLed4 = filtfilt(gfBo, gfA, denoiseLed4);

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

makePlot = lower(input('Want to see Gaussian filter plots? (y/n): ','s'));
switch(makePlot)
    case 'y'
        plotIt = true;
    otherwise
        plotIt = false;
end

if (plotIt)
    
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

closeFig = lower(input('Close Gaussian window filter plot? (y/n): ','s'));
switch(closeFig)
    case 'y'
         if (exist('figGaussKin','var'))
            close(figGaussOpt);
        end
    otherwise
end

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

closeFig = lower(input('Close Gaussian window filter plot? (y/n): ','s'));
switch(closeFig)
    case 'y'
        if (exist('figGaussKin','var'))
            close(figGaussKin);
        end
    otherwise
end

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
[AkDCHP, BkDCHP, CkDCHP, DkDCHP] = ...
    butter(dcoFiltOrder, ( dcoPassbandFreqHz / (kNyquistFs) ) );

[oSOS,oG] = ss2sos(AoDCHP, BoDCHP, CoDCHP, DoDCHP);

% design the high-pass filter

% dcoOpticalFilter with using oNyquistFs as samplingRate for filter
%
% dcoOpticalFilter = designfilt('highpassiir','FilterOrder',dcoFiltOrder, ...
%     'PassbandFrequency', dcoPassbandFreqHz, ... 
%     'PassbandRipple', 0.2, ...
%     'SampleRate', oNyquistFs);

% dcoOpticalFilter using samplingRate = 2, which is the designfilt default
% resulting in normalized frequencies
dcoOpticalFilter = designfilt('highpassiir','FilterOrder',dcoFiltOrder, ...
    'PassbandFrequency', dcoPassbandFreqHz, ... 
    'PassbandRipple', 0.2, ...
    'SampleRate', 2);

% % dcoKinematicFilter using samplingRate of kNyquistFs
% dcoKinematicFilter = designfilt('highpassiir','FilterOrder',dcoFiltOrder, ...
%     'PassbandFrequency', dcoPassbandFreqHz, ... 
%     'PassbandRipple', 0.2, ...
%     'SampleRate', kNyquistFs);

% dcoKinematicFilter using samplingRate = 2, which is the designfilt default
% resulting in normalized frequencies
dcoKinematicFilter = designfilt('highpassiir','FilterOrder',dcoFiltOrder, ...
    'PassbandFrequency', dcoPassbandFreqHz, ... 
    'PassbandRipple', 0.2, ...
    'SampleRate', 2);

% take a look at these
fvt = fvtool(oSOS, dcoOpticalFilter, 'Fs', oFs);
legend(fvt,'butter','designfilt');

dcoLed1 = filtfilt(dcoOpticalFilter', led1);
dcoLed2 = filtfilt(dcoOpticalFilter', led2);
dcoLed3 = filtfilt(dcoOpticalFilter', led3);
dcoLed4 = filtfilt(dcoOpticalFilter', led4);
% 
% dcodLed1 = filtfilt(dcoOpticalFilter', denoiseLed1);
% dcodLed2 = filtfilt(dcoOpticalFilter', denoiseLed2);
% dcodLed3 = filtfilt(dcoOpticalFilter', denoiseLed3);
% dcodLed4 = filtfilt(dcoOpticalFilter', denoiseLed4);
 
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

makePlot = lower(input('Want to see DCO/HP 0.01 Hz filter plots? (y/n): ','s'));
switch(makePlot)
    case 'y'
        plotIt = true;
    otherwise
        plotIt = false;
end

if (plotIt)
    
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

closeFig = lower(input('Close DC-removal filter plot? (y/n): ','s'));
switch(closeFig)
    case 'y'
        if (exist('figDco','var'))
            close(figDco);
        end
    otherwise
end

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

sos = ss2sos(Aolp,Bolp,Colp,Dolp);
fvt = fvtool(sos, lowPassOpticalFilter, 'Fs', oNyquistFs);
legend(fvt,'butter','designfilt');

% design the high-pass filter for kinematics
lowPassKinematicFilter = designfilt('lowpassiir','FilterOrder', lpFiltOrder, ...
    'PassbandFrequency', kPassbandFreqHz, ... 
    'PassbandRipple', 0.2, ...
    'SampleRate', kNyquistFs);

sos = ss2sos(Aklp, Bklp, Cklp, Dklp);
fvt = fvtool(sos, lowPassKinematicFilter, 'Fs', kNyquistFs);
legend(fvt,'butter','designfilt');
 
% proceedStr = lower(input('Are these good filters (y/n) [y]: ','s'));
% if (isempty(proceedStr))
%     proceedStr = 'y';
% end
% this filter works great, so just proceed, no need to look at it
proceedStr = 'y';

switch (proceedStr)
    case 'y'
        
        lpLed1 = filtfilt(lowPassOpticalFilter', dcoLed1);
        lpLed2 = filtfilt(lowPassOpticalFilter', dcoLed2);
        lpLed3 = filtfilt(lowPassOpticalFilter', dcoLed3);
        lpLed4 = filtfilt(lowPassOpticalFilter', dcoLed4);
        
%         lpdLed1 = filtfilt(lowPassOpticalFilter', dcodLed1);
%         lpdLed2 = filtfilt(lowPassOpticalFilter', dcodLed2);
%         lpdLed3 = filtfilt(lowPassOpticalFilter', dcodLed3);
%         lpdLed4 = filtfilt(lowPassOpticalFilter', dcodLed4);
        
        lpAx = filtfilt(lowPassKinematicFilter', dcoAx);
        lpAy = filtfilt(lowPassKinematicFilter', dcoAy);
        lpAz = filtfilt(lowPassKinematicFilter', dcoAz);
        lpGx = filtfilt(lowPassKinematicFilter', dcoGx);
        lpGy = filtfilt(lowPassKinematicFilter', dcoGy);
        lpGz = filtfilt(lowPassKinematicFilter', dcoGz);
        
        lpOdba = sqrt( lpAx .^ 2 + lpAy .^ 2 + lpAz .^ 2);
        lpOdav = sqrt( lpGx .^ 2 + lpGy .^ 2 + lpGz .^ 2);
        
        %close(fvt);
        
    otherwise
        fprintf('Re-run this section after specifying new filter params.\n');
        return;
        
end

%% plot the differences after low-pass (high-frequency) removal

makePlot = lower(input('Want to see low-pass filter plots? (y/n): ','s'));
switch(makePlot)
    case 'y'
        plotIt = true;
    otherwise
        plotIt = false;
end

if (plotIt)
    
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

closeFig = lower(input('Close low-pass (high frequency) filter plot? (y/n): ','s'));
switch(closeFig)
    case 'y'
        if (exist('figLowPass','var'))
            close(figLowPass);
        end
    otherwise
end

end         % end if (plotIt == true)

%% Try using a 30th order, 401 frame Sovitzy-Golay filter, where
%       a) separate the slow-varying motion wander
%       b) leave the out the more rapidly varying heart signal component

fprintf('Sovitzy-Golay filtering, 30th order, 401-frame.\n');
fprintf('\tsgoLed* = Sovitzy-Golay filter of raw OPTICS.led*\n');
fprintf('\tsgoA*/G* = Sovitzy-Golay filter of raw KINEMATICS.*\n');
fprintf('\tsgfLed* = Sovitzy-Golay filter of filtered Led* (-DC +LP optics)\n');
fprintf('\tsgfA*/G* = Sovitzy-Golay filter of filtered A*/G* (-DC +LP kinematics)\n');
fprintf('\tsgr* = Sovitzy-Golay filter subtraction from raw signals\n');

sgoLed1 = sgolayfilt(led1, 30, 401);
sgoLed2 = sgolayfilt(led2, 30, 401);
sgoLed3 = sgolayfilt(led3, 30, 401);
sgoLed4 = sgolayfilt(led4, 30, 401);

% sgodLed1 = sgolayfilt(denoiseLed1, 30, 401);
% sgodLed2 = sgolayfilt(denoiseLed2, 30, 401);
% sgodLed3 = sgolayfilt(denoiseLed3, 30, 401);
% sgodLed4 = sgolayfilt(denoiseLed4, 30, 401);

% high-pass (DC-removed) + low-pass optical signals with an SG filter

sgfLed1 = sgolayfilt(lpLed1, 30, 401);
sgfLed2 = sgolayfilt(lpLed2, 30, 401);
sgfLed3 = sgolayfilt(lpLed3, 30, 401);
sgfLed4 = sgolayfilt(lpLed4, 30, 401);

% sgfdLed1 = sgolayfilt(lpdLed1, 30, 401);
% sgfdLed2 = sgolayfilt(lpdLed2, 30, 401);
% sgfdLed3 = sgolayfilt(lpdLed3, 30, 401);
% sgfdLed4 = sgolayfilt(lpdLed4, 30, 401);

% SG-removed optics

sgrLed1 = led1 - sgoLed1;
sgrLed2 = led2 - sgoLed2;
sgrLed3 = led3 - sgoLed3;
sgrLed4 = led4 - sgoLed4;

sgfrLed1 = lpLed1 - sgfLed1;
sgfrLed2 = lpLed2 - sgfLed2;
sgfrLed3 = lpLed3 - sgfLed3;
sgfrLed4 = lpLed4 - sgfLed4;



% sgrdLed1 = denoiseLed1 - sgodLed1;
% sgrdLed2 = denoiseLed2 - sgodLed2;
% sgrdLed3 = denoiseLed3 - sgodLed3;
% sgrdLed4 = denoiseLed4 - sgodLed4;


sgoAx = sgolayfilt(kAx, 30, 401);
sgoAy = sgolayfilt(kAy, 30, 401);
sgoAz = sgolayfilt(kAz, 30, 401);
sgoGx = sgolayfilt(kGx, 30, 401);
sgoGy = sgolayfilt(kGy, 30, 401);
sgoGz = sgolayfilt(kGz, 30, 401);

% now do the kinematics

sgfAx = sgolayfilt(lpAx, 30, 401);
sgfAy = sgolayfilt(lpAy, 30, 401);
sgfAz = sgolayfilt(lpAz, 30, 401);
sgfGx = sgolayfilt(lpGx, 30, 401);
sgfGy = sgolayfilt(lpGy, 30, 401);
sgfGz = sgolayfilt(lpGz, 30, 401);

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

makePlot = lower(input('Want to see Savitzky-Golay filter plots? (y/n): ','s'));
switch(makePlot)
    case 'y'
        plotIt = true;
    otherwise
        plotIt = false;
end

if (plotIt)
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

closeFig = lower(input('Close Sovitzy-Golay filter plot? (y/n): ','s'));
switch(closeFig)
    case 'y'
        if (exist('figSgolay','var'))
            close(figSgolay);
        end
    otherwise
end

end         % end if (plotIt == true)

%% plot the differences after low-pass (high-frequency) removal

makePlot = lower(input('Want to see filter difference plots? (y/n): ','s'));
switch(makePlot)
    case 'y'
        plotIt = true;
    otherwise
        plotIt = false;
end

if (plotIt)
    
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

% %% Plot sgSubtraction from raw A* and G* signals
% 
% figure;
% 
% plot(oT_s, hrLed1,'k');
% hold on;
% plot(oT_s, hrLed2,'b');
% plot(oT_s, hrLed3,'r');
% plot(oT_s, hrLed4,'g');
% xlabel('Time, seconds');
% ylabel('Optical intensity');
% title('Optics after HP+LP+SG subtraction');
% legend('Ambient1','1050 nm','1200 nm','Ambient2');
% grid;
% 
% 
% figure; 
% 
% s4a = subplot(211);
% plot(kT_s, hrAx, kT_s, hrAy, kT_s, hrAz);
% xlabel('Time, seconds');
% ylabel('A* - sgA*, m/s^2');
% title('Accelerometry HR estimator');
% legend('ax', 'ay', 'az');
% grid;
% 
% s4b = subplot(212);
% 
% plot(kT_s, hrGx, kT_s, hrGy, kT_s, hrGz); grid;
% xlabel('Time, seconds');
% ylabel('G* - sgG*, °/s');
% title('Gyroscopy HR estimator');
% grid;
% legend('gx', 'gy', 'gz');
% grid;
% 
% linkaxes([s4a s4b],'x');


closeFig = lower(input('Close DCO vs. SGO filter plot? (y/n): ','s'));
switch(closeFig)
    case 'y'
        if (exist('figDifferences','var'))
            close(figDifferences);
        end
    otherwise
end

end         % end if (plotIt == true)

%% Band-pass only components between 0.333 Hz and 3.5 Hz

bandpassFiltOrder = 6;        % 6th order butterworth filter

oLowCut            = 0.3333;   % 0.3333 Hz = ~20 beats per minute
oHighCut           = 3.5;      % 3.5 Hz    = 210 beats per minute
opticalTestFs      = 8;        % 8 produces big pulses around odba HR in bpLed2                      

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

fprintf('Band-pass butterworth filter for seismocardiography energies:\n');
fprintf('\tFilter order: %d\n', bandpassFiltOrder);
fprintf('\tLow cut frequency: %d Hz\n', kLowCut);
fprintf('\tHigh cut frequency: %d Hz\n', kHighCut);

[Akbp, Bkbp, Ckbp, Dkbp] = butter(bandpassFiltOrder/2, ...
    [kLowCut kHighCut] / (kNyquistFs) );

bandPassHrOpticalFilter = designfilt('bandpassiir','FilterOrder',bandpassFiltOrder, ...
    'HalfPowerFrequency1',oLowCut,'HalfPowerFrequency2', oHighCut, ...
	'SampleRate', opticalTestFs);

bandPassHrKinematicFilter = designfilt('bandpassiir','FilterOrder',bandpassFiltOrder, ...
    'HalfPowerFrequency1', kLowCut,'HalfPowerFrequency2', kHighCut, ...
	'SampleRate', kNyquistFs);

[sosBpo,gBpo] = ss2sos(Aobp,Bobp,Cobp,Dobp);
fvtO = fvtool(sosBpo, bandPassHrOpticalFilter, 'Fs', opticalTestFs);
legend(fvtO,'butter','designfilt');

% [sosBpk,gBpk] = ss2sos(Akbp,Bkbp,Ckbp,Dkbp);
% fvtK = fvtool(sosBpk, kHeartRateFilter, 'Fs', kNyquistFs);
% legend(fvtK,'butter','designfilt');

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

% bpSgrdLed1 = filtfilt(bandPassHrOpticalFilter', sgrdLed1);
% bpSgrdLed2 = filtfilt(bandPassHrOpticalFilter', sgrdLed2);
% bpSgrdLed3 = filtfilt(bandPassHrOpticalFilter', sgrdLed3);
% bpSgrdLed4 = filtfilt(bandPassHrOpticalFilter', sgrdLed4);

bpAx   = filtfilt(bandPassHrKinematicFilter',kAx);
bpAy   = filtfilt(bandPassHrKinematicFilter',kAy);
bpAz   = filtfilt(bandPassHrKinematicFilter',kAz);
bpGx   = filtfilt(bandPassHrKinematicFilter',kGx);
bpGy   = filtfilt(bandPassHrKinematicFilter',kGy);
bpGz   = filtfilt(bandPassHrKinematicFilter',kGz);

bpOdba = sqrt( bpAx.^2 + bpAy.^2 + bpAz.^2);
bpOdav = sqrt( bpGx.^2 + bpGy.^2 + bpGz.^2);

% odav comes in now with kOdav
% odav = sqrt(kGx.^2 + kGy.^2 + kGz.^2);
% dcoOdav = filtfilt(dcoKinematicFilter', odav);
% lpOdav = filtfilt(lowPassKinematicFilter', odav);
% dcoLpOdav = filtfilt(dcoKinematicFilter', lpOdav);
% bpOdavAlt = filtfilt(bandPassHrKinematicFilter', odav);






%% Now all the heart rate investigation variables should be in place...

% make some bpLed plots as a contrast to the bpSgrLed plots

makePlot = lower(input('Want to see band-pass filter plots? (y/n): ','s'));
switch(makePlot)
    case 'y'
        plotIt = true;
    otherwise
        plotIt = false;
end

if (plotIt) 
    
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

figBpKin = figure; 

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


%% Okay, everything is in place for heart rate analysis...

%  Pause here and let the user select which analysis to run

fprintf('\n\n');
fprintf('Signal processing is complete on this trial data.\n');
fprintf('Run wsst/wsstridge & pulse-square-wave manually...\n');
fprintf('This is currently starting around line 1225.\n');

waitForInput = lower(input('Press [enter] to begin your analyses...','s'));
return;


%% wsst for some of these

% for tt21_141a, 100 seconds into the trial is a great start time
% for tt21_142d, 563 seconds into the trial is a great start time

% entire time series, or just a time slice?

useEntireTimeSeries = lower(input('Analyze entire time series? (y/n): ','s'));
switch(useEntireTimeSeries)
    case 'y'
        oStart = oT(1);
        oEnd = oT(end);
        oS = find(oT == oStart);
        oE = find(oT == oEnd);
        
        kStart = kT(1);
        kEnd = kT(end);
        kS = find(kT == kStart);
        kE = find(kT == kEnd);
    
    otherwise
        
        % timeStart = 563;    % good pulse data on tt21_142d
        % timeStart = 110;    % good pulse data on tt21_141a
        timeStart = str2double('Enter time in seconds to look at: ','s');
        
        %timeSlice = 30;     % look at 30 second windows
        timeSlice = str2double('Number of seconds for timeSlice: ','s');
        
        oS = oFs * timeStart;
        oE = oS + (oFs * timeSlice);

        kS = kFs * timeStart;
        kE = kS + (kFs * timeSlice);

end


%% do kinematic WSSTs (or not)

% specify which kinematics to use for this...

% testAx = sgrAx;
% testAy = sgrAy;
% testAz = sgrAz;
% testGx = sgrGx;
% testGy = sgrGy;
% testGz = sgrGz;

testAx = bpAx;
testAy = bpAy;
testAz = bpAz;
testGx = bpGx;
testGy = bpGy;
testGz = bpGz;

testOdba = bpOdba;
testOdav = bpOdav;

% check out lpGx

doKinematicWssts = lower(input('WSST for accelXYZ & gyroXYZ? (y/n): ','s'));

if (strcmp(doKinematicWssts,'y'))

    fprintf('Starting wsst compute for accelXYZ...\n');
    tic
    [wsstAx,fAx] = wsst(testAx(kS:kE), kFs);
    [wsstAy,fAy] = wsst(testAy(kS:kE), kFs);
    [wsstAz,fAz] = wsst(testAz(kS:kE), kFs);
    fprintf('Completed wsst compute for accelXYZ:\n');
    toc
    
    fprintf('Starting wsst compute for gyroXYZ...\n');
    [wsstGx,fGx] = wsst(testGx(kS:kE), kFs);
    [wsstGy,fGy] = wsst(testGy(kS:kE), kFs);
    [wsstGz,fGz] = wsst(testGz(kS:kE), kFs);
    fprintf('Completed wsst compute for gyroXYZ:\n');
    toc
    
    FigKinematicWssts = figure;

    wsstAccelX = subplot(321);
    pcolor(kT(kS:kE), fAx, abs(wsstAx) );
    shading interp
    xlabel('Seconds');
    ylabel('Frequency (Hz)');
    title('Synchrosqueezed Transform - Low-passed Accel X');
    ylim([0 4]);

    wsstGyroX = subplot(322);
    pcolor(kT(kS:kE), fGx, abs(wsstGx) );
    shading interp
    xlabel('Seconds');
    ylabel('Frequency (Hz)');
    title('Synchrosqueezed Transform - Low-passed Gyro X');
    ylim([0 4]);

    wsstAccelY = subplot(323);
    pcolor(kT(kS:kE), fAy, abs(wsstAy) );
    shading interp
    xlabel('Seconds');
    ylabel('Frequency (Hz)');
    title('Synchrosqueezed Transform - Low-passed Accel Y');
    ylim([0 4]);

    wsstGyroY = subplot(324);
    pcolor(kT(kS:kE), fGy, abs(wsstGy) );
    shading interp
    xlabel('Seconds');
    ylabel('Frequency (Hz)');
    title('Synchrosqueezed Transform - Low-passed Gyro Y');
    ylim([0 4]);

    wsstAccelZ = subplot(325);
    pcolor(kT(kS:kE), fAz, abs(wsstAz) );
    shading interp
    xlabel('Seconds');
    ylabel('Frequency (Hz)');
    title('Synchrosqueezed Transform - Low-passed Accel Z');
    ylim([0 4]);

    wsstGyroZ = subplot(326);
    pcolor(kT(kS:kE), fGz, abs(wsstGz) );
    shading interp
    xlabel('Seconds');
    ylabel('Frequency (Hz)');
    title('Synchrosqueezed Transform - Low-passed Gyro Z');
    ylim([0 4]);

    linkaxes( [wsstAccelX wsstGyroX wsstAccelY wsstGyroY wsstAccelZ wsstGyroZ], 'x');

end         % end for doKinematicsWssts


% Gaussian smooth ODBA and ODAV to make nice smooth peaks

% work backwards to get an integer divisble window size

kWinSizeDivisor = kFs / 10;                  % == 10, for a 1/10 second window
kGaussWindowSize = kFs / kWinSizeDivisor;   % = 0.100 s
gfBpKin = ones(1, kGaussWindowSize) / kGaussWindowSize;

gfTestOdba = filtfilt(gfBpKin, gfA, testOdba);
gfTestOdav = filtfilt(gfBpKin, gfA, testOdav);

figCompTestBpVsGfBp = figure;

pTestOdba1 = subplot(221);
plot(kT, testOdba); 
xlabel('Time, local');
ylabel('testOdba');
grid;
legend('testODBA');

pTestOdav1 = subplot(222);
plot(kT, testOdav); 
xlabel('Time, local');
ylabel('testOdav');
grid;
legend('testODAV');

pGfTestOdba = subplot(223);
plot(kT, gfTestOdba); 
xlabel('Time, local');
ylabel('GfTestOdba');
grid;
legend('gfBp ODBA');

pGfTestOdav = subplot(224);
plot(kT, gfTestOdav); 
xlabel('Time, local');
ylabel('GfTestOdav');
grid;
legend('gfBp ODAV');

linkaxes([pTestOdba1 pTestOdav1 pGfTestOdba pGfTestOdav],'x');

% prompt user to replace testOd* with gfTestOd* (if it improves signal)

replaceTestWithGfTest = lower(input('Replace test with gfTest? (y/n): ','s'));
if (strcmp(replaceTestWithGfTest,'y'))
    testOdba = gfTestOdba;
    testOdav = gfTestOdav;
end
    

doOdbaOdavWssts = lower(input('WSST for ODBA & ODAV? (y/n): ','s'));
if (strcmp(doOdbaOdavWssts,'y'))
    
    tic
    fprintf('Computing ODBA WSST...\n');
    [wsstOdba,fOdba] = wsst(testOdba(kS:kE), kFs); 
    fprintf('Completed ODBA WSST compute:\n');
    toc
    
    tic
    fprintf('Computing ODAV WSST...\n');
    [wsstOdav,fOdav] = wsst(testOdav(kS:kE), kFs);
    fprintf('Completed ODAV WSST compute:\n');
    toc

    FigOdbaOdavWssts = figure;

    pWsstOdba = subplot(211);
    pcolor(kT(kS:kE), fOdba, abs(wsstOdba) );
    shading interp;
    xlabel('Time, seconds');
    ylabel('Frequency, Hz');
    title('Synchrosqueezed Transform - band-passed ODBA');
    ylim([0 2]);

    pWsstOdav = subplot(212);
    pcolor(kT(kS:kE), fOdav, abs(wsstOdav) );
    shading interp;
    xlabel('Time, seconds');
    ylabel('Frequency, Hz');
    title('Synchrosqueezed Transform - band-passed ODAV');
    ylim([0 2]);

    linkaxes([pWsstOdba pWsstOdav ],'x');

end         % end doOdbaOdavWssts

% partition those ODBA and ODAV into frequencies 0-4 Hz only!

reduceTxt = 'Reduce ODBA & ODAV WSSTs to 0-4 Hz prior to ridging? (y/n): ';
bandPassOdbaOdavWssts = lower(input(reduceTxt,'s'));
if (strcmp(bandPassOdbaOdavWssts,'y'))

    hrOdbaUpperLimit = numel(find(fOdba < 4));
    hrOdbaF = fOdba(1:hrOdbaUpperLimit);
    hrOdbaWSST = wsstOdba(1:hrOdbaUpperLimit,:);
    
    hrOdavUpperLimit = numel(find(fOdav < 4));
    hrOdavF = fOdav(1:hrOdavUpperLimit);
    hrOdavWSST = wsstOdav(1:hrOdavUpperLimit,:);

    FigBpOdbaOdavWsst = figure;
    
    pBpOdbaWsst = subplot(211);
    pcolor(kT(kS:kE), hrOdbaF, abs(hrOdbaWSST) );
    shading interp;
    xlabel('Time, local');
    ylabel('Frequency, Hz');
    title('Synchrosqueezed Transform - band-passed ODBA');
    ylim([0 4]);
    
    pBpOdavWsst = subplot(212);
    pcolor(kT(kS:kE), hrOdavF, abs(hrOdavWSST) );
    shading interp;
    xlabel('Time, local');
    ylabel('Frequency, Hz');
    title('Synchrosqueezed Transform - band-passed ODAV');
    ylim([0 4]);
    
end         %	end 0-4 Hz band-pass of ODBA & ODAV WSST

findOdbaOdavWsstRidges = lower(input('Find 0-4 Hz WSST ridges? (y/n): ','s'));

if (strcmp(findOdbaOdavWsstRidges,'y'))

    [hrOdbaRidge, iOdavRidge] = wsstridge(hrOdbaWSST, 20, hrOdbaF, ...
        'NumRidges', 2);
    [hrOdavRidge, iOdavRidge] = wsstridge(hrOdavWSST, 20, hrOdavF, ...
        'NumRidges', 2);

    figOdbaOdavWsstRidges = figure;

    pHrOdba1 = subplot(221);
    plot(kT(kS:kE), hrOdbaRidge(:,1), 'Color', blue);
    hold on;
    plot(kT(kS:kE), hrOdbaRidge(:,2), '--', 'Color', blue);
    hold off;
    xlabel('Time, local');
    ylabel('Frequency, Hz');
    xlim([ kT(kS) kT(kE)]);
    ylim([ 0 round(max(hrOdbaRidge(:,1))) ]);
    title('BpOdba Synchrosqueezed Wavelet Heart Rate Estimation, Hz');
    grid;    
    
    pHrOdav1 = subplot(222);
    plot(kT(kS:kE), hrOdavRidge(:,1), 'Color', red);
    hold on; 
    plot(kT(kS:kE), hrOdavRidge(:,2), '--', 'Color', red);
    hold off;
    xlabel('Time, local');
    ylabel('Frequency, Hz');
    xlim([ kT(kS) kT(kE)]);
    ylim([ 0 round(max(hrOdavRidge(:,1))) ]);
    title('BpOdav Synchrosqueezed Wavelet Heart Rate Estimation, Hz');
    grid;

    pHrOdba2 = subplot(223);
    plot(kT(kS:kE), hrOdbaRidge(:,1) * 60, 'Color', blue);
    hold on;
    plot(kT(kS:kE), hrOdbaRidge(:,2) * 60, '--', 'Color', blue);
    hold off;
    xlabel('Time, local');
    ylabel('Frequency, BPM');
    xlim([ kT(kS) kT(kE) ]);
    ylim([ 0 round(max(hrOdbaRidge(:,1)) * 60) ]);
    title('BpOdba Synchrosqueezed Wavelet Heart Rate Estimation, BPM');
    grid;
    
    pHrOdav2 = subplot(224);
    plot(kT(kS:kE), hrOdavRidge(:,1) * 60, 'Color', red);
    hold on;
    plot(kT(kS:kE), hrOdavRidge(:,2) * 60, '--', 'Color', red);
    hold off;
    xlabel('Time, local');
    ylabel('Frequency, BPM');
    xlim([ kT(kS) kT(kE)] );
    ylim([ 0 round(max(hrOdavRidge(:,1)) * 60) ]);
    title('BpOdav Synchrosqueezed Wavelet Heart Rate Estimation, BPM');
    grid;

    linkaxes([ pHrOdba1 pHrOdav1 pHrOdba2 pHrOdav2],'x');

end     % end findOdbaOdavWsstRidges


%% Create some aggregate band-passed optics, then Gaussian filter them...

%   specify which optical treatment you want to use
%	options include:
%       1 = gfLed* (Gaussian smoothed @ 0.1s, raw reflectance-adjusted LEDs)
%       2 = dcoLed* (HP + 0.01 Hz, raw reflectance-adjusted LEDs)
%       3 = lpLed* (LP high-freq noise-removal of dcoLed*)
%       4 = sgoLed* (Savitzky-Golay filtered raw reflectance-adjusted LEDs)
%       5 = sgfLed* (Savitzky-Golay filtered lpLed*)
%       6 = sgrLed* (Savitzky-Golay filter LED* subtracted from raw LED*)
%       7 = sgfrLed* (Savitzky-Golay filter lpLed* subtracted from lpLed*)
%       8 = bpSgrLed* (band-passed SGR LED*)

fprintf('\n\n');
fprintf('WSST and WSST ridge finding for optical data:\n');
fprintf('Specify one of the optical data sets to use for WSST processing:\n');
fprintf('\t0 = led* (raw reflectance-adjusted LEDs)\n');
fprintf('\t1 = gfLed* (Gaussian smoothed @ 0.1s, raw reflectance-adjusted LEDs)\n');
fprintf('\t2 = dcoLed* (HP + 0.01 Hz, raw reflectance-adjusted LEDs)\n');
fprintf('\t3 = lpLed* (LP high-freq noise-removal of dcoLed*)\n');
fprintf('\t4 = sgoLed* (Savitzky-Golay filtered raw reflectance-adjusted LEDs)\n');
fprintf('\t5 = sgfLed* (Savitzky-Golay filtered lpLed*)\n');
fprintf('\t6 = sgrLed* (Savitzky-Golay filter LED* subtracted from raw LED*)\n');
fprintf('\t7 = sgfrLed* (Savitzky-Golay filter lpLed* subtracted from lpLed*)\n');
fprintf('\t8 = bpSgrLed* (band-passed SGR LED*)\n');

ledChoiceTxt = 'Enter the number for your optical data set choice: ';
opticalDataChoice = str2double(input(ledChoiceTxt,'s'));

switch(opticalDataChoice)
    case 0
        fprintf('Using led* (raw reflectance-adjusted optical data)\n');
        testLed1 = led1;
        testLed2 = led2;
        testLed3 = led3;
        testLed4 = led4;
        ledChoice = 'led*';    
    case 1
        fprintf('Using gfLed*\n');
        testLed1 = gfLed1;
        testLed2 = gfLed2;
        testLed3 = gfLed3;
        testLed4 = gfLed4;
        ledChoice = 'gfLed*';
    case 2
        fprintf('Using dcoLed*\n');
        testLed1 = dcoLed1;
        testLed2 = dcoLed2;
        testLed3 = dcoLed3;
        testLed4 = dcoLed4;
        ledChoice = 'dcoLed*';
    case 3
        fprintf('Using lpLed*\n');
        testLed1 = lpLed1;
        testLed2 = lpLed2;
        testLed3 = lpLed3;
        testLed4 = lpLed4;
        ledChoice = 'lpLed*';
    case 4
        fprintf('Using sgoLed*\n');
        testLed1 = sgoLed1;
        testLed2 = sgoLed2;
        testLed3 = sgoLed3;
        testLed4 = sgoLed4;
        ledChoice = 'sgoLed*';
    case 5
        fprintf('Using sgfLed*\n');
        testLed1 = sgfLed1;
        testLed2 = sgfLed2;
        testLed3 = sgfLed3;
        testLed4 = sgfLed4;
        ledChoice = 'sgfLed*';
    case 6
        fprintf('Using sgrLed*\n');
        testLed1 = sgrLed1;
        testLed2 = sgrLed2;
        testLed3 = sgrLed3;
        testLed4 = sgrLed4;
        ledChoice = 'sgrLed*';
    case 7
        fprintf('Using sgfrLed*\n');
        testLed1 = sgfrLed1;
        testLed2 = sgfrLed2;
        testLed3 = sgfrLed3;
        testLed4 = sgfrLed4;
        ledChoice = 'sgfrLed*';
    case 8
        fprintf('Using bpSgrLed*\n');
        testLed1 = bpSgrLed1;
        testLed2 = bpSgrLed2;
        testLed3 = bpSgrLed3;
        testLed4 = bpSgrLed4;
        ledChoice = 'bpSgrLed*';
    otherwise
        fprintf('Defaulting to bpSgrLed*\n');
        testLed1 = bpSgrLed1;
        testLed2 = bpSgrLed2;
        testLed3 = bpSgrLed3;
        testLed4 = bpSgrLed4;
        ledChoice = 'bpSgrLed*';
end

testLed23 = sqrt(testLed2 .^ 2 + testLed3 .^ 2);
testLedAll = sqrt(testLed1 .^ 2 + testLed2 .^ 2 + testLed3 .^ 2 + ...
    testLed4 .^ 2);

plotComboOptics = lower(input('Plot combo optics? (y/n): ','s'));
if(strcmp(plotComboOptics,'y'))
    
    figComboOptics = figure;
    pComboOpt1 = subplot(211);
    plot(oT, testLed23, 'Color', blue); 
    xlabel('Time, local');
    ylabel('Intensity (A.U.)');
    title('Active optical channels (1050 nm and 1200 nm only)');
    grid;
    legend(ledChoice);
    
    pComboOpt2 = subplot(212);
    plot(oT, testLedAll, 'Color', red);
    xlabel('Time, local');
    ylabel('Intensity (A.U.)');
    title('All optical channels (active plus ambient)');
    grid;
    legend(ledChoice);
    
    linkaxes([pComboOpt1 pComboOpt2],'x');
    
end

% Ask if user wants to Gaussian filter these to make nice peaks

gfComboTxt = 'Apply Gaussian filter to emphasize peaks? (y/n): ';
gaussFilterComboOptics = lower(input(gfComboTxt,'s'));
if (strcmp(gaussFilterComboOptics,'y'))

    % work backwards to get an integer divisble window sizer
    oWinSizeDivisor = oFs / 25;
    optGaussHrWindowSize = oFs / oWinSizeDivisor;      % = 0.300 s
    gfHrBo = ones(1, optGaussHrWindowSize) / optGaussHrWindowSize;

    %kinGaussianWindowSize = kFs / 10;    % = 0.100 s
    %gfBk = ones(1, kinGaussianWindowSize) / kinGaussianWindowSize;

    gfTestLed23 = filtfilt(gfHrBo, gfA, testLed23);
    gfTestLedAll = filtfilt(gfHrBo, gfA, testLedAll);

    pGfTestLedPlot = figure;
    
    pGfTestLed23 = subplot(211);
    plot(oT, gfTestLed23, 'Color', blue); 
    xlabel('Time, local');
    ylabel('Intensity (A.U.)');
    title('Gaussian smoothed testLedActive');
    grid;
    legend(ledChoice);
    
    pGfTestLedAll = subplot(212);
    plot(oT, gfTestLedAll, 'Color', red);
    xlabel('Time, local');
    ylabel('Intensity (A.U.)');
    title('Gaussian smoothed testLedAll');
    grid;
    legend(ledChoice);

end

doLed2 = lower(input('WSST for testLed2? (y/n): ','s'));
if (strcmp(doLed2,'y'))
    FigLed2Wsst = figure;
    [wsstLed2, fLed2] = wsst(testLed2(oS:oE), oFs);
    pcolor(oT(oS:oE), fLed2, abs(wsstLed2) );
    shading interp;
    ylim([0 4]);
    xlabel('Seconds');
    ylabel('Frequency (Hz)');
    title('Synchrosqueezed Transform SG-subtracted LED2');
end

doLed3 = lower(input('WSST for testLed3? (y/n): ','s'));
if (strcmp(doLed3,'y'))
    FigLed3Wsst = figure;
    [wsstLed3, fLed3] = wsst(testLed3(oS:oE), oFs);
    pcolor(oT(oS:oE), fLed3, abs(wsstLed3) );
    shading interp;
    ylim([0 4]);
    xlabel('Seconds');
    ylabel('Frequency (Hz)');
    title('Synchrosqueezed Transform SG-subtracted LED3');
end

% wsst of denoised Led2 with Savitzky-Golay filter removal

[wsstSgrdLed2,fSgrdLed2] = wsst(sgrdLed2,oFs);
fourHzCutOff = max(find(fSgrdLed2 < 4));
fNew = fSgrdLed2(1,1:fourHzCutOff);
wsstNew = wsstSgrdLed2(1:fourHzCutOff,:);

figure;
pcolor(oT, fNew, abs(wsstNew));
shading interp;

% doLed23;

% doLedAll;


%% Turn heart rate energy into simple square wave pulses

if (exist('oPulse','var'))
    % confirm numel of oPulse and kPulse equal numel of oT_s and kT_s
    if ( (numel(oPulse) == numel(oT_s) ) & (numel(kPulse) == numel(kT_s)))
        fprintf('oPulse and kPulse are correctly sized. Proceeding.\n');
    else
        fprintf('oPulse and kPulse are not correctly sized. Check this out.\n');
    end
else
    oPulse = zeros(1,numel(oT_s(oS:oE)));
    kPulse = zeros(1,numel(kT_s(kS:kE)));
end

% This should be deprecated now. Getting oS:oE and kS:kE from the
% conditional input above, i.e.: useEntireTimeSeries test after all the
% signal processing is completed.
%
% entireTimeSeries = true;
% 
% switch(entireTimeSeries)
%     case true 
%         fprintf('Using entire time series (probably a bad idea).\n');
%         oStart = 1;
%         oEnd = numel(oT_s);
%         kStart = 1;
%         kEnd = numel(kT_s);
%     case false
%         timeStart = 400;
%         windowSeconds = 10;
%         timeEnd = timeStart + windowSeconds;
%         fprintf('Starting at %d seconds using a %d second window.\n', ...
%             timeStart, windowSeconds);
%         oStart = find(oT_s == timeStart);
%         oEnd = find(oT_s == timeEnd);
%         kStart = find(kT_s == timeStart);
%         kEnd = find(kT_s == timeEnd);
% end

% specify the signal(s) to use for oTestSignal and kTestSignal

% oSignal = bpLed23;
oSignal = gfTestLed23;

%oTestSignal = normalize(abs(diff(sgrLed3(oStart:oEnd))));

% NOTE: using this method results in a slight linear trend that ends up
% making the signal wander away from a common axis around y = 0; it's
% probably worth trying to use bpLed23 and/or bpSgrLed23 by themselves, as
% they don't seem to have this linear trend moving away from zero by the
% end of the time series

% changed by Dave Haas on 26 July 2021 to experiment with bp and bpSgr
% a coarse aggregated sgrLed2 + sgrLed3
%
% oTestSignal = normalize(abs(diff(sgrLed2(oStart:oEnd)))) + ...
%    normalize(abs(diff(sgrLed3(oStart:oEnd))));
% gaussian sliding window filter to smooth
% oTestSignal = detrend(filtfilt(gfBo, gfA, oTestSignal));

oTestSignal = oSignal(oS:oE);


% sgrOdba
% kTestSignal = normalize(abs(diff(sgrOdba(kStart:kEnd))));
%
% a coarse aggregated ax + ay + az summation
%kTestSignal = normalize(abs(diff(sgrAx(kStart:kEnd)))) + ...
%    normalize(abs(diff(sgrAy(kStart:kEnd)))) + ...
%    normalize(abs(diff(sgrAz(kStart:kEnd))));
% gaussian sliding window filter to smooth
%kTestSignal = detrend(filtfilt(gfBk, gfA, kTestSignal));

kSignal = bpOdba;
%kSignal = gfBpOdav

kTestSignal = kSignal;

% make some easy time scales

oTimeComparator = numel(oS:oE) - numel(oTestSignal);
switch(oTimeComparator)
    case 0
        fprintf('Optical time scale and sample scale have equal numels.\n');
        oTimeSeries = oT_s(oS:oE);
    case 1
        fprintf('Adjusting optical time scale (using diff == 1).\n');
        oTimeSeries = oT_s(oS:oE-1);
    case 2
        fprintf('Adjusting optical time scale (using diff == 2).\n');
        oTimeSeries = oT_s(oS:oE-2);
    otherwise
        fprintf('Optical time scale mismatch is too large. Stopping...\n');
        return;
end

kTimeComparator = numel(kS:kE) - numel(kTestSignal);
switch(kTimeComparator)
    case 0
        fprintf('Kinematic time scale and sample scale have equal numels.\n');
        kTimeSeries = kT_s(kS:kE);
    case 1
        fprintf('Adjusting kinematic time scale (using diff == 1).\n');
        kTimeSeries = kT_s(kS:kE-1);
    case 2
        fprintf('Adjusting kinematic time scale (using diff == 2).\n');
        kTimeSeries = kT_s(kS:kE-2);
    otherwise
        fprintf('Kinematic time scale mismatch is too large. Stopping...\n');
        return;
end


% gets some quantiles to check thresholding
% 1: 50% || 2: 55% || 3: 60% || 4: 65% || 5: 70% || 6: 75% || 7: 80%
oQuants = quantile(oTestSignal,[0.50 0.55 0.60 0.65 0.70 0.75 0.80]);
kQuants = quantile(kTestSignal,[0.50 0.55 0.60 0.65 0.70 0.75 0.80]);

figure; 
subplot(211);
plot(oTimeSeries, oTestSignal);
hold on;
yline(oQuants(1), 'Color', black, 'LineWidth', 1.5);
yline(oQuants(2), 'Color', cyan, 'LineWidth', 1.5);
yline(oQuants(3), 'Color', maroon, 'LineWidth', 1.5);
yline(oQuants(4), 'Color', blue, 'LineWidth', 1.5);
yline(oQuants(5), 'Color', green, 'LineWidth', 1.5);
yline(oQuants(6), 'Color', goldenrod, 'LineWidth', 1.5);
yline(oQuants(7), 'Color', red, 'LineWidth', 1.5);
hold off;
xlabel('Time, local');
ylabel('Intensity');
title('oTestSignal');
grid;
legend('oTestSignal','50%','55%','60%','65%','70%','75%','80%');

subplot(212);
plot(kTimeSeries, kTestSignal);
hold on;
yline(kQuants(1), 'Color', black, 'LineWidth', 1.5);
yline(kQuants(2), 'Color', cyan, 'LineWidth', 1.5);
yline(kQuants(3), 'Color', maroon, 'LineWidth', 1.5);
yline(kQuants(4), 'Color', blue, 'LineWidth', 1.5);
yline(kQuants(5), 'Color', green, 'LineWidth', 1.5);
yline(kQuants(6), 'Color', goldenrod, 'LineWidth', 1.5);
yline(kQuants(7), 'Color', red, 'LineWidth', 1.5);
hold off;
xlabel('Time, local');
ylabel('Intensity');
title('kTestSignal');
grid;
legend('kTestSignal','50%','55%','60%','65%','70%','75%','80%');

fprintf('Choose an optical threshold for peak vs "noise"...\n');
fprintf('\t1 = 50th percentile\n');
fprintf('\t2 = 55th percentile\n');
fprintf('\t3 = 60th percentile\n');
fprintf('\t4 = 65th percentile\n');
fprintf('\t5 = 70th percentile\n');
fprintf('\t6 = 75th percentile\n');
fprintf('\t7 = 80th percentile\n');
oQuartChoiceText = 'Which optical quartile threshold do you want to use -> ';
oQuartChoice = str2double(input(oQuartChoiceText,'s'));

% 1: 50% || 2: 55% || 3: 60% || 4: 65% || 5: 70% || 6: 75% || 7: 80%
switch(oQuartChoice)
    case 1
        oLowSignal = oTestSignal < oQuants(1);
    case 2
        oLowSignal = oTestSignal < oQuants(2);
    case 3
        oLowSignal = oTestSignal < oQuants(3);
    case 4
        oLowSignal = oTestSignal < oQuants(4);
    case 5
        oLowSignal = oTestSignal < oQuants(5);
    case 6
        oLowSignal = oTestSignal < oQuants(6);
    case 7
        oLowSignal = oTestSignal < oQuants(7);
    otherwise
        fprintf('Defaulting to 75th percentile for thresholding\n');
        oLowSignal = oTestSignal < oQuants(2);
end


fprintf('Choose a kinematic threshold for peak vs "noise"...\n');
fprintf('\t1 = 50th percentile\n');
fprintf('\t2 = 55th percentile\n');
fprintf('\t3 = 60th percentile\n');
fprintf('\t4 = 65th percentile\n');
fprintf('\t5 = 70th percentile\n');
fprintf('\t6 = 75th percentile\n');
fprintf('\t7 = 80th percentile\n');
kQuartChoiceText = 'Which kinematic quartile threshold do you want to use -> ';
kQuartChoice = str2double(input(kQuartChoiceText,'s'));

% 1: 50% || 2: 55% || 3: 60% || 4: 65% || 5: 70% || 6: 75% || 7: 80%
switch(kQuartChoice)
    case 1
        kLowSignal = kTestSignal < kQuants(1);
    case 2
        kLowSignal = kTestSignal < kQuants(2);
    case 3
        kLowSignal = kTestSignal < kQuants(3);
    case 4
        kLowSignal = kTestSignal < kQuants(4);
    case 5
        kLowSignal = kTestSignal < kQuants(5);
    case 6
        kLowSignal = kTestSignal < kQuants(6);
    case 7
        kLowSignal = kTestSignal < kQuants(7);
    otherwise
        fprintf('Defaulting to 75th percentile for thresholding\n');
        kLowSignal = kTestSignal < kQuants(2);
end


oProposedSignal = ~bwareaopen(oLowSignal, (oFs/5) );
kProposedSignal = ~bwareaopen(kLowSignal, (kFs/5) ); 

figure;
subplot(211);
plot(oTimeSeries, oProposedSignal, 'b', 'LineWidth', 2); 
hold on; 
plot(oTimeSeries, oTestSignal ); 
hold off; 
grid;
subplot(212);
plot(kTimeSeries, kProposedSignal, 'b', 'LineWidth', 2); 
hold on; 
plot(kTimeSeries, kTestSignal ) ; 
hold off; 
grid;

[wsstO, fO] = wsst(double(oProposedSignal), oFs);
[wsstK, fK] = wsst(double(kProposedSignal), kFs);

figure;

p6a = subplot(211);
pcolor(oTimeSeries, fO, abs(wsstO));
shading interp;
ylim([0 4]);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Optical heart rate estimate - square-wave proposal');
grid;

p6b = subplot(212);
pcolor(kTimeSeries, fK, abs(wsstK));
shading interp;
ylim([0 4]);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('Kinematic (Gy) heart rate estimate - square-wave proposal');
grid;

linkaxes([p6a p6b],'x');

[hrRidgeO, iRidgeO] = wsstridge(wsstO, 20, fO, 'NumRidges', 3);
[hrRidgeK, iRidgeK] = wsstridge(wsstK, 10, fK, 'NumRidges', 3);

figure;

p7a = subplot(211);
plot(oTimeSeries, hrRidgeO(:,1) * 60, 'b', 'LineWidth', 1.2 );
hold on;
plot(oTimeSeries, hrRidgeO(:,2) * 60, 'b--', 'LineWidth', 1.2 );
hold off;
xlabel('Time, seconds');
ylabel('Heart rate, bpm');
title('Optical heart rate estimate - square-wave proposal');
grid;
legend('r1','r2');

p7b = subplot(212);
plot(kTimeSeries, hrRidgeK(:,1) * 60, 'r', 'LineWidth', 1.2 ); 
hold on;
plot(kTimeSeries, hrRidgeK(:,2) * 60, 'r--', 'LineWidth', 1.2 );
hold off;
xlabel('Time, seconds');
ylabel('Heart rate, bpm');
title('Kinematic heart rate estimate - square-wave proposal');
grid;
legend('r1','r2');

linkaxes([p7a p7b],'x');


% greatest common divisor calc for decimation factors

oDec = oFs / (gcd(oFs,kFs));
kDec = kFs / (gcd(oFs,kFs));

% fork out hrO and hrK from 1st ridges

oHR1 = ( decimate( hrRidgeO(:,1) * 60, oDec) )';
oHR2 = ( decimate( hrRidgeO(:,2) * 60, oDec) )';
oT_sDec = decimate(oT_s(oS:oE-1), oDec);

kHR1 = ( decimate( hrRidgeK(:,1) * 60, kDec) )';
kHR2 = ( decimate( hrRidgeK(:,2) * 60, kDec) )';
kT_sDec = decimate(kT_s(kStart:kEnd-1), kDec);

figure;
plot(oT_sDec, oHR1,'b');
hold on;
plot(oT_sDec, oHR2,'b--');
plot(kT_sDec, kHR1, 'r');
plot(kT_sDec, kHR2, 'r--');
hold off;
grid;

% add to oPulse and kPulse

addPulsesStr = lower(input('Add these to oPulse and kPulse (y/n) [y]: ','s'));
if (isempty(addPulsesStr))
    addPulsesStr = 'y';
end
switch(addPulsesStr)
    case 'y'
        oPulse(oStart:oEnd-1) = double(oProposedSignal);
        oPulse(end) = 0;
        kPulse(kStart:kEnd-1) = double(kProposedSignal); 
        kPulse(end) = 0;
    otherwise
        fprintf('Skipping addition of this set of pulses to oPulse and kPulse.\n');
end

%% do some bland altman plots of this

label = {'oHR','kHR'};
figTitle = 'Bland-Altman for oHR:kHR';
gnames = 'heart rate';
corrinfo = {'n','SSE','r2','eq'}; % stats to display of correlation scatter plot
% specify this next line for Gaussian distributed data
%BAinfo = {'RPC(%)','ks'}; % stats to display on Bland-ALtman plot
% specify this line for non-Gaussian distributed data
BAinfo = {'RPCnp','ks'};
limits = 'auto'; % how to set the axes limits
if 1 % colors for the data sets may be set as:
	colors = 'br';      % character codes
else
	colors = [0 0 1;... % or RGB triplets
		      1 0 0];
end
symbols = 'o';

BlandAltman(oHR1, kHR1, label, ...
    [figTitle ' (numbers, forced 0 intercept, fixed BA y-axis limits)'], ...
    gnames, 'corrInfo', corrinfo, 'axesLimits', limits, ...
    'colors', colors, 'symbols', 'Num', ...
    'baYLimMode', 'square', 'forceZeroIntercept', 'on')
    

BlandAltman(oHR1, kHR1, label, ...
    [figTitle ' (numbers, forced 0 intercept, fixed BA y-axis limits)'], ...
    gnames, 'corrInfo', corrinfo, 'axesLimits', limits, ...
    'colors', colors, 'symbols', 'Num', ...
    'baYLimMode', 'square', 'forceZeroIntercept', 'on',...
    'baStatsMode','non-parametric')

BlandAltman(oHR2, kHR2, label, ...
    [figTitle ' (numbers, forced 0 intercept, fixed BA y-axis limits)'], ...
    gnames, 'corrInfo', corrinfo, 'axesLimits', limits, ...
    'colors', colors, 'symbols', 'Num', ...
    'baYLimMode', 'square', 'forceZeroIntercept', 'on',...
    'baStatsMode','non-parametric')

%% investigate some CWT and other frequency analysis stuff



% Look at CWT and WSST for LED2 and LED3

% these show great pulsatile energy at ~10-12 Hz
cwt(sgrLed2(oStart:oEnd), oFs);    
cwt(sgrLed3(oStart:oEnd), oFs);

[wsstGfLed2, fGfLed2] = wsst(gfLed2(oStart:oEnd),oFs);
pcolor(opticsTime(oStart:oEnd), fGfLed2, abs(wsstGfLed2));
ylim([0 15]);
shading interp;


[wsstSgrLed2, fSgrLed2] = wsst(sgrLed2(oStart:oEnd),oFs);
pcolor(opticsTime(oStart:oEnd), fSgrLed2, abs(wsstSgrLed2))
ylim([0 15]);
shading interp;

[sgrLed3SST, sgrLed3F] = ...
    wsst(sgrLed3(oStart:oEnd), oFs);

% do CWT for sgrOdba and sgrOdav

figure;
cwt(sgrOdba(kStart:kEnd), kFs);
title('Savitsky-Golay removed ODBA');
grid;

figure;
cwt(sgrOdav(kStart:kEnd), kFs);
title('Savitsky-Golay removed ODAV');
grid;

% do -SG ODBA

[sgrOdbaSST, sgrOdbaF] = ...
    wsst(sgrOdba(kStart:kEnd), kFs);

[fSgrOdbaRidge,iSgrOdbaRidge] = ...
    wsstridge(sgrOdbaSST, 5, sgrOdbaF, 'NumRidges', 2);

figure;

subplot(211);
contour( kT_s(kStart:kEnd), sgrOdbaF, abs(sgrOdbaSST));
% ylim([0 20]);
grid;

subplot(212);
plot(kT_s(kStart:kEnd), sgrOdbaF, 'Color', kxColor, 'LineWidth', 1.2);
grid;


plot(kT_s(kStart:kEnd), fSgrOdbaRidge(:,2), ...
    'Color',kxColor,'LineWidth',1.2); 
xlabel('Time, seconds');
ylabel('Frequency, Hz');
title('WSST Savitsky-Golay subtracted ODBA');
grid;

% do -SG ODAV

[sgrOdavSST, sgrOdavF] = ...
    wsst(sgrOdav(kStart:kEnd), kFs);
[fSgrOdavRidge,iSgrOdavRidge] = ...
    wsstridge(sgrOdavSST, 5, sgrOdavF, 'NumRidges',2);

figure;

kWSST1 = subplot(221);
plot(kT_s(kStart:kEnd), fSgrOdbaRidge(:,1), ...
    'Color',kxColor,'LineWidth',1.2); 
xlabel('Time, seconds');
ylabel('Frequency, Hz');
title('WSST Savitsky-Golay subtracted ODBA');
grid;


kWSST2 = subplot(222);

grid;

kWSST3 = subplot(223);
plot(kT_s(kStart:kEnd), fSgrOdbaRidge(:,1)*60, ...
    'Color',kxColor,'LineWidth',1.2); 
xlabel('Time, seconds');
ylabel('Frequency, Hz');
title('WSST Savitsky-Golay subtracted ODBA');
grid;


kWSST4 = subplot(224);

grid;

%% Great plot of heart rate in tt21_134b 

% @ 178-192 seconds (free-breathe kStart:kEnd = 17800:19200)
%kStart = 17800; kEnd = 19200;

% @ 300-310 seconds (breath-hold kStart:kEnd = 30000:31000)
kStart = 29000; kEnd = 29500;

% @ 560-570 seconds (recover kStart:kEnd = 56000:57000)
% kStart = 56000; kEnd = 57000;

figOdbaOdav = figure; 

plotOdba = subplot(221); 
plot(kT_s(kStart:kEnd),kOdba(kStart:kEnd)); 
xlim([kT_s(kStart) kT_s(kEnd)]);
xlabel('Time, seconds');
ylabel('ODBA, m/s^2');
title('Overall Dynamic Body Acceleration');
grid; 

plotOdav = subplot(222);
plot(kT_s(kStart:kEnd),odav(kStart:kEnd)); 
xlim([kT_s(kStart) kT_s(kEnd)]);
xlabel('Time, seconds');
ylabel('ODAV, °/s');
title('Overall Dynamic Angular Velocity');
grid; 

plotWsstOdba = subplot(223); 
[wsstOdba, fOdba] = wsst(kOdba(kStart:kEnd), kFs);
pcolor(kT_s(kStart:kEnd), fOdba, abs(wsstOdba));
xlabel('Time, seconds');
ylabel('Frequency, Hz');
title('Synchrosqueezed wavelet transform · ODBA');
ylim([0 50]);
shading interp;

plotWsstOdav = subplot(224);
[wsstOdav, fOdav] = wsst(odav(kStart:kEnd), kFs);
pcolor(kT_s(kStart:kEnd), fOdav, abs(wsstOdav));
xlabel('Time, seconds');
ylabel('Frequency, Hz');
title('Synchrosqueezed wavelet transform · ODAV');
ylim([0 50]);
shading interp;

linkaxes([plotOdba, plotOdav, plotWsstOdba, plotWsstOdav],'x');

[fOdbaRidge, iOdbaRidge] = wsstridge(wsstOdba, fOdba, 'NumRidges', 4); 
[fOdavRidge, iOdavRidge] = wsstridge(wsstOdav, fOdav, 'NumRidges', 4);

figure;
plot(kT_s(kStart:kEnd),fOdbaRidge(:,1),'b'); 
hold on;
plot(kT_s(kStart:kEnd),fOdbaRidge(:,2),'r');
plot(kT_s(kStart:kEnd),fOdbaRidge(:,3),'g');
plot(kT_s(kStart:kEnd),fOdbaRidge(:,4),'k');
hold off;
grid;
xlim([kT_s(kStart) kT_s(kEnd)]);
xlabel('Time, seconds');
ylabel('Frequency, Hz');
title('WSST Ridges · ODBA');
legend('Ridge1', 'Ridge2');

figure;
plot(kT_s(kStart:kEnd),fOdavRidge(:,1),'b'); 
hold on;
plot(kT_s(kStart:kEnd),fOdavRidge(:,2),'r');
plot(kT_s(kStart:kEnd),fOdavRidge(:,3),'g');
plot(kT_s(kStart:kEnd),fOdavRidge(:,4),'k');
hold off;
grid;
xlim([kT_s(kStart) kT_s(kEnd)]);
xlabel('Time, seconds');
ylabel('Frequency, Hz');
title('WSST Ridges · ODAV');
legend('Ridge1', 'Ridge2');


%% investigate power spectrum of band-passed -SG LED2 and LED3

startTime = 563;        % time_s = 563 seconds

oStart = find(oT_s == startTime);
oEnd = oStart + (10 * oFs);

kStart = find(KINEMATICS.time_s == startTime);
kEnd = kStart + (10 * kFs);

bpSgrLed2 = filtfilt(heartRateFilter', sgrLed2);
bpSgrLed3 = filtfilt(heartRateFilter', sgrLed3);

[powerSpecLed2, fLed2] = pspectrum(gfLed2(oStart:oEnd), oFs);
[powerSpecLed3, fLed3] = pspectrum(gfLed3(oStart:oEnd), oFs);
[powerSpecGy, fGy] = pspectrum(testGy(kStart:kEnd), kFs);

figure;
subplot(211);
[wsstGfLed2, wsstFLed2] = wsst(gfLed2(oStart:oEnd),oFs);
pcolor(opticsTime(oStart:oEnd),wsstFLed2,abs(wsstGfLed2))
ylim([0 15]);
shading interp;
subplot(212);
[wsstGfLed3, wsstFLed3] = wsst(gfLed3(oStart:oEnd),oFs);
pcolor(opticsTime(oStart:oEnd),wsstFLed3,abs(wsstGfLed3))
ylim([0 15]);
shading interp;


pcolor(opticsTime(oStart:oEnd),wsstFLed2,abs(wsstGfLed2))
ylim([0 4]);
shading interp;

figure;
subplot(311);
plot(fLed2, pow2db(powerSpecLed2));
xlabel('Frequency, Hz');
ylabel('Power, dB');
title('-DC+HP Led2 1050 nm');
xlim([0 5]);
grid;
subplot(312);
plot(fLed3, pow2db(powerSpecLed3));
xlabel('Frequency, Hz');
ylabel('Power, dB');
title('-DC+HP Led3 1050 nm');
xlim([0 5]);
grid;
subplot(313);
plot(fGy, pow2db(powerSpecGy));
xlabel('Frequency, Hz');
ylabel('Power, dB');
title('-DC+HP Led3 1050 nm');
xlim([0 5]);
grid;

figure;
plot(fGy,pow2db(powerSpecGy))
grid on
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Default Frequency Resolution')

figure;
pspectrum(sgrGy(1,kStart:(kStart+1000)), kFs, 'spectrogram', ...
    'TimeResolution', 0.01);
grid;

figure;
subplot(211)
instfreq(hrLed2(1,140750:143250), oFs, 'Method','hilbert');
subplot(212);
instfreq(hrLed3(1,140750:143250), oFs, 'Method','hilbert');


%% Do some FFT analysis

fftWinStart = 6;
fftWinStop = 12;
fftSize = 128;  


for n = 1 : ( fix(length(bpOdba) / (2 * kFs)) - 2)

    bpOdbaSlice = bpOdba( (n * kFs):(n * kFs + fftSize - 1) );   
    hannWindow = hann(size(bpOdbaSlice, 2))';            % Hann window
    bpOdbaFFT = abs(fft(bpOdbaSlice .* hannWindow));     % FFT calculation
    
end

% calculate the single sided positive frequency
fftAx = abs(bpOdbaFFT(fftWinStart:fftWinStop));

%%

% [gyPeaks, gyLocs, gyWidths, gyProms] = findpeaks(sgrGy, kFs, ...
%     'MinPeakDistance',0.50);

% [gyPeaks, gyLocs, gyWidths, gyProms] = findpeaks(sgrGy, kT_s, ...
%     'MinPeakDistance', 0.50);

% this works for producing good indices for doing kT_s plots w/ gyLocs
[gyPeaks, gyLocs, gyWidths, gyProms] = findpeaks(sgrGy, ...
    'MinPeakDistance', 50);

figure;
plot(kT, sgrGy, 'Color', red);
hold on;
plot(kT(gyLocs), gyPeaks, 'd', 'Color', blue, 'MarkerFace', blue);
hold off;
xlabel('Time, local');
ylabel('°/s');
title('kGy -DC+LP-SGF');
grid;

dLevel = 3;
iLpLed2 = lpLed2 * -1;
ihLpLed2 = hampel(iLpLed2);
dIhLed2 = zeros(1,numel(lpLed2));
dIhLed2(1:end-dLevel) = cmddenoise( diff(ihLpLed2,dLevel), 'db3', 3);
figure; plot(oT, dIhLed2); grid;

iSgrLed2 = sgrLed2 * -1;
ihSgrLed2 = hampel(iSgrLed2);
dIhSgrLed2 = zeros(1,numel(sgrLed2));
dIhSgrLed2(1:end-dLevel) = cmddenoise( diff(ihSgrLed2,dLevel), 'db3', 3);
figure; plot(oT, dIhSgrLed2); grid;

% Gaussian filter the optical signal, to help minimize peak confusion
oWinSizeDivisor = oFs / 25;
optGaussHrWindowSize = oFs / oWinSizeDivisor;      % = 0.300 s
gfHrBo = ones(1, optGaussHrWindowSize) / optGaussHrWindowSize;
gfLpLed2 = filtfilt(gfHrBo, gfA, dIhLed2);
figure; plot(oT, gfLpLed2); grid;

figure;findpeaks(gfLpLed2,'MinPeakDistance',125);
[o2Peaks, o2Locs, o2Widths, o2Proms] = findpeaks(gfLpLed2, ...
    'MinPeakDistance', 125);

figure;

pHrLed2 = subplot(211);
plot(oT, lpLed2,'Color',blue);
hold on;
plot(oT, gfLpLed2, 'Color', red);
plot(oT(o2Locs), o2Peaks, 'd', 'Color', goldenrod, 'MarkerFace', goldenrod);
hold off;
xlabel('Time, local');
ylabel('Intensity (A.U.)');
title('1050 nm -DC+LP+GF');
grid;
legend('lpLed2','gfLpLed2','peaks');

pHrGy = subplot(212);
plot(kT, sgrGy, 'Color', red);
hold on;
plot(kT(gyLocs), gyPeaks, 'd', 'Color', blue, 'MarkerFace', blue);
hold off;
xlabel('Time, local');
ylabel('°/s');
title('kGy -DC+LP-SGF');
grid;

linkaxes([pHrLed2 pHrGy],'x');


figure;
findpeaks(d2IhLed2, 'MinPeakDistance', 125);

[o2Peaks, o2Locs, o2Widths, o2Proms] = findpeaks(d2IhLed2, ...
    'MinPeakDistance', 125);



%% End

fprintf('Ending ftHeartRate.m.\n');