%% FaunaLabs dataAnalysis
%
%       Created and written by Dave Haas, 18 June 2021
%       Most recent version: 0.5 on 29 June 2021
%
%       Based on work from FaunaTagDataWorkup.m by Dave Haas 
%           created on 17 June 2020
%           last updated on 18 June 2021
%       
%       FaunaData_Analysis fuses pulseOx/ppgPlayground.m (also written by 
%       Dave Haas in 2017) with figure plotting, peak-finding for heart
%       rate detection and calculation, frequency analysis tools, and other
%       analytics methods for describing FaunaTag data. This uses a
%       synthesis of bio-optical spatially-resolved reflectance data
%       collected by the FaunaTag's near-infrared sensor package as well as
%       kinematic signals (both raw and synthetic, e.g.: ODBA, to detect
%       patterns in cardio-pulmonary functioning of live animals)

clc;
clear;

%%  Load global TAG_PATHS, specify tag data set, and load RAW & PRH files

global TAG_PATHS;

if (exist(TAG_PATHS.RAW,'dir'))
    fprintf('RAW data folder specified in TAG_PATHS is present.\n');
else
    fprintf('RAW data folder specified in TAG_PATHS is not present.\n');
    fprintf('Fix missing RAW data folder and retry.\n');
    return;
end

if (exist(TAG_PATHS.PRH,'dir'))
    fprintf('PRH data folder specified in TAG_PATHS is present.\n');
else
    fprintf('PRH data folder specified in TAG_PATHS is not present.\n');
    fprintf('Fix missing PRH data folder and retry.\n');
    return;
end

if (exist(TAG_PATHS.METADATA,'dir'))
    fprintf('Meta data folder specified in TAG_PATHS is present.\n');
else
    fprintf('Meta data folder specified in TAG_PATHS is not present.\n');
    fprintf('Fix missing meta data folder and retry.\n');
    return;
end

%% Get tag data file of interest and read in the RAW and PRH data

tagStr = input('Which tag file do you want to work with? ','s');

rawTargetFile = sprintf('%s/%sraw.mat', TAG_PATHS.RAW, tagStr);
prhTargetFile = sprintf('%s/%sprh.mat', TAG_PATHS.PRH, tagStr);

if (exist(rawTargetFile,'file'))
    fprintf('RAW file %s is present. Loading here for analysis...\n', ...
        tagStr);
    load(rawTargetFile);
else
    fprintf('RAW file %s is not present at %s. Check and retry.\n', ...
        tagStr);
    return;
end

if (exist(prhTargetFile,'file'))
    fprintf('PRH file %s is present. Loading here for analysis...\n', ...
        tagStr);
    load(prhTargetFile);
else
    fprintf('PRH file %s is not present at %s. Check and retry.\n', ...
        tagStr);
    return;
end



%% Place TIME struct values into time scales commonly used in analysis

time = TIME.time;
time_s = TIME.time_s;

%% Always have a look at raw optical data before doing anything

figure;

rawLed1 = subplot(411);
plot(OPTICS.time, OPTICS.led1, 'k-');
xlabel('Time (local)');
ylabel('Intensity');
title('LED1 (amb1)');
grid;

rawLed2 = subplot(412);
plot(OPTICS.time, OPTICS.led2, 'b-');
xlabel('Time, seconds');
ylabel('Intensity');
title('LED2 (\lambda = 1050nm)');
grid;

rawLed3 = subplot(413);
plot(OPTICS.time, OPTICS.led3, 'r-');
xlabel('Time (local)');
ylabel('LED3 (1200nm) Intensity');
title('LED3 (\lambda = 1200 nm)');
grid;

rawLed4 = subplot(414);
plot(OPTICS.time, OPTICS.led4, 'g-');
xlabel('Time (local)');
ylabel('Intensity');
title('LED4 (amb2)');
grid;

linkaxes([ rawLed1, rawLed2, rawLed3, rawLed4], 'x');

%% plot some log plots of active LEDs...

figure;

optPlot1 = subplot(311);
plot(OPTICS.time, OPTICS.log10Led2, 'b-'); 
hold on; 
plot(OPTICS.time, OPTICS.log10Led3, 'r-'); 
hold off;
xlabel('Time (local)');
ylabel('Log10 Active LEDs')
title('Log10 LED2 (\lambda = 1050 nm) and LED3 (\lambda = 1200 nm)');
grid;

optPlot2 = subplot(312);
plot(OPTICS.time, OPTICS.log10Active * -1, 'k-');
xlabel('Time (local)');
ylabel('Tissue saturation index');
title('Log10 LED2 (\lambda = 1050 nm) / Log10 LED3 (\lambda = 1200 nm)');
grid;

optPlot3 = subplot(313);
plot(OPTICS.time, OPTICS.log10Ambient * -1, 'g-');
hold off;
xlabel('Time (local)');
ylabel('Log10 ambients');
title('Log10 Ambient Optical Channels');
grid;

linkaxes([optPlot1 optPlot2 optPlot3], 'x');

%% Also have a look at pitch roll & heading... 

% should probably do tag2whale on this somehow automatically

figure;

figPitch = subplot(311);
plot(KINEMATICS.time, KINEMATICS.pitch, 'b-');
xlabel('Time (local)');
ylabel('Pitch, °');
grid;

figRoll = subplot(312);
plot(KINEMATICS.time, KINEMATICS.roll, 'r-');
xlabel('Time (local)');
ylabel('Roll, °');
grid;

figHeading = subplot(313);
plot(KINEMATICS.time, KINEMATICS.heading, 'k-');
xlabel('Time (local)');
ylabel('Heading, °');
grid;

linkaxes([ figPitch, figRoll, figHeading], 'x');

%% overplot undecimated norm'd rawLed2 & rawLed3

figure;

plot(OPTICS.time, OPTICS.nrLed1, 'k-');
hold on;
plot(OPTICS.time, OPTICS.nrLed2, 'b-');
plot(OPTICS.time, OPTICS.nrLed3, 'r-');
plot(OPTICS.time, OPTICS.nrLed4, 'g-');
hold off;
xlabel('Time, seconds');
ylabel('Normalized reflectance LED Intensities');
grid;


%% Plot some decimated data to know where to manually trim out bad data

figure;

dec1 = subplot(411);
plot(time, O.rLed1, 'k-');
hold on;
plot(time, O.rLed2, 'b-');
plot(time, O.rLed3, 'r-');
plot(time, O.rLed4, 'g-');
hold off;
xlabel('Time (local)');
ylabel('Reflectance Intensity');
grid;
title('Decimated optical data');

dec2 = subplot(412);
plot(time, K.pitch, 'b-');
hold on;
plot(time, K.roll, 'r-');
plot(time, K.heading, 'g-');
xlabel('Time, seconds');
ylabel('Degrees');
title('Pitch · Roll · Heading');
ylim([-200 200]);
grid;
hold off;

dec3 = subplot(413);
plot(time, K.normOdba, 'k-');
xlabel('Time (local)');
ylabel('ODBA');
grid;

dec4 = subplot(414);
plot(time, P.depth, 'b-');
ylabel('Depth, meters');
xlabel('Time (local)');
set(gca, 'YDir','reverse');
grid;

linkaxes([dec1 dec2 dec3 dec4], 'x');

%% Resize sensor data to get rid of pre-tag-on + post-tag-off

%   TO-DO: think about inserting a function here to trim off the noisy data
%   bits, e.g.: the data immediately before and immediately after the tag
%   is placed on and taken off the animal (where sunlight overwhelms the
%   signals).


%% normalize decimated accel & gyro data

normAx   = normalize(K.ax);
normAy   = normalize(K.ay);
normAz   = normalize(K.az);

normGx   = normalize(K.gx);
normGy   = normalize(K.gy);
normGz   = normalize(K.gz);

figure;

normKin1 = subplot(211);
plot(time_s, normAx, 'b-');
hold on;
plot(time_s, normAy, 'r-');
plot(time_s, normAz, 'g-');
hold off;
xlabel('Time, seconds');
ylabel('Normalized accel');
title('Norm accel');
legend('ax','ay','az');
grid;

normKin2 = subplot(212);
plot(time_s, normGx, 'b-');
hold on;
plot(time_s, normGy, 'r-');
plot(time_s, normGz, 'g-');
hold off;
xlabel('Time, seconds');
ylabel('Normalized gyro');
title('Norm gyro');
legend('gx','gy','gz');
grid;

linkaxes([ normKin1 normKin2 ], 'x');

%% Is there an ODBA-like version for gyroscope?

%   there is now... call it Overall Dynamic Gyroscopic Signal (ODGS)

odav = sqrt(KINEMATICS.gx.^2 + KINEMATICS.gy.^2 + KINEMATICS.gz.^2);
normOdav = normalize(KINEMATICS.gx.^2 + KINEMATICS.gy.^2 + KINEMATICS.gz.^2);

bwFiltOrder = 4;        % specify the order butterworth filter

%lowCut      = 0.2;
%lowCut      = 0.3333;   % 0.3333 Hz = ~20 beats per minute
lowCut      = 0.5;      % 0.5 Hz    = ~30 beats per minute
highCut     = 3.5;      % 3.5 Hz    = 210 beats per minute

testFs = 8;             % 8 produces big pulses around odba HR in bpLed2
% testFs = 9;             % 9 produces big pulses around odba HR in bpLed2
% testFs = 25;            % 25 produces some overlaps in bpLed3 but washes
                          % out the bpLed2 signal. Maybe this suggests
                          % using separate filters for LED2 and LED3?                      
nyquistFs   = testFs;        

fprintf('\t\tFilter order: %d\n', bwFiltOrder);
fprintf('\t\tLow cut frequency: %d Hz | low HR: %d BPM\n', ...
    lowCut, (lowCut/(1/60)) );
fprintf('\t\tHigh cut frequency: %d Hz | high HR: %d BPM\n', ...
    highCut, (highCut/(1/60)) );

[A, B, C, D] = butter(bwFiltOrder/2, [lowCut highCut] / (nyquistFs/2) );

d = designfilt('bandpassiir','FilterOrder',bwFiltOrder, ...
    'HalfPowerFrequency1',lowCut,'HalfPowerFrequency2',highCut, ...
	'SampleRate', nyquistFs);

sos = ss2sos(A,B,C,D);
fvt = fvtool(sos, d, 'Fs', nyquistFs);
legend(fvt,'butter','designfilt');

%% If you like that filter, then use it here!

KINEMATICS.bpAz = filtfilt(d, KINEMATICS.az);
KINEMATICS.bpGx = filtfilt(d, KINEMATICS.gx);
KINEMATICS.bpGy = filtfilt(d, KINEMATICS.gy);
KINEMATICS.bpGz = filtfilt(d, KINEMATICS.gz);

figure;
plot(KINEMATICS.time_s, normalize(KINEMATICS.bpAz), 'k-');
hold on;
plot(KINEMATICS.time_s, normalize(KINEMATICS.bpGx), 'b-');
plot(KINEMATICS.time_s, normalize(KINEMATICS.bpGy), 'r-');
plot(KINEMATICS.time_s, normalize(KINEMATICS.bpGz), 'g-');
hold off;
xlabel('Time, seconds');
ylabel('normBp.K');
title('Normalized band-passed kinematics');
legend('Az','Gx','Gy','Gz');
grid;

K.bpAz = filtfilt(d, K.az);
K.bpGx = filtfilt(d, K.gx);
K.bpGy = filtfilt(d, K.gy);
K.bpGz = filtfilt(d, K.gz);

figure;
plot(time_s, normalize(bpAz), 'k-');
hold on;
plot(time_s, normalize(bpGx), 'b-');
plot(time_s, normalize(bpGy), 'r-');
plot(time_s, normalize(bpGz), 'g-');
hold off;
xlabel('Time, seconds');
ylabel('normBp.K');
title('Normalized band-passed decimated kinematics');
legend('Az','Gx','Gy','Gz');
grid;


bpOdba = filtfilt(d, KINEMATICS.odba);
bpOdav = filtfilt(d, odav);
bpNormOdav = filtfilt(d, normOdav);

odav25 = sqrt(K.gx.^2 + K.gy.^2 + K.gz.^2);
normOdav25 = normalize(odav25);
bpOdba25 = filtfilt(d, K.odba);
bpOdav25 = filtfilt(d, odav25);
bpNormOdav25 = filtfilt(d, normOdav25);

altOdav25 = K.gx + K.gy + K.gz;
bpAltOdav25 = filtfilt(d, altOdav25);
normAltOdav25 = normalize(altOdav25);
bpNormAltOdav25 = filtfilt(d, normAltOdav25);

d2AltOdav25 = diff(altOdav25, 2); d2AltOdav25(end+1:end+2) = 0;

%% Windowed-look at ODAV and 2nd derivative ODAV

sampleStart = 13550;
windowLength = 10 * K.fs;  % 2 seconds expressed in K.fs samples
sampleEnd = sampleStart + windowLength;

QQQ = normalize(altOdav25(sampleStart:sampleEnd));
d2QQQ = normalize(d2AltOdav25(sampleStart:sampleEnd));
bpQQQ = filtfilt(d, QQQ);
bpD2QQQ = filtfilt(d, d2QQQ);

figure;

figODGS1 = subplot(211);
plot(time(sampleStart:sampleEnd), bpQQQ, 'b-'); 
xlabel('Time (local)');
ylabel('Additive ODGS');
grid;

figODGS2 = subplot(212);
plot(time(sampleStart:sampleEnd), bpD2QQQ, 'b-'); 
xlabel('Time (local)');
ylabel('2° Additive ODGS');
grid;

linkaxes([figODGS1 figODGS2],'x');

%% EMD and Hilbert-Huang transforms the gyro and accel data

sampleStart = 13625;
windowLength = 10 * K.fs;  % 2 seconds expressed in K.fs samples
sampleEnd = sampleStart + windowLength;

testSample = altOdav25(sampleStart:sampleEnd);
bpTestSample = bpAltOdav25(sampleStart:sampleEnd);

log10TestSample = real(log10(testSample));
log10BpTestSample = real(log10(bpTestSample));

imfTestSample   = emd(testSample);
imfBpTestSample = emd(bpTestSample);

figure;
hht(imfTestSample, K.fs, 'FrequencyLimits',[0 5]);
title('IMF: altOdav (gx + gy + gz)');

figure;
hht(imfBpTestSample, K.fs, 'FrequencyLimits',[0 5]);
title('IMF: Band-passed altOdav (gx + gy + gz)');

figure; vmd(testSample);
figure; hht(vmd(testSample),K.fs);
figure; vmd(bpTestSample);
figure; hht(vmd(bpTestSample), K.fs, 'FrequencyLimits',[0 4]);

%% 
figure; 
figA = subplot(211);
plot(KINEMATICS.time, bpOdba, 'b-');
xlabel('Time (local)');
ylabel('ODBA');
title('Overall dynamic body accelerometer');
grid;

figB = subplot(212);
plot(KINEMATICS.time, bpOdav, 'k-');
xlabel('Time (local)');
ylabel('ODGS');
title('Overall dynamic gyroscopic signal');
grid;

linkaxes([figA figB],'x');

figure; 
figA = subplot(211);
plot(time, bpOdba25, 'b-');
xlabel('Time (local)');
ylabel('ODBA');
title('Overall dynamic body accelerometer 25 Hz');
grid;

figB = subplot(212);
plot(time, bpOdav25, 'k-');
xlabel('Time (local)');
ylabel('ODGS');
title('Overall dynamic gyroscopic signal 25 Hz');
grid;

linkaxes([figA figB],'x');

%%

sampleStart = 13550;
windowLength = 10 * K.fs;
sampleEnd = sampleStart + windowLength;

[imf_rLed2, imfResidual_rLed2] = vmd(O.rLed2,'NumIMF',9);
[imf_bpOdba, imfResidual_bpOdba] = vmd(K.odba,'NumIMF',9);

vmd(K.gx(sampleStart:sampleEnd), 'NumIMF', 9);

[imf_gx, imfResidual_gx] = vmd(K.gx, 'NumIMF', 9);
[imf_gy, imfResidual_gy] = vmd(K.gy, 'NumIMF', 9);

signalOfInterest = imf_gx;

t = tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
for n = 1:9
    ax(n) = nexttile(t);
    plot(time_s,signalOfInterest(:,n)')
    xlim([time_s(1) time_s(end)])
    txt = ['IMF',num2str(n)];
    title(txt)
    xlabel('Time (seconds)')
end
title(t,'Variational Mode Decomposition')

cleanSignal = sum(signalOfInterest(:,2:8),2);
cleanSignal = cleanSignal';
figure;
plot(time_s, K.gx, 'k-');
hold on;
plot(time_s, cleanSignal, 'b-');
hold off;
legend('gx','cleanGx','Location','southeast');
grid;

%% seismocardiogram scratchspace

inertialMoment = 10; % dirtyEnableGyro guess value

ensembleGxyz = 0.5 * ( (inertialMoment * K.gx(sampleStart:sampleEnd) .^ 2 ) + ...
    ( inertialMoment * K.gy(sampleStart:sampleEnd) .^ 2 ) + ...
    ( inertialMoment * K.gz(sampleStart:sampleEnd) .^2  ) );

ensembleGxy = 0.5 * ( (inertialMoment * K.gx(sampleStart:sampleEnd) .^ 2 ) + ...
    ( inertialMoment * K.gy(sampleStart:sampleEnd) .^ 2 ) + ...
    0); % ( inertialMoment * K.az(sampleStart:sampleEnd) .^2  ) );

ensembleGx = 0.5 * ( (inertialMoment * K.gx(sampleStart:sampleEnd) .^ 2 ) + ...
    0 + ... % ( inertialMoment * K.gy(sampleStart:sampleEnd) .^ 2 ) + ...
    0); % ( inertialMoment * K.az(sampleStart:sampleEnd) .^2  ) );

figure;

figSCGa = subplot(211);
plot(time_s(sampleStart:sampleEnd), ensembleGx, 'b-');grid;
xlabel('Time, seconds');
ylabel('SCG Gxyz');
title('SCG Test | Gx·Gy·Gz');
grid;

figSCGb = subplot(212);
plot(time_s(sampleStart:sampleEnd), ensembleGxy, 'k-');grid;
xlabel('Time, seconds');
ylabel('SCG Gxy');
title('SCG Test | Gx·Gy');
grid;
linkaxes([figSCGa figSCGb],'x');

%% create an optical vector and compute means

activeOpticsVector = [O.rLed2; O.rLed3];
ambientOpticsVector = [O.rLed1; O.rLed4];

normActiveOpticsVector = [O.nrLed2; O.nrLed3];
normAmbientOpticsVector = [O.nrLed1; O.nrLed4];

activeOpticsMean = mean(activeOpticsVector);
ambientOpticsMean = mean(ambientOpticsVector);

normActiveOpticsMean = mean(normActiveOpticsVector);
normAmbientOpticsMean = mean(normAmbientOpticsVector);


%% create an accel and gyro vector and compute the mean using all axes

accelVector = [K.ax; K.ay; K.az];
normAccelVector = [ normalize(K.ax); normalize(K.ay); normalize(K.az) ];

gyroVector  = [K.gx; K.gy; K.gz];
normGyroVector = [ normalize(K.gx); normalize(K.gy); normalize(K.gz) ];

accelMean   = mean(accelVector);
gyroMean    = mean(gyroVector);

normAccelMean = mean(normAccelVector);
normGyroMean = mean(normGyroVector);

figure;

subplot(211);
plot(time_s, accelMean, 'r-');
title('Normalized decimated accelerometer & gyroscope data');
xlabel('Time, seconds');
ylabel('Mean norm accel');
grid;

subplot(212);
plot(time_s, gyroMean, 'm-');
xlabel('Time, seconds');
ylabel('Mean norm gyro'); 
grid;


%% put it all together...

%  plot means norm reflectance optics, mean accel & gyro, and natural odba

figure;

multiSensorA = subplot(411);
plot(time_s, normActiveOpticsMean, 'b-', 'LineWidth', 2);
hold on;
plot(time_s, normAmbientOpticsMean, 'k--');
hold off;
grid;
xlabel('Time, seconds');
ylabel('Norm Reflectance');
title('Mean normalized active and ambient optical intensity');

multiSensorB = subplot(412);
plot(time_s, accelMean, 'r-');
grid;
xlabel('Time, seconds');
ylabel('m/s^2');
title('Mean accelerometer signal');

multiSensorC = subplot(413);
plot(time_s, gyroMean, 'm-');
grid;
xlabel('Time, seconds');
ylabel('°/s');
title('Mean gyroscopic signal');

multiSensorD = subplot(414);
plot(time_s, K.odba, 'k-');
grid;
xlabel('Time, seconds');
ylabel('ODBA');
title('Normalized overall dynamic body acceleration');

linkaxes([multiSensorA, multiSensorB, multiSensorC, multiSensorD], 'x');


%% Diagnostic optical plots

fprintf('\tPlotting raw LED1-4 signals...\n');

figure;

led1Plot = subplot(411);
plot(OPTICS.time_s, OPTICS.led1, 'k-');
title('LED1');
xlabel('Time, seconds');
ylabel('Intensity');
grid;

led2Plot = subplot(412);
plot(OPTICS.time_s, OPTICS.led2, 'b-');
title('LED2');
xlabel('Time, seconds');
ylabel('Intensity');
grid;

led3Plot = subplot(413);
plot(OPTICS.time_s, OPTICS.led3, 'r-');
title('LED3');
xlabel('Time, seconds');
ylabel('Intensity');
grid;

led4Plot = subplot(414);
plot(OPTICS.time_s, OPTICS.led4, 'g-');
title('LED4');
xlabel('Time, seconds');
ylabel('Intensity');
grid;

linkaxes([led1Plot led2Plot led3Plot led4Plot], 'x');

%% Diagnostic optical plots with theoretical ambient (LED1) subtraction

figure;

led1Plot = subplot(211);
plot(OPTICS.time_s, OPTICS.rLed2 - OPTICS.rLed4, 'b-');
title('LED1');
xlabel('Time, seconds');
ylabel('Reflectane Intensity');
grid;

led2Plot = subplot(212);
plot(OPTICS.time_s, OPTICS.rLed3 - OPTICS.rLed1, 'r-');
title('LED2 - LED1');
xlabel('Time, seconds');
ylabel('Intensity');
grid;

linkaxes([led1Plot led2Plot], 'x');

%% Do a bandpass filter on 250 Hz optical data and decimate a copy to 25 Hz

fprintf('\tApplying a band-pass filter with the following properties:\n');

% Bandpass example is a 500 to 560 Hz bandpass, 
% [A,B,C,D] = butter(10,[500 560]/750); 
% d = designfilt('bandpassiir','FilterOrder',20, ...
%    'HalfPowerFrequency1',500,'HalfPowerFrequency2',560, ...
%    'SampleRate',1500);
%    //  where 10 = 1/2 filter order?
%    //       750 = 1/2 sampling frequency   


% these seem good for heart beat detection
bwFiltOrder = 6;    % 6th order works well for HR peak detection
%lowCut      = 0.3333;   % 0.3333 Hz = ~20 beats per minute
lowCut      = 0.5;      % 0.5 Hz = 30 beats per minute
highCut     = 3.5;      % 3.5 Hz    = 210 beats per minute
testFs = 8;             % 8 produces big pulses around odba HR in bpLed2
% testFs = 9;             % 9 produces big pulses around odba HR in bpLed2
% testFs = 25;            % 25 produces some overlaps in bpLed3 but washes
                          % out the bpLed2 signal. Maybe this suggests
                          % using separate filters for LED2 and LED3?
% testFs = 250;           % the default value  


% supposedly good for retaining high resolution physio signals
bwFiltOrder = 2;    % supposed preserves valuable physio (valve open/close)
lowCut = 1;
highCut = 10;
testFs = 25;

% numbers suggested by the fieldtriptoolbox folks - good for O2 sat/desat
% lowCut      = 0.01;       % 0.01 Hz   = fluctuations over 10 seconds
% highCut     = 0.1;        % 0.1 Hz    = fluctuations over 100 seconds

% try some different values for nyquistFs...
% here's an explainer of some prelim results / looks at different values
                        
nyquistFs   = testFs;        

fprintf('\t\tFilter order: %d\n', bwFiltOrder);
fprintf('\t\tLow cut frequency: %d Hz | low HR: %d BPM\n', ...
    lowCut, (lowCut/(1/60)) );
fprintf('\t\tHigh cut frequency: %d Hz | high HR: %d BPM\n', ...
    highCut, (highCut/(1/60)) );


% experiment with a lower sampling frequency, try for 10 Hz to get around a 
% rapid roll-off at 5 Hz

[A, B, C, D] = butter(bwFiltOrder/2, [lowCut highCut] / (nyquistFs/2) );

d = designfilt('bandpassiir','FilterOrder',bwFiltOrder, ...
    'HalfPowerFrequency1',lowCut,'HalfPowerFrequency2',highCut, ...
	'SampleRate', nyquistFs);

% Convert the state-space representation to second-order sections. 
% Visualize the frequency responses using fvtool.

sos = ss2sos(A,B,C,D);
fvt = fvtool(sos, d, 'Fs', nyquistFs);
legend(fvt,'butter','designfilt');

%% if filter is what's wanted / needed, proceed with data filtering!

bpRLed1 = filtfilt(d, OPTICS.rLed1');
bpRLed2 = filtfilt(d, OPTICS.rLed2');
bpRLed3 = filtfilt(d, OPTICS.rLed3');
bpRLed4 = filtfilt(d, OPTICS.rLed4');

bpRLed1_25Hz = filtfilt(d, O.rLed1');
bpRLed2_25Hz = filtfilt(d, O.rLed2');
bpRLed3_25Hz = filtfilt(d, O.rLed3');
bpRLed4_25Hz = filtfilt(d, O.rLed4');

bpOdba_25Hz = filtfilt(d, K.odba');

bpGx_25Hz = filtfilt(d, K.gx');
bpGy_25Hz = filtfilt(d, K.gy');

bpAz_25Hz = filtfilt(d, K.az');

%% plot band-passed high resolution optical data

fprintf('\tPlotting band-passed optical data for LEDs 1-4\n');

figure;

bpFigA = subplot(411);
plot(OPTICS.time_s, bpRLed1, 'k-');
xlabel('Time, seconds');
ylabel('LED1 ambient');
title('Bandpass filtered 250 Hz optical data');
grid;

bpFigB = subplot(412);
plot(OPTICS.time_s, bpRLed2, 'b-');
xlabel('Time, seconds');
ylabel('LED2');
grid;

bpFigC = subplot(413);
plot(OPTICS.time_s, bpRLed3, 'r-');
xlabel('Time, seconds');
ylabel('LED3');
grid;

bpFigD = subplot(414);
plot(OPTICS.time_s, bpRLed4, 'g-');
xlabel('Time, seconds');
ylabel('LED4');
grid;

linkaxes( [bpFigA bpFigB bpFigC bpFigD], 'x');

% plot all of the bpLeds overlaid

figure;
plot(OPTICS.time_s, bpRLed1, 'k-');
hold on;
plot(OPTICS.time_s, bpRLed2, 'b-');
plot(OPTICS.time_s, bpRLed3, 'r-');
plot(OPTICS.time_s, bpRLed4, 'g-');
hold off;
xlabel('Time, seconds');
ylabel('Intensity, LEDs 1-4');
title('All Bandpassed LED Samples Overlaid');
grid;

%% plot decimated overlaid optical channels and band-passed ODBA

figure;
bpDec1 = subplot(311);
plot(time_s, bpRLed1_25Hz, 'k-');
hold on;
plot(time_s, bpRLed2_25Hz * -1, 'b-');
plot(time_s, bpRLed3_25Hz * -1, 'r-');
plot(time_s, bpRLed4_25Hz, 'g-');
hold off;
xlabel('Time, seconds');
ylabel('Band-passed optics');
title('Overlaid band-passed optical channels');
grid;

bpDec2 = subplot(312);
plot(time_s, K.normOdba * -1, 'k-');
xlabel('Time, seconds');
ylabel('Norm ODBA');
title('Normalized ODBA');
grid;

bpDec3 = subplot(313);
plot(time_s, bpOdba_25Hz * -1, 'k-');
xlabel('Time, seconds');
ylabel('Band-passed ODBA');
title('Experimental ODBA HR plot');
grid;

linkaxes([bpDec1 bpDec2, bpDec3], 'x');

%% ambient subtraction from band-passed LED2-4

bpLed2NoAmb = bpLed2 - bpLed4;
bpLed3NoAmb = bpLed3 - bpLed1;

decBpLed2NoAmb = decimate(bpLed2NoAmb, 10);
decBpLed3NoAmb = decimate(bpLed3NoAmb, 10);

%% plot ambient subtraction overlay

fprintf('\tPlotting ambient light subtraction from LEDs 2 & 3\n');

figure;

noAmb1 = subplot(211);

plot(time_s, decBpLed2NoAmb, 'b-');
hold on;
plot(time_s, decBpLed3NoAmb, 'r-');
hold off;
xlabel('Time, seconds');
ylabel('Intensity, LEDs 2-4');
title('All Bandpassed LED Samples Overlaid with Ambient Subraction');
grid;

noAmb2 = subplot(212);
plot(time_s, normOdba_25Hz, 'k-');
xlabel('Time, seconds');
ylabel('ODBA');
title('Bandpassed Overall Dynamic Body Acceleration');
grid;

linkaxes([noAmb1 noAmb2], 'x');

%% Experiment with abs and/or squaring these results

figure;

p1a = subplot(311);
plot( time_s, decBpLed2NoAmb * -1, 'b'); 
xlabel('Time, seconds');
ylabel('Intensity');
ylim('auto');
title('Band-passed ambient-filtered LED2 @ 25 Hz');
grid;

p1b = subplot(312);
plot( time_s, decBpLed3NoAmb * -1, 'r-'); 
xlabel('Time, seconds');
ylabel('Intensity');
ylim('auto');
title('Band-passed ambient-filtered LED3 @ 25 Hz');
grid;

p1c = subplot(313);
plot(time_s, normOdba_25Hz, 'k-');
xlabel('Time, seconds');
ylabel('ODBA');
ylim('auto');
title('Normalized Overall Dynamic Body Acceleration @ 25 Hz');
grid;

linkaxes([p1a p1b p1c], 'x');

%% Do or do not do derivatives. There is no try.

doDerivatives = false;

%% Generate the 1st & 2nd derivatives of LEDs 1-4

if (doDerivatives)
    
fprintf('\tGenerating 1st and 2nd derivatives of band-passed LEDs 1-4...\n');

d1led1 = diff(bpLed1) / (opticsFs * 0.1) ; d1led1(end + 1, 1) = 0;
d1led2 = diff(bpLed2) / (opticsFs * 0.1) ; d1led2(end + 1, 1) = 0;
d1led3 = diff(bpLed3) / (opticsFs * 0.1) ; d1led3(end + 1, 1) = 0;
d1led4 = diff(bpLed4) / (opticsFs * 0.1) ; d1led4(end + 1, 1) = 0;

d2led1 = diff(bpLed1)/ (opticsFs * 0.1) ; d2led1(end + 1, 1) = 0;
d2led2 = diff(bpLed2)/ (opticsFs * 0.1) ; d2led2(end + 1, 1) = 0;
d2led3 = diff(bpLed3)/ (opticsFs * 0.1) ; d2led3(end + 1, 1) = 0;
d2led4 = diff(bpLed4)/ (opticsFs * 0.1) ; d2led4(end + 1, 1) = 0;

%% plot 1st & 2nd derivatives of band-pass filtered optical data

fprintf('\tPlotting 1st and 2nd derivatives of band-passed LEDs 1-4...\n');

figure; 

fig4a = subplot(311);

plot(opticsTime_s, bpLed1, 'k-');
xlabel('Time, seconds');
ylabel('LED1 Intensity');
% ylim([-200, 200]);
tempTitle = sprintf('Filtered LED1 PPG w/ 1° and 2nd Derivatives');
title(tempTitle);
grid;

fig4b = subplot(312);
plot(opticsTime_s', d1led1, 'k-');
xlabel('Time, seconds');
ylabel('1° LED1 Intensity');
grid;

fig4c = subplot(313);
plot(opticsTime_s', d2led1,'k-');
xlabel('Time, s');
ylabel('2° Red Intensity');
grid;

linkaxes( ([fig4a, fig4b, fig4c]), 'x');

%% plot 1st and 2nd ppg derivatives for led2

figure; 

fig4a = subplot(311);

plot(opticsTime_s, bpLed2, 'b-');
xlabel('Time, seconds');
ylabel('LED2 Intensity');
tempTitle = sprintf('Filtered LED2 PPG w/ 1° and 2nd Derivatives');
title(tempTitle);
grid;

fig4b = subplot(312);
plot(opticsTime_s, d1led2, 'b-');
xlabel('Time, seconds');
ylabel('1° LED2 Intensity');
grid;

fig4c = subplot(313);
plot(opticsTime_s, d2led2,'b-');
xlabel('Time, s');
ylabel('2° Red Intensity');
grid;

linkaxes( ([fig4a, fig4b, fig4c]), 'x');

%% plot 1st and 2nd ppg derivatives for led3

figure; 

fig4a = subplot(311);

plot(opticsTime_s, bpLed3, 'r-');
xlabel('Time, seconds');
ylabel('LED3 Intensity');
tempTitle = sprintf('Filtered LED3 PPG w/ 1° and 2nd Derivatives');
title(tempTitle);
grid;

fig4b = subplot(312);
plot(opticsTime_s, d1led3, 'r-');
xlabel('Time, seconds');
ylabel('1° LED3 Intensity');
grid;

fig4c = subplot(313);
plot(opticsTime_s, d2led3,'r-');
xlabel('Time, s');
ylabel('2° LED3 Intensity');
grid;

linkaxes( ([fig4a, fig4b, fig4c]), 'x');

%% plot 1st and 2nd ppg derivatives for led4

figure; 

fig4a = subplot(311);

plot(opticsTime_s, bpLed4, 'g-');
xlabel('Time, seconds');
ylabel('LED4 Intensity');
tempTitle = sprintf('Filtered LED4 PPG w/ 1° and 2nd Derivatives');
title(tempTitle);
grid;

fig4b = subplot(312);
plot(opticsTime_s, d1led4, 'g-');
xlabel('Time, seconds');
ylabel('1° LED4 Intensity');
grid;

fig4c = subplot(313);
plot(opticsTime_s, d2led4,'g-');
xlabel('Time, s');
ylabel('2° LED4 Intensity');
grid;

linkaxes( ([fig4a, fig4b, fig4c]), 'x');

end


%% Diagnostic kinematics plots

fprintf('\tPlotting kinematic data...\n');

figure;

accelPlot = subplot(411);
plot(moveTime_s, ax, 'b-');
hold on;
plot(moveTime_s, ay, 'r-');
plot(moveTime_s, az, 'g-');
hold off;
title('Accelerometer');
xlabel('Time, seconds');
ylabel('m/s\^2');
grid;

gyroPlot = subplot(412);
plot(moveTime_s, gx, 'b-');
hold on;
plot(moveTime_s, gy, 'r-');
plot(moveTime_s, gz, 'g-');
hold off;
title('Gyroscope');
xlabel('Time, seconds');
ylabel('DPS');
grid;

magPlot = subplot(413);
plot(magTime_s, mx, 'b-');
hold on;
plot(magTime_s, my, 'r-');
plot(magTime_s, mz, 'g-');
hold off
title('Magnetometer');
xlabel('Time, seconds');
ylabel('Degrees');
grid;

prhPlot = subplot(414);
plot(moveTime_s, pitch, 'b-');
hold on;
plot(moveTime_s, roll, 'r-');
plot(moveTime_s, heading, 'g-');
hold off
title('Pitch Roll and Heading');
xlabel('Time, seconds');
ylabel('Degrees');
grid;

linkaxes([accelPlot gyroPlot magPlot prhPlot], 'x');

%% Diagnostic depth and temperature plots

fprintf('\tPlotting depth and temperature data...\n');

figure; 

depthPlot = subplot(211);
plot(moveTime_s, depth, 'b-');
set(gca, 'YDir','reverse');
title('Depth');
xlabel('Time, seconds');
ylabel('Meters');
grid;

tempPlot = subplot(212);
plot(moveTime_s, temperature, 'r-');
title('Temperature');
xlabel('Time, seconds');
ylabel('Deg Celsius');
grid;

linkaxes([depthPlot tempPlot], 'x');

%% Diagnostic power measurement plots

fprintf('\tPlotting FaunaTag power usage...\n');

figure;

voltagePlot = subplot(511);
plot(powerTime_s, voltage, 'b-');
title('Voltage');
xlabel('Time, seconds');
ylabel('mV');
grid;

currentPlot = subplot(512);
plot(powerTime_s, current, 'g-');
title('Current');
xlabel('Time, seconds');
ylabel('mA');
grid;

powerPlot = subplot(513);
plot(powerTime_s, power, 'r-');
title('Power Consumption');
xlabel('Time, seconds');
ylabel('mW');
grid;

chargePlot = subplot(514);
plot(powerTime_s, chargeState, 'm-');
title('State of Charge');
xlabel('Time, seconds');
ylabel('%');
grid;

capacityPlot = subplot(515);
plot(powerTime_s, remainingCapacity, 'k-');
title('Remaining Capacity');
xlabel('Time, seconds');
ylabel('mAh');
grid;

linkaxes([voltagePlot currentPlot powerPlot chargePlot capacityPlot], 'x');

%% AC/DC and FFT with decimated optical signals

% use decimated LED1 as signal length
sigLength = length(led1_25Hz);

% deprecate Ts, since we always know Fs
% Ts     = mean(diff(time_s));          % average time sample

% declare Fs explicitly, since we always decimate to 25 Hz here
% Fs     = 1 / Ts;                      % effective sampling rate
% Fs = 25;                              <-- previously declared

% now try some FFT magic...
nFFT   = 2 ^ nextpow2(sigLength);     % next power of 2 from sig length
ambFFT = fft(led1_25Hz, nFFT) / sigLength;      % LED1 FFT
led2FFT = fft(led2_25Hz, nFFT) / sigLength;     % LED2 FFT
led3FFT = fft(led3_25Hz, nFFT) / sigLength;     % LED3 FFT
led4FFT = fft(led4_25Hz, nFFT) / sigLength;     % LED4 FFT
odbaFFT = fft(odba_25Hz, nFFT) / sigLength;     % ODBA FFT

decBpLed2NoAmbFFT = fft(decBpLed2NoAmb', nFFT) / sigLength;
decBpLed3NoAmbFFT = fft(decBpLed3NoAmb', nFFT) / sigLength;

f      = Fs/2 * linspace(0,1,nFFT/2+1);
freq   = -Fs/2 + Fs/nFFT:Fs/nFFT:Fs/2;


ambCenter  = fftshift(ambFFT);
led2Center = fftshift(led2FFT);
led3Center = fftshift(led3FFT);
led4Center = fftshift(led4FFT);
odbaCenter = fftshift(odbaFFT);

decBpLed2Center = fftshift(decBpLed2NoAmbFFT);
decBpLed3Center = fftshift(decBpLed3NoAmbFFT);

ambAC      = fftshift(ambFFT(2:end));
led2AC     = fftshift(led2FFT(2:end));
led3AC     = fftshift(led3FFT(2:end));
led4AC     = fftshift(led4FFT(2:end));
odbaAC     = fftshift(odbaFFT(2:end));

decBpLed2NoAmbAC = fftshift(decBpLed2NoAmbFFT(2:end));
decBpLed3NoAmbAC = fftshift(decBpLed3NoAmbFFT(2:end));

ambDC      = ambFFT(1);
led2DC     = led2FFT(1);
led3DC     = led3FFT(1);
led4DC     = led4FFT(1);
odbaDC     = odbaFFT(1);

%% try to plot some of this FFT magic...

figure;

figFFT1a = subplot(411);
plot(freq(2:end), abs(ambAC),'k-');
xlim([ 0.166666667 4.3 ]);
ylabel('abs(amb(f))');
title('Ambient AC Single-Sided Amplitude Spectrum');

figFFT1b = subplot(412);
plot(freq(2:end), abs(led2AC),'b-');
xlim([ 0.166666667 4.3 ]);
ylabel('abs(led2(f))');
title('LED2 AC Single-Sided Amplitude Spectrum');

figFFT1c = subplot(413);
plot(freq(2:end), abs(led3AC),'r-');
xlim([ 0.166666667 4.3 ]);
ylabel('abs(LED3(f))');
title('LED3 AC Single-Sided Amplitude Spectrum');

figFFT1d = subplot(414);
plot(freq(2:end), abs(led3AC),'g-');
xlim([ 0.166666667 4.3 ]);
ylabel('abs(LED4(f))');
xlabel('Frequency (Hz)');
title('LED4 AC Single-Sided Amplitude Spectrum');

linkaxes( [figFFT1a, figFFT1b, figFFT1c, figFFT1d], 'x');

figure;
plot(freq(2:end), abs(odbaAC),'k-');
xlim([ 0.166666667 4.3 ]);
ylabel('abs(ODBA(f))');
xlabel('Frequency (Hz)');
title('ODBA AC Single-Sided Amplitude Spectrum');

% decBpLed2 & 3 noAmb FFTs

figure;

fft1a = subplot(211);
plot(freq(2:end), abs(decBpLed2NoAmbAC),'b-');
xlim([ 0.166666667 4.3 ]);
ylabel('abs(decBpLed2NoAmb(f))');
xlabel('Frequency (Hz)');
title('Decimated Band-passed LED2 AC Single-Sided Amplitude Spectrum');

fft1b = subplot(212);
plot(freq(2:end), abs(decBpLed3NoAmbAC),'r-');
xlim([ 0.166666667 4.3 ]);
ylabel('abs(decBpLed3NoAmb(f))');
xlabel('Frequency (Hz)');
title('Decimated Band-passed LED3 AC Single-Sided Amplitude Spectrum');


%% apply an x-th order butterworth filter & appropriate biological cutoff

sampleInterval    = Fs;                     
samplingRate      = 1 / Fs;
% At ~50 Hz, cutoffValues tranlate to cutoffFreq and BPM, as follows:
% cutoffValue =  100 -> cutOffFreq = 0.5013 Hz -> BPM = 30.078
% cutoffValue =  125 -> cutOffFreq = 0.4010 Hz -> BPM = 24.063
% cutoffValue =  150 -> cutOffFreq = 0.3342 Hz -> BPM = 20.052
% cutoffValue =  200 -> cutOffFreq = 0.2507 Hz -> BPM = 15.039
% cutoffValue =  250 -> cutOffFreq = 0.2005 Hz -> BPM = 12.030
% cutoffValue =  333 -> cutOffFreq = 0.1505 Hz -> BPM = 9.0300
% cutoffValue =  500 -> cutOffFreq = 0.1003 Hz -> BPM = 6.0180
% cutoffValue = 1000 -> cutOffFreq = 0.0501 Hz -> BPM = 3.0060
% cutoffValue = 2000 -> cutOffFreq = 0.0251 Hz -> BPM = 1.5060
% cutoffValue = 3000 -> cutOffFreq = 0.0167 Hz -> BPM = 1.0020                             
% cutoffValue = 4000 -> cutOffFreq = 0.0125 Hz -> BPM = 0.7500                               
cutoffValue       = 200;  % was 4000, but trying with 3000 = 1 BPM          
cutoffFrequency   = 1 / (samplingRate * cutoffValue) ;    
butterOrder       = 6;
butterType        = 'high';
% butterNyquistNormalizedCutoff = cutoffFrequency / (samplingRate / 2);
butterNyquistNormalizedCutoff = cutoffFrequency;

[b, a] = butter(butterOrder, butterNyquistNormalizedCutoff, butterType);

filtered_led1_25Hz = filtfilt( b, a, led1_25Hz );
filtered_led2_25Hz = filtfilt( b, a, led2_25Hz );
filtered_led3_25Hz = filtfilt( b, a, led3_25Hz );
filtered_led4_25Hz = filtfilt( b, a, led4_25Hz );


%% Plot individual filtered ppgWaveforms

figure;

filtPPG1 = subplot(411);
plot(time_s, filtered_led1_25Hz, 'k-');
xlabel('Time, s');
ylabel('Filtered intensity LED1');
tempTitle = sprintf('Filtered LED1-4 (%d° BW w %0.3g Hz Cutoff)', ...
    butterOrder, cutoffFrequency);
title(tempTitle);
grid;

filtPPG2 = subplot(412);
plot(time_s, filtered_led2_25Hz, 'b-');
xlabel('Time, s');
ylabel('Filtered intensity LED2');
grid;

filtPPG3 = subplot(413);
plot(time_s, filtered_led3_25Hz, 'r-');
xlabel('Time, s');
ylabel('Filtered intensity LED3');
grid;

filtPPG4 = subplot(414);
plot(time_s, filtered_led4_25Hz, 'g-');
xlabel('Time, s');
ylabel('Filtered intensity LED4');
grid;

linkaxes([ filtPPG1 filtPPG2 filtPPG3 filtPPG4], 'x');

%% Plot overlaid filtered ppgWaveforms 

figure;

overlay1a = subplot(211);
plot(time_s, filtered_led1_25Hz * -1,'k-');
xlabel('Time, s');
ylabel('LED1-4 Reflectance');
tempTitle = sprintf('Filtered LED1-4 (%d° BW w %0.3g Hz Cutoff)', ...
    butterOrder, cutoffFrequency);
title(tempTitle);
hold on;
plot(time_s, filtered_led2_25Hz * -1, 'b-');
plot(time_s, filtered_led3_25Hz * -1, 'r-');
plot(time_s, filtered_led4_25Hz * -1, 'g-');
hold off;
grid;

overlay1b = subplot(212);
plot(time_s, normOdba_25Hz, 'k-');
xlabel('Time, s');
ylabel('ODBA');
title('Normalized ODBA');
grid;

linkaxes([overlay1a overlay1b], 'x');

%% Plot overlaid filtered ppgWaveforms - ambient w/ normOdba

figure;

overlay1a = subplot(211);
plot(time_s, (filtered_led1_25Hz - filtered_led4_25Hz) * -1,'b-');
xlabel('Time, s');
ylabel('LED2 & LED3');
tempTitle = sprintf('Filtered LED2 & LED3 (%d° BW w %0.3g Hz Cutoff)', ...
    butterOrder, cutoffFrequency);
title(tempTitle);
hold on;
plot(time_s, (filtered_led3_25Hz - filtered_led1_25Hz) * -1, 'r-');
hold off;
grid;

overlay1b = subplot(212);
plot(time_s, normOdba_25Hz, 'k-');
xlabel('Time, s');
ylabel('ODBA');
title('Normalized ODBA');
grid;

linkaxes([overlay1a overlay1b], 'x');

%% Plot the filtered ppgWaveforms in individual panels

figure;

fig3a = subplot(411);
plot(time_s, filtered_led1_25Hz, 'k-');
xlabel('Time, seconds');
ylabel('Ambient Reflectance');
% ylim([-200 200]);
tempTitle = sprintf('Filtered ambient ((%d° BW w %0.3g Hz Cutoff)', ...
    butterOrder, cutoffFrequency);
title(tempTitle);
grid;

fig3b = subplot(412);
plot(time_s, filtered_led2_25Hz,'b-');
xlabel('Time, seconds');
ylabel('LED2 1050 nm Reflectance');
% ylim([-200 200]);
grid;

fig3c = subplot(413);
plot(time_s, filtered_led3_25Hz, 'r-');
xlabel('Time, s');
ylabel('LED3 1200·1300 nm Reflectance');
% ylim([-200 200]);
grid;

fig3d = subplot(414);
plot(time_s, filtered_led4_25Hz, 'g-');
xlabel('Time, s');
ylabel('LED4 980 m, Reflectance');
% ylim([-2000 2000]);
grid;

linkaxes( [fig3a, fig3b, fig3c, fig3d], 'x');

%% Plot the filtered ppgWaveforms minus ambient LED1 in individual panels

figure;

fig3a = subplot(411);
plot(time_s, filtered_led2_25Hz - filtered_led1_25Hz, 'b-');
xlabel('Time, seconds');
ylabel('Filtered LED2-LED1 (1050nm)');
xlim([0 120]);
ylim([-5000 5000]);
tempTitle = sprintf('Filtered LED2-4 - LED1 ((%d° BW w %0.3g Hz Cutoff)', ...
    butterOrder, cutoffFrequency);
title(tempTitle);
grid;

fig3b = subplot(412);
plot(time_s, filtered_led3_25Hz - filtered_led1_25Hz ,'r-');
xlabel('Time, seconds');
ylabel('Filt LED3-LED1 (12/13)');
ylim([-5000 5000]);
grid;

fig3c = subplot(413);
plot(time_s, filtered_led4_25Hz - filtered_led1_25Hz, 'g-');
xlabel('Time, s');
ylabel('Filt LED4-LED1 (980nm)');
ylim([-5000 5000]);
grid;


fig3d = subplot(414);
plot(time_s, odba_25Hz, 'k-');
xlabel('Time, s');
ylabel('ODBA');
ylim([-10 10]);
grid;

linkaxes( [fig3a, fig3b, fig3c, fig3d], 'x');


%% Make some spectral analysis plots

welchNFFT     = 2^7;        % try segment lengths of 2^8 = 256
winSize  = hanning(welchNFFT);   % set hanning window shape
nOverlap = welchNFFT / 2;        % set 50% overlap between segments

[Pled2, Fled2] = periodogram(led2_25Hz, [], welchNFFT, Fs, 'power');
figure; 

psd1 = subplot(211);
plot(Fled2,10*log10(Pled2),'b-');
xlabel('Frequency, Hz'); ylabel('Power spectrum (dBW) LED2 1050 nm');
title('Periodogram for LED2 Samples');

[Pled2Power,Fled2Power] = pwelch(led2_25Hz,ones(welchNFFT,1),0,welchNFFT,Fs,'power');
ps2 = subplot(212);
plot(Fled2Power,10*log10(Pled2Power),'b-');
xlabel('Frequency, Hz'); ylabel('Power spectrum (dBW)');
title('Welch PSD for LED2 Samples');


[psd_led1, f_led1] = pwelch(filtered_led1_25Hz, winSize, nOverlap, welchNFFT);
% [psd_led2, f_led2] = pwelch(filtered_led2_25Hz, winSize, nOverlap, welchNFFT);
% [psd_led3, f_led3] = pwelch(filtered_led3_25Hz, winSize, nOverlap, welchNFFT);
[psd_led2, f_led2] = pwelch(decBpLed2NoAmb', winSize, nOverlap, welchNFFT);
[psd_led3, f_led3] = pwelch(decBpLed3NoAmb', winSize, nOverlap, welchNFFT);
[psd_led4, f_led4] = pwelch(filtered_led4_25Hz, winSize, nOverlap, welchNFFT);


%% Try to make some spectrograms and periodograms of pulse ox data

welchNFFT     = 2^6;        % try segment lengths of 2^8 = 256
winSize  = hanning(welchNFFT);   % set hanning window shape
nOverlap = welchNFFT / 2;        % set 50% overlap between segments

[S, F, T, P] = spectrogram(decBpLed2NoAmb(450*25:470*25,1), winSize, ... 
    nOverlap, nFFT, Fs);

% make surface3d plot
figure;
colormap(parula(5));
% surface with edgecolor
% surfc(T,F,10*log10(abs(P)));
% surface with no edgecolor
surfc(T,F,10*log10(abs(P)), 'EdgeColor','none');
xlabel('Time,s'); ylabel('Frequency, Hz'); zlabel('Magnitude, dB');
title('Spectrogram of LED2 Samples');
axis tight;
%view(-45,45);

%% make contour3d plot

figure;
colormap(parula(5));
% surfc(T,F,10*log10(abs(P)), 'EdgeColor','none');
contour3(T,F,10*log10(abs(P)));
xlabel('Time,s'); ylabel('Frequency, Hz'); 
title('3D Contour plot of LED2 sample spectra');
c = colorbar;
c.Label.String = 'Intensity, dB';
axis tight;

%% Try a static spectrogram plot

figure;
colormap(parula(5));
colorLimit = [40 90] ;  % color axis limits in dB for specgram
%figure; 
imagesc(T, F, 10*log10(abs(P)));
% contourf(10*log10(abs(P)));
xlabel('Time,seconds'); ylabel('Frequency, Hz'); 
title('Static spectrogram of LED2 samples');
c = colorbar;
c.Label.String = 'Intensity, dB';
axis tight;
axis xy;

%% Make some PSDs and plots of PSDs

doPSDPlots = 1;

if (doPSDPlots == 1)
        
    figure;

    figPSD1 = subplot(411);
    loglog(f_led1, psd_led1, '-', 'Color', 'Black');
    xlabel ('Frequency (Hz)');
    ylabel ('PSD');
    tempTitle = sprintf('Filtered Ambient PPG PSD');
    title(tempTitle);
    grid;

    figPSD2 = subplot(412);
    loglog(f_led2, psd_led2, 'b-');
    xlabel ('Frequency (Hz)');
    ylabel ('PSD');
    % draw the cutoff line on filter plot
    hold on;
    loglog( [-log10(cutoffFrequency) -log10(cutoffFrequency)], [10^2 10^8], 'k:');
    hold off;
    tempTitle = sprintf('Filtered LED2 PPG PSD');
    title(tempTitle);
    grid;

    figPSD3 = subplot(413);
    loglog(f_led3, psd_led3, '-', 'Color', 'red');
    xlabel ('Frequency (Hz)');
    ylabel ('PSD');
    hold on;
    loglog( [-log10(cutoffFrequency) -log10(cutoffFrequency)], [10^2 10^8], 'k:');
    hold off;
    tempTitle = sprintf('Filtered LED3 PPG PSD');
    title(tempTitle);
    grid;
    
    figPSD4 = subplot(414);
    loglog(f_led4, psd_led4, '-', 'Color', 'green');
    xlabel ('Frequency (Hz)');
    ylabel ('PSD');
    hold on;
    loglog( [-log10(cutoffFrequency) -log10(cutoffFrequency)], [10^2 10^8], 'k:');
    hold off;    
    tempTitle = sprintf('Filtered LED4 PPG PSD');
    title(tempTitle);
    grid;
    
    linkaxes([ figPSD1, figPSD2, figPSD3, figPSD4], 'x');
    
end

%% Get 1st and 2nd derivatives of the filtered PPG waveforms

d1led1 = diff(filtered_led1_25Hz) / (samplingRate * 0.1) ; d1led1(1, end+1) = 0;
d1led2 = diff(filtered_led2_25Hz) / (samplingRate * 0.1) ; d1led2(1, end+1) = 0;
d1led3 = diff(filtered_led3_25Hz) / (samplingRate * 0.1) ; d1led3(1, end+1) = 0;
d1led4 = diff(filtered_led4_25Hz) / (samplingRate * 0.1) ; d1led4(1, end+1) = 0;

d2led1 = diff(d1led1)/ (samplingRate * 0.1) ; d2led1(1, end+1) = 0;
d2led2 = diff(d1led2)/ (samplingRate * 0.1) ; d2led2(1, end+1) = 0;
d2led3 = diff(d1led3)/ (samplingRate * 0.1) ; d2led3(1, end+1) = 0;
d2led4 = diff(d1led4)/ (samplingRate * 0.1) ; d2led4(1, end+1) = 0;

%% plot 1st and 2nd ppg derivatives for led1

figure; 

fig4a = subplot(311);

plot(time_s, filtered_led1_25Hz, 'k-');
% [maxtab, mintab] = peakdet(filteredPPGred, 0.5);
%hold on; plot(mintab(:,1), mintab(:,2), 'g*');
%plot(maxtab(:,1), maxtab(:,2), 'r*');
xlabel('Time, seconds');
ylabel('LED1 Intensity');
% ylim([-200, 200]);
tempTitle = sprintf('Filtered LED1 PPG w/ 1° and 2nd Derivatives');
title(tempTitle);
grid;

fig4b = subplot(312);
plot(time_s, d1led1, 'k-');
xlabel('Time, seconds');
ylabel('1° LED1 Intensity');
grid;

fig4c = subplot(313);
plot(time_s, d2led1,'k-');
xlabel('Time, s');
ylabel('2° Red Intensity');
grid;

linkaxes( ([fig4a, fig4b, fig4c]), 'x');

%% plot 1st and 2nd ppg derivatives for led2

figure; 

fig4a = subplot(311);

plot(time_s, filtered_led2_25Hz, 'b-');
% [maxtab, mintab] = peakdet(filteredPPGred, 0.5);
%hold on; plot(mintab(:,1), mintab(:,2), 'g*');
%plot(maxtab(:,1), maxtab(:,2), 'r*');
xlabel('Time, seconds');
ylabel('LED2 Intensity');
% ylim([-200, 200]);
tempTitle = sprintf('Filtered LED2 PPG w/ 1° and 2nd Derivatives');
title(tempTitle);
grid;

fig4b = subplot(312);
plot(time_s, d1led2, 'b-');
xlabel('Time, seconds');
ylabel('1° LED2 Intensity');
grid;

fig4c = subplot(313);
plot(time_s, d2led2,'b-');
xlabel('Time, s');
ylabel('2° Red Intensity');
grid;

linkaxes( ([fig4a, fig4b, fig4c]), 'x');

%% plot 1st and 2nd ppg derivatives for led3

figure; 

fig4a = subplot(311);

plot(time_s, filtered_led3_25Hz, 'r-');
% [maxtab, mintab] = peakdet(filteredPPGred, 0.5);
%hold on; plot(mintab(:,1), mintab(:,2), 'g*');
%plot(maxtab(:,1), maxtab(:,2), 'r*');
xlabel('Time, seconds');
ylabel('LED3 Intensity');
% ylim([-200, 200]);
tempTitle = sprintf('Filtered LED3 PPG w/ 1° and 2nd Derivatives');
title(tempTitle);
grid;

fig4b = subplot(312);
plot(time_s, d1led3, 'r-');
xlabel('Time, seconds');
ylabel('1° LED3 Intensity');
grid;

fig4c = subplot(313);
plot(time_s, d2led3,'r-');
xlabel('Time, s');
ylabel('2° LED3 Intensity');
grid;

linkaxes( ([fig4a, fig4b, fig4c]), 'x');

%% plot 1st and 2nd ppg derivatives for led4

figure; 

fig4a = subplot(311);

plot(time_s, filtered_led4_25Hz, 'g-');
% [maxtab, mintab] = peakdet(filteredPPGred, 0.5);
%hold on; plot(mintab(:,1), mintab(:,2), 'g*');
%plot(maxtab(:,1), maxtab(:,2), 'r*');
xlabel('Time, seconds');
ylabel('LED4 Intensity');
% ylim([-200, 200]);
tempTitle = sprintf('Filtered LED4 PPG w/ 1° and 2nd Derivatives');
title(tempTitle);
grid;

fig4b = subplot(312);
plot(time_s, d1led4, 'g-');
xlabel('Time, seconds');
ylabel('1° LED4 Intensity');
grid;

fig4c = subplot(313);
plot(time_s, d2led4,'g-');
xlabel('Time, s');
ylabel('2° LED4 Intensity');
grid;

linkaxes( ([fig4a, fig4b, fig4c]), 'x');

%% do some experimental frequency domain analysis using LED2

nFFT = length(filtered_led2_25Hz);
Y    = fft(filtered_led2_25Hz, nFFT); 
F    = ((0:1/nFFT:1-1/nFFT)*Fs);

magnitudeLed2FFT = abs(Y);        % Magnitude of the FFT
phaseLed2FFT     = unwrap(angle(Y));  % Phase of the FFT

figure;
fig9a = subplot(211);
plot(F, magnitudeLed2FFT,'b-');
xlabel('Frequency, Hz');
ylabel('Magnitude, dB');
tempTitle = sprintf('LED2 FFT');
title(tempTitle);
grid;

fig9b = subplot(212);
plot(F, phaseLed2FFT, 'b-');
xlabel('Frequency, Hz');
ylabel('Phase, radians');
tempTitle = sprintf('LED2 Phase FFT');
title(tempTitle);
grid;

%% do sliding / moving window FFt/wavelet analysis for HR estimation

doSlidingWindowWaveletAnalysis = false;

if (doSlidingWindowWaveletAnalysis)
    fprintf('\tBeginning sliding window wavelet HR estimation...\n');
    
    % subsample some 10 second windows of bpLed2

    bpOpticsTime_s = opticsTime_s';

    lenLedSamples = length(bpLed2);

    if ( length(bpOpticsTime_s) == lenLedSamples )
        fprintf('\tVector time and sample sizes match. Continuing!\n');
    else
        fprintf('\tVector time and sample sizes mismatch. Stopping!\n');
        return;
    end

    windowLength = opticsFs * 10;   % i.e.: ten second window

    for i=1:lenLedSamples-windowLength

        data=bpLed2(i:i+windowLength); % Select data range
        thisMean = mean(data);

    %     nFFT = windowLength;
    %     Y    = fft(data, nFFT); 
    %     F    = ((0:1/nFFT:1-1/nFFT)*opticsFs);
    % 
    %     magnitudeLed2FFT = abs(Y);        % Magnitude of the FFT
    %     phaseLed2FFT     = unwrap(angle(Y));  % Phase of the FFT

        [wsstLed2, wsstLed2F]   = wsst(bpLed2, opticsFs);

        % build wsst ridges for signals of interest, using 'NumRidges'
        % set to 2, to get primary signal and 1st harmonic
        % and run a loop (for now) to explore the effect of 'penalty'

        penalty = 10;

        wsstLed2Ridge = wsstridge(wsstLed2, penalty, wsstLed2F, ...
            'NumRidges',2);

        led2PulseF       = wsstLed2Ridge(:,1);          
        led2Pulse1stH    = wsstLed2Ridge(:,2);
        led2PulseBPM     = 60 * led2PulseF;
        led2Pulse1stHBPM = 60 * led2Pulse1stH;

        % fprintf('\taverage bpLed2 HR estimate in window #%d = %d\n', ...
        %    i, led2PulseF );

    %     figFFT = figure;
    %     plot(F,magnitudeLed2FFT,'b-');
    %     xlim([0 4]);
    %     xlabel('Frequency, Hz');
    %     ylabel('Magnitude');
    %     title('FFT');
    %     input('Return to see next window\n');
    %     close(figFFT);


    end

    
    % nFFT = length(bpLed2);
    % Y    = fft(bpLed2, nFFT); 
    % F    = ((0:1/nFFT:1-1/nFFT)*opticsFs);
    % 
    % magnitudeLed2FFT = abs(Y);        % Magnitude of the FFT
    % phaseLed2FFT     = unwrap(angle(Y));  % Phase of the FFT
    % 
    % figure;
    % fig9a = subplot(211);
    % plot(F, magnitudeLed2FFT,'b-');
    % xlabel('Frequency, Hz');
    % ylabel('Magnitude, dB');
    % tempTitle = sprintf('LED2 FFT');
    % title(tempTitle);
    % xlim([0 4]);
    % grid;
    % 
    % fig9b = subplot(212);
    % plot(F, phaseLed2FFT, 'b-');
    % xlabel('Frequency, Hz');
    % ylabel('Phase, radians');
    % tempTitle = sprintf('LED2 Phase FFT');
    % title(tempTitle);
    % xlim([0 4]);
    % grid;    
    
else
    fprintf('\tSkipping sliding window wavelet HR estimation...\n');
end


%% do entire time-series wavelet HR estimatation (only use w/ high quality data)
   
entireTimeSeriesWaveletAnalysis = true;

% specify decimated or non-decimated data for use in wavelet analysis

useDecimated = true;

if(entireTimeSeriesWaveletAnalysis)

    doCWT = false;

    if (doCWT)

        fprintf('\tWorking on cwt...\n');

    %     % do cwt for LED1-4
    %     
    %     figure;
    %     cwt(led1_25Hz,Fs);
    %     caxis([0 400]);
    %     title('cwt(LED2)');
    %     grid;
    %     
    %     figure;
    %     cwt(led2_25Hz,Fs);
    %     caxis([0 400]);
    %     title('cwt(LED2)');
    %     grid;
    %     
    %     figure;
    %     caxis([0 400]);
    %     cwt(led3_25Hz,Fs);
    %     title('cwt(LED3)');
    %     grid;
    %     
    %     figure;
    %     caxis([0 400]);
    %     cwt(led3_25Hz,Fs);
    %     title('cwt(LED4)');
    %     grid;

        % do cwt for filtered LED1-4

        figure;
        if (useDecimated)
            cwt(filtered_led1_25Hz, Fs);
        else
            cwt(bpLed1, opticsFs);
        end
        caxis([0 200]);
        title('cwt(filteredLED1)');

        figure;
        if (useDecimated)
            cwt(filtered_led2_25Hz, Fs);
        else
            cwt(bpLed2, opticsFs);
        end
        caxis([0 200]);
        title('cwt(filteredLED2)');
        

        if (useDecimated)
            cwt(filtered_led3_25Hz, Fs);
        else
            cwt(bpLed3, opticsFs);
        end
        caxis([0 200]);
        title('cwt(filteredLED3)');
        
        
        if (useDecimated)
            cwt(filtered_led4_25Hz, Fs);
        else
            cwt(bpLed4, opticsFs);
        end
        caxis([0 200]);
        title('cwt(filteredLED4)');

    end

    % do wsst for for raw led2 and led3

    doWSST = true;

    if (doWSST)

        fprintf('\tWorking on wsst...\n');

        % do wsst for for raw LED 1-4

        figure;
        subplot(221);
        if (useDecimated)
            wsst(led1_25Hz,Fs);
        else
            wsst(led1, opticsFs);            
        end
        ylim([0 5]);
        caxis([0 10]);
        title('wsst(led1)');
        
        subplot(222);
        if (useDecimated)
            wsst(led2_25Hz,Fs);
        else
            wsst(led2, opticsFs);            
        end
        ylim([0 5]);
        caxis([0 10]);
        title('wsst(led2)');
        
        subplot(223);
        if (useDecimated)
            wsst(led3_25Hz,Fs);
        else
            wsst(led3, opticsFs);            
        end
        ylim([0 5]);
        caxis([0 10]);
        title('wsst(led3)');
        
        subplot(224);
        if (useDecimated)
            wsst(led4_25Hz,Fs);
        else
            wsst(led4, opticsFs);            
        end 
        ylim([0 5]);
        caxis([0 10]);
        title('wsst(led4)');

        % do wsst for for filtered LED 1-4

        figure;
        figWSSTled1 = subplot(221);
        if (useDecimated)
            wsst(filtered_led1_25Hz, Fs);
        else
            wsst(bpLed1, opticsFs);
        end
        ylim([0 5]);
        caxis([0 5]);
        title('wsst(filteredLed1)');

        figWSSTled2 = subplot(222);
        if (useDecimated)
            wsst(filtered_led2_25Hz, Fs);
        else
            wsst(bpLed2, opticsFs);
        end
        ylim([0 5]);
        caxis([0 10]);
        title('wsst(filteredLed2)');


        figWSSTled3 = subplot(223);
        if (useDecimated)
            wsst(filtered_led3_25Hz, Fs);
        else
            wsst(bpLed3, opticsFs);
        end
        ylim([0 5]);
        caxis([0 5]);
        title('wsst(filteredLed3)');

        figWSSTled4 = subplot(224);
        if (useDecimated)
            wsst(filtered_led4_25Hz, Fs);
        else
            wsst(bpLed4, opticsFs);
        end
        ylim([0 5]);
        caxis([0 10]);
        title('wsst(filteredLed4)');

    end
    
    doManualWSST = true;
    
    if (doManualWSST)
        
        % do some manual wsst stuff for plotting contours of least
        % penalized WSST results

        % build wsst for signals of interest

        if (useDecimated)
            
            % do both raw and filtered (band-passed) data
            [wsstLed1, wsstLed1F]         = wsst(led1_25Hz,Fs);
            [wsstLed2, wsstLed2F]         = wsst(led2_25Hz,Fs);
            [wsstLed3, wsstLed3F]         = wsst(led3_25Hz,Fs);
            [wsstLed4, wsstLed4F]         = wsst(led4_25Hz,Fs);    
            [wsstFiltLed1, wsstFiltLed1F] = wsst(filtered_led1_25Hz,Fs);
            [wsstFiltLed2, wsstFiltLed2F] = wsst(filtered_led2_25Hz,Fs);
            [wsstFiltLed3, wsstFiltLed3F] = wsst(filtered_led3_25Hz,Fs);
            [wsstFiltLed4, wsstFiltLed4F] = wsst(filtered_led4_25Hz,Fs);
            
        else
            
            [wsstLed1, wsstLed1F]         = wsst(led1,opticsFs);
            [wsstLed2, wsstLed2F]         = wsst(led2,opticsFs);
            [wsstLed3, wsstLed3F]         = wsst(led3,opticsFs);
            [wsstLed4, wsstLed4F]         = wsst(led4,opticsFs); 
            [wsstFiltLed1, wsstFiltLed1F] = wsst(bpLed1,opticsFs);
            [wsstFiltLed2, wsstFiltLed2F] = wsst(bpLed2,opticsFs);
            [wsstFiltLed3, wsstFiltLed3F] = wsst(bpLed3,opticsFs);
            [wsstFiltLed4, wsstFiltLed4F] = wsst(bpLed4,opticsFs);            
            
        end

        % build wsst ridges for signals of interest, using 'NumRidges'
        % set to 2, to get primary signal and 1st harmonic
        % and run a loop (for now) to explore the effect of 'penalty'

        doMultiPenaltyTest = 0;

        if doMultiPenaltyTest == 0

            penalty = 10;

            wsstLed2Ridge = wsstridge(wsstLed2, penalty, wsstLed2F, ...
                'NumRidges',2);
            wsstFiltLed2Ridge = wsstridge(wsstFiltLed2, penalty, ...
                wsstFiltLed2F, 'NumRidges',2);
            wsstLed3Ridge = wsstridge(wsstLed3 , penalty, wsstLed3F, ...
                'NumRidges',2);
            wsstFiltLed3Ridge = wsstridge(wsstFiltLed3, penalty, ...
                wsstFiltLed3F, 'NumRidges',2);

            led2PulseF = wsstFiltLed2Ridge(:,1);
            led3PulseF = wsstFiltLed3Ridge(:,1);
            led2Pulse1stH = wsstFiltLed2Ridge(:,2);
            led3Pulse1stH = wsstFiltLed3Ridge(:,2);

            led2PulseBPM   = 60 * led2PulseF;
            led3PulseBPM   = 60 * led3PulseF;
            led2PulseBPMh1 = 60 * led2Pulse1stH;
            led3PulseBPMh1 = 60 * led3Pulse1stH;

            % do contour plots with ridge overlays

            if (useDecimated)
               waveTime_s = time_s; 
            else
                waveTime_s = opticsTime_s;
            end
            
            figure;
            contour(waveTime_s, wsstLed2F, abs(wsstLed2));
            hold on;
            plot(waveTime_s, wsstLed2Ridge, 'b.');
            hold off;
            xlabel('Time, s'); ylabel('Frequency, Hz');
            grid on;
            ylim([0 4]);
            caxis([0 5]);
            title('Synchrosqueezed Transform - Raw LED2 PPG');

            figure;
            contour(waveTime_s,wsstFiltLed2F,abs(wsstFiltLed2));
            hold on;
            plot(waveTime_s, wsstFiltLed2Ridge, 'b.');
            hold off;
            xlabel('Time, s'); ylabel('Frequency, Hz');
            grid on;
            ylim([0 4]);
            caxis([0 5]);
            title('Synchrosqueezed Transform for Filtered LED2 PPG');            

            figure;
            contour(waveTime_s,wsstLed3F,abs(wsstLed3));
            hold on;
            plot(waveTime_s, wsstLed3Ridge, 'r.');
            hold off;
            xlabel('Time, s'); ylabel('Frequency, Hz');
            grid on;
            ylim([0 4]);
            caxis([0 10]);
            title('Synchrosqueezed Transform for Raw LED3 PPG');  

            figure;
            contour(waveTime_s,wsstFiltLed3F,abs(wsstFiltLed3));
            hold on;
            plot(waveTime_s, wsstFiltLed3Ridge, 'r.');
            hold off;
            xlabel('Time, s'); ylabel('Frequency, Hz');
            grid on;
            ylim([0 4]);
            caxis([0 10]);
            title('Synchrosqueezed Transform for Filtered LED3 PPG'); 

            % plot 1st wsstridge for filtered red and nir as -
            % plot 2nd wsstridge for filtered red and nir as -.

            figure;
            plot(waveTime_s, wsstFiltLed2Ridge(:,1), 'b-');
            hold on; 
            plot(waveTime_s, wsstFiltLed2Ridge(:,2), 'b--');
            plot(waveTime_s, wsstFiltLed3Ridge(:,1), 'r-');
            plot(waveTime_s, wsstFiltLed3Ridge(:,2), 'r--');
            hold off;
            xlabel('Time, s'); ylabel('Frequency, Hz');
            ylim([0 2.5]);
            title('Synchrosqueezed Ridges - Raw and Filtered PPGs');
            legend('LED2 Pulse','LED3 Pulse','LED2 1st Harmonic', ...
                'LED3 1st Harmonic');

            figure;
            plot(waveTime_s, led2PulseBPM,'b-');
            hold on;
            plot(waveTime_s, led3PulseBPM,'r.');
            hold off;
            xlabel('Time, s'); ylabel('Estimated Heart Rate, BPM');
            title('Synchrosqueeze Ridges - Estimated Heart Rate');
            grid on;

            figure;
            plot(waveTime_s,led2PulseBPMh1,'b-');
            hold on;
            plot(waveTime_s,led3PulseBPMh1,'r.');
            hold off;
            xlabel('Time, s'); ylabel('Estimated 1H Heart Rate, BPM');
            title('Synchrosqueezed Ridges - Estimated HR 1st Harmonic');
            grid on;

        else

            for penalty = 6:4:11

                wsstLed2Ridge = wsstridge(wsstLed2, penalty, wsstLed2F, ...
                    'NumRidges',2);
                wsstFiltLed2Ridge = wsstridge(wsstFiltLed2, penalty, ...
                    wsstFiltLed2F, 'NumRidges',2);
                wsstLed3Ridge = wsstridge(wsstLed3, penalty, wsstLed3F, ...
                    'NumRidges',2);
                wsstFiltLed3Ridge = wsstridge(wsstFiltLed3, penalty, ...
                    wsstFiltLed3F, 'NumRidges',2);

                led2PulseF = wsstFiltLed2Ridge(:,1);
                led3PulseF = wsstFiltLed3Ridge(:,1);
                led2Pulse1stH = wsstFiltLed3Ridge(:,2);
                led3Pulse1stH = wsstFiltLed3Ridge(:,2);

                led2PulseBPM   = 60 * led2PulseF;
                led3PulseBPM   = 60 * led3PulseF;
                led2PulseBPMh1 = 60 * led2Pulse1stH;
                led3PulseBPMh1 = 60 * led2Pulse1stH;

                % do contour plots with ridge overlays

                figure;
                contour(waveTime_s,wsstLed2F,abs(wsstLed2));
                hold on;
                plot(waveTime_s, wsstLed2Ridge, 'b.');
                hold off;
                xlabel('Time, s'); ylabel('Frequency, Hz');
                grid on;
                ylim([0 4]);
                caxis([0 5]);
                title('Synchrosqueezed Transform - Raw LED2 PPG');

                figure;
                contour(waveTime_s,wsstFiltLed2F,abs(wsstFiltLed2));
                hold on;
                plot(waveTime_s, wsstFiltLed2Ridge, 'b.');
                hold off;
                xlabel('Time, s'); ylabel('Frequency, Hz');
                grid on;
                ylim([0 4]);
                caxis([0 5]);
                title('Synchrosqueezed Transform for Filtered LED2 PPG');            

                figure;
                contour(waveTime_s,wsstLed3F,abs(wsstLed3));
                hold on;
                plot(waveTime_s, wsstLed3Ridge, 'r.');
                hold off;
                xlabel('Time, s'); ylabel('Frequency, Hz');
                grid on;
                ylim([0 4]);
                caxis([0 10]);
                title('Synchrosqueezed Transform for Raw LED3 PPG');  

                figure;
                contour(waveTime_s,wsstFiltLed3F,abs(wsstFiltLed3));
                hold on;
                plot(waveTime_s, wsstFiltLed3Ridge, 'r.');
                hold off;
                xlabel('Time, s'); ylabel('Frequency, Hz');
                grid on;
                ylim([0 4]);
                caxis([0 10]);
                title('Synchrosqueezed Transform for Filtered LED3 PPG'); 

                % plot 1st wsstridge for filtered LED2 and LED3 as -
                % plot 2nd wsstridge for filtered LED2 and LED3 as -.

                figure;
                plot(waveTime_s, wsstFiltLed2Ridge(:,1), 'b-');
                hold on; 
                plot(waveTime_s, wsstFiltLed2Ridge(:,2), 'b--');
                plot(waveTime_s, wsstFiltLed3Ridge(:,1), 'r-');
                plot(waveTime_s, wsstFiltLed3Ridge(:,2), 'r--');
                hold off;
                xlabel('Time, s'); ylabel('Frequency, Hz');
                ylim([0 2.5]);
                title('Synchrosqueezed Ridges - Raw and Filtered PPGs');
                legend('LED2 Pulse','LED3 Pulse','LED2 1st Harmonic', ...
                    'LED3 1st Harmonic');

                figure;
                plot(waveTime_s,led2PulseBPM,'b-');
                hold on;
                plot(waveTime_s,led3PulseBPM,'r--');
                hold off;
                xlabel('Time, s'); ylabel('Estimated Heart Rate, BPM');
                title('Synchrosqueeze Ridges - Estimated Heart Rate');
                grid on;

                figure;
                plot(waveTime_s,led2PulseBPMh1,'b-');
                hold on;
                plot(waveTime_s,led3PulseBPMh1,'r--');
                hold off;
                xlabel('Time, s'); ylabel('Estimated 1H Heart Rate, BPM');
                title('Synchrosqueezed Ridges - Estimated HR 1st Harmonic');
                grid on;
        end

        end

    end

else
    
    fprintf('\tSkipping entire time-series wavelet HR estimation.\n');

end


%% try dddtree cleanup on a segment of unfiltered red ppg

doDDDTree = false;

if (doDDDTree)
    % next to be an even multiple of 2 ^ J, e.g.: 8192, 16384
    newLed2  = led2_noAmb_lessAccel(150:150+2047);         
    dt      = 1/fix(Fs);
    t       = 0:dt:(length(newLed2)*dt)-dt;
    J       = 6;
    dtcplx1 = dddtree('cplxdt', newLed2, J, 'dtf3');
    QQQ     = dtcplx1.cfs{2}(:,1,1) + 1i * dtcplx1.cfs{2}(:,1,2);
    QQQQ    = dtcplx1.cfs{3}(:,1,1) + 1i * dtcplx1.cfs{3}(:,1,2);

    figure;
    figD3Tree1b = subplot(211);
    stem(real(QQQ),'b-');
    figD3Tree1c = subplot(212);
    stem(real(QQQQ),'b-');
    % figD3Tree1a = subplot(311);
    % plot(t, newLed2,'b-');

    % linkaxes([ figD3Tree1a, figD3Tree1b, figD3Tree1c], 'x');
    linkaxes([ figD3Tree1b, figD3Tree1c], 'x');
else
    fprintf('\tSkipping dddtree cleanup.\n\n');
end

%% try some adaptive filters...


runAdaptiveFilter = 0;

if(runAdaptiveFilter)

    % a least mean squares adaptive filter...

    blms = dsp.BlockLMSFilter(250,250);
    blms.StepSize = 0.01;
    blms.WeightsOutputPort = false;
    filt = dsp.FIRFilter;
    filt.Numerator = fir1(13,[.5, .75]);

    [y, err] = blms(led2', led4');


    % a freq domain active filter
    mu  = 0.1;
    fdaf = dsp.FrequencyDomainAdaptiveFilter('Length', 250, 'StepSize', mu);
    [ led2AdaptiveFilter , err ] = fdaf(led2Norm, led3Norm);

end

%% creating a windowSlice routine

doWindowSliceWork = false;

if (doWindowSliceWork)
    
    windowSlice(1)=1;

    % specify the size of window slice in seconds

    windowSize = 5;        % this is five seconds

    % samplingIncrement = decimated fraction of a second each sample represents

    slidingWindowIncrement = windowSize / samplingIncrement;

    % compute the number of windowSlices in this data set

    for i = 1+slidingWindowIncrement:slidingWindowIncrement:numSamples; 
        windowSlice(end+1)=i;
        fprintf('\tComputing number of window slices...%d\n',i);
    end


    % compute the number of window slices to work with
    numWindowSlices = length(windowSlice);

    % put some console space between what happened and what's gonna happen
    fprintf('\n');

    % now do the work on each of those slices...
    for i = 1:1:numWindowSlices

        if ( (windowSlice(i) + slidingWindowIncrement ) > numSamples)

            fprintf('\tGoing past end of data. Stopping here...\n');
            return;

        else

            fprintf('\tWorking on +%d samples starting at windowSlice(%d)\n', ...
            slidingWindowIncrement, windowSlice(i) );

            % do more work on this set of samples... e.g.: FFT, wavelets, etc.

            sampleStart = windowSlice(i);
            sampleEnd = windowSlice(i) + slidingWindowIncrement;

            % wavelet start here

            doManualWSST = true;
            useDecimated = true;

            if (doManualWSST)

                % do some manual wsst stuff for plotting contours of least
                % penalized WSST results

                % build wsst for signals of interest

                if (useDecimated)

                    % do both raw and filtered (band-passed) data
                    [wsstLed1, wsstLed1F]         = ...
                        wsst(led1_25Hz(1,sampleStart:sampleEnd),Fs);
                    [wsstLed2, wsstLed2F]         = ...
                        wsst(led2_25Hz(1,sampleStart:sampleEnd),Fs);
                    [wsstLed3, wsstLed3F]         = ...
                        wsst(led3_25Hz(1,sampleStart:sampleEnd),Fs);
                    [wsstLed4, wsstLed4F]         = ...
                        wsst(led4_25Hz(1,sampleStart:sampleEnd),Fs);  

                    [wsstFiltLed1, wsstFiltLed1F] = ....
                        wsst(filtered_led1_25Hz(1,sampleStart:sampleEnd),Fs);
                    [wsstFiltLed2, wsstFiltLed2F] = ...
                        wsst(filtered_led2_25Hz(1,sampleStart:sampleEnd),Fs);
                    [wsstFiltLed3, wsstFiltLed3F] = ...
                        wsst(filtered_led3_25Hz(1,sampleStart:sampleEnd),Fs);
                    [wsstFiltLed4, wsstFiltLed4F] = ...
                        wsst(filtered_led4_25Hz(1,sampleStart:sampleEnd),Fs);

                else

                    [wsstLed1, wsstLed1F]         = wsst(led1,opticsFs);
                    [wsstLed2, wsstLed2F]         = wsst(led2,opticsFs);
                    [wsstLed3, wsstLed3F]         = wsst(led3,opticsFs);
                    [wsstLed4, wsstLed4F]         = wsst(led4,opticsFs); 
                    [wsstFiltLed1, wsstFiltLed1F] = wsst(bpLed1,opticsFs);
                    [wsstFiltLed2, wsstFiltLed2F] = wsst(bpLed2,opticsFs);
                    [wsstFiltLed3, wsstFiltLed3F] = wsst(bpLed3,opticsFs);
                    [wsstFiltLed4, wsstFiltLed4F] = wsst(bpLed4,opticsFs);            

                end

                % build wsst ridges for signals of interest, using 'NumRidges'
                % set to 2, to get primary signal and 1st harmonic
                % and run a loop (for now) to explore the effect of 'penalty'

                doMultiPenaltyTest = 0;

                if doMultiPenaltyTest == 0

                    penalty = 10;

                    wsstLed2Ridge = wsstridge(wsstLed2, penalty, wsstLed2F, ...
                        'NumRidges',2);
                    wsstFiltLed2Ridge = wsstridge(wsstFiltLed2, penalty, ...
                        wsstFiltLed2F, 'NumRidges',2);
                    wsstLed3Ridge = wsstridge(wsstLed3 , penalty, wsstLed3F, ...
                        'NumRidges',2);
                    wsstFiltLed3Ridge = wsstridge(wsstFiltLed3, penalty, ...
                        wsstFiltLed3F, 'NumRidges',2);

                    led2PulseF = wsstFiltLed2Ridge(:,1);
                    led3PulseF = wsstFiltLed3Ridge(:,1);
                    led2Pulse1stH = wsstFiltLed2Ridge(:,2);
                    led3Pulse1stH = wsstFiltLed3Ridge(:,2);

                    led2PulseBPM   = 60 * led2PulseF;
                    led3PulseBPM   = 60 * led3PulseF;
                    led2PulseBPMh1 = 60 * led2Pulse1stH;
                    led3PulseBPMh1 = 60 * led3Pulse1stH;

                    % do contour plots with ridge overlays

                    if (useDecimated)
                       waveTime_s = time_s; 
                    else
                        waveTime_s = opticsTime_s;
                    end

                    figure;
                    % contour(waveTime_s, wsstLed2F, abs(wsstLed2));
                    contour(sampleStart:sampleEnd, wsstLed2F, abs(wsstLed2));
                    hold on;
                    plot(sampleStart:sampleEnd, wsstLed2Ridge, 'b.');
                    hold off;
                    xlabel('Time, s'); ylabel('Frequency, Hz');
                    grid on;
                    ylim([0 4]);
                    caxis([0 5]);
                    title('Synchrosqueezed Transform - Raw LED2 PPG');

                    figure;
                    contour(sampleStart:sampleEnd,wsstFiltLed2F,abs(wsstFiltLed2));
                    hold on;
                    plot(sampleStart:sampleEnd, wsstFiltLed2Ridge, 'b.');
                    hold off;
                    xlabel('Time, s'); ylabel('Frequency, Hz');
                    grid on;
                    ylim([0 4]);
                    caxis([0 5]);
                    title('Synchrosqueezed Transform for Filtered LED2 PPG');            

                    figure;
                    contour(sampleStart:sampleEnd,wsstLed3F,abs(wsstLed3));
                    hold on;
                    plot(sampleStart:sampleEnd, wsstLed3Ridge, 'r.');
                    hold off;
                    xlabel('Time, s'); ylabel('Frequency, Hz');
                    grid on;
                    ylim([0 4]);
                    caxis([0 10]);
                    title('Synchrosqueezed Transform for Raw LED3 PPG');  

                    figure;
                    contour(sampleStart:sampleEnd,wsstFiltLed3F,abs(wsstFiltLed3));
                    hold on;
                    plot(sampleStart:sampleEnd, wsstFiltLed3Ridge, 'r.');
                    hold off;
                    xlabel('Time, s'); ylabel('Frequency, Hz');
                    grid on;
                    ylim([0 4]);
                    caxis([0 10]);
                    title('Synchrosqueezed Transform for Filtered LED3 PPG'); 

                    % plot 1st wsstridge for filtered red and nir as -
                    % plot 2nd wsstridge for filtered red and nir as -.

                    figure;
                    plot(sampleStart:sampleEnd, wsstFiltLed2Ridge(:,1), 'b-');
                    hold on; 
                    plot(sampleStart:sampleEnd, wsstFiltLed2Ridge(:,2), 'b--');
                    plot(sampleStart:sampleEnd, wsstFiltLed3Ridge(:,1), 'r-');
                    plot(sampleStart:sampleEnd, wsstFiltLed3Ridge(:,2), 'r--');
                    hold off;
                    xlabel('Time, s'); ylabel('Frequency, Hz');
                    ylim([0 2.5]);
                    title('Synchrosqueezed Ridges - Raw and Filtered PPGs');
                    legend('LED2 Pulse','LED3 Pulse','LED2 1st Harmonic', ...
                        'LED3 1st Harmonic');

                    figure;
                    plot(sampleStart:sampleEnd, led2PulseBPM,'b-');
                    hold on;
                    plot(sampleStart:sampleEnd, led3PulseBPM,'r.');
                    hold off;
                    xlabel('Time, s'); ylabel('Estimated Heart Rate, BPM');
                    title('Synchrosqueeze Ridges - Estimated Heart Rate');
                    grid on;

                    figure;
                    plot(sampleStart:sampleEnd,led2PulseBPMh1,'b-');
                    hold on;
                    plot(sampleStart:sampleEnd,led3PulseBPMh1,'r.');
                    hold off;
                    xlabel('Time, s'); ylabel('Estimated 1H Heart Rate, BPM');
                    title('Synchrosqueezed Ridges - Estimated HR 1st Harmonic');
                    grid on;

                else

                    for penalty = 6:4:11

                        wsstLed2Ridge = wsstridge(wsstLed2, penalty, wsstLed2F, ...
                            'NumRidges',2);
                        wsstFiltLed2Ridge = wsstridge(wsstFiltLed2, penalty, ...
                            wsstFiltLed2F, 'NumRidges',2);
                        wsstLed3Ridge = wsstridge(wsstLed3, penalty, wsstLed3F, ...
                            'NumRidges',2);
                        wsstFiltLed3Ridge = wsstridge(wsstFiltLed3, penalty, ...
                            wsstFiltLed3F, 'NumRidges',2);

                        led2PulseF = wsstFiltLed2Ridge(:,1);
                        led3PulseF = wsstFiltLed3Ridge(:,1);
                        led2Pulse1stH = wsstFiltLed3Ridge(:,2);
                        led3Pulse1stH = wsstFiltLed3Ridge(:,2);

                        led2PulseBPM   = 60 * led2PulseF;
                        led3PulseBPM   = 60 * led3PulseF;
                        led2PulseBPMh1 = 60 * led2Pulse1stH;
                        led3PulseBPMh1 = 60 * led2Pulse1stH;

                        % do contour plots with ridge overlays

                        figure;
                        contour(waveTime_s,wsstLed2F,abs(wsstLed2));
                        hold on;
                        plot(waveTime_s, wsstLed2Ridge, 'b.');
                        hold off;
                        xlabel('Time, s'); ylabel('Frequency, Hz');
                        grid on;
                        ylim([0 4]);
                        caxis([0 5]);
                        title('Synchrosqueezed Transform - Raw LED2 PPG');

                        figure;
                        contour(waveTime_s,wsstFiltLed2F,abs(wsstFiltLed2));
                        hold on;
                        plot(waveTime_s, wsstFiltLed2Ridge, 'b.');
                        hold off;
                        xlabel('Time, s'); ylabel('Frequency, Hz');
                        grid on;
                        ylim([0 4]);
                        caxis([0 5]);
                        title('Synchrosqueezed Transform for Filtered LED2 PPG');            

                        figure;
                        contour(waveTime_s,wsstLed3F,abs(wsstLed3));
                        hold on;
                        plot(waveTime_s, wsstLed3Ridge, 'r.');
                        hold off;
                        xlabel('Time, s'); ylabel('Frequency, Hz');
                        grid on;
                        ylim([0 4]);
                        caxis([0 10]);
                        title('Synchrosqueezed Transform for Raw LED3 PPG');  

                        figure;
                        contour(waveTime_s,wsstFiltLed3F,abs(wsstFiltLed3));
                        hold on;
                        plot(waveTime_s, wsstFiltLed3Ridge, 'r.');
                        hold off;
                        xlabel('Time, s'); ylabel('Frequency, Hz');
                        grid on;
                        ylim([0 4]);
                        caxis([0 10]);
                        title('Synchrosqueezed Transform for Filtered LED3 PPG'); 

                        % plot 1st wsstridge for filtered LED2 and LED3 as -
                        % plot 2nd wsstridge for filtered LED2 and LED3 as -.

                        figure;
                        plot(waveTime_s, wsstFiltLed2Ridge(:,1), 'b-');
                        hold on; 
                        plot(waveTime_s, wsstFiltLed2Ridge(:,2), 'b--');
                        plot(waveTime_s, wsstFiltLed3Ridge(:,1), 'r-');
                        plot(waveTime_s, wsstFiltLed3Ridge(:,2), 'r--');
                        hold off;
                        xlabel('Time, s'); ylabel('Frequency, Hz');
                        ylim([0 2.5]);
                        title('Synchrosqueezed Ridges - Raw and Filtered PPGs');
                        legend('LED2 Pulse','LED3 Pulse','LED2 1st Harmonic', ...
                            'LED3 1st Harmonic');

                        figure;
                        plot(waveTime_s,led2PulseBPM,'b-');
                        hold on;
                        plot(waveTime_s,led3PulseBPM,'r--');
                        hold off;
                        xlabel('Time, s'); ylabel('Estimated Heart Rate, BPM');
                        title('Synchrosqueeze Ridges - Estimated Heart Rate');
                        grid on;

                        figure;
                        plot(waveTime_s,led2PulseBPMh1,'b-');
                        hold on;
                        plot(waveTime_s,led3PulseBPMh1,'r--');
                        hold off;
                        xlabel('Time, s'); ylabel('Estimated 1H Heart Rate, BPM');
                        title('Synchrosqueezed Ridges - Estimated HR 1st Harmonic');
                        grid on;
                    end
                end
            end

            % wavelet end here

        end

    end

end     % doWindowSliceWork

%% hey program, let the whole world know you're just done!

fprintf('\n---=== Ending FaunaTagDataWorkup! ===---\n\n');

return;

%% do some input tests

prompt = '\n(b)ack (f)orward (m)ark (r)ange (c)omment (e)xclude)';

%% Experimental kinematic heart rate & optics comparison

% partition the data into a time segment, 
% here 600 samples / 25 Hz = 24 seconds

timeStartSeconds    = 525.0;
timeEndSeconds      = 535.00;
sampleStart         = round(timeStartSeconds * 25);
sampleEnd           = round(timeEndSeconds * 25);

odbaWork = odba_25Hz(1,sampleStart:sampleEnd);
led2Work = led2_25Hz(1,sampleStart:sampleEnd);
led3Work = led3_25Hz(1,sampleStart:sampleEnd);
led1Work = led1_25Hz(1,sampleStart:sampleEnd);
led4Work = led4_25Hz(1,sampleStart:sampleEnd);

% run some loops to get a window'd version of norm LED2 LED3 and ODBA
% this avoids the problem with ±2,000,000 for optics, better accounts for
% local mins and maxes

for i = 1:length(led2Work) 
    led2Norm(1,i) = ( ( (led2Work(1,i) ) - mean(led2Work(1,:))) / std2(led2Work(1,:)) ); 
end

for i = 1:length(led3Work) 
    led3Norm(1,i) = ( ( (led3Work(1,i) ) - mean(led3Work(1,:))) / std2(led3Work(1,:)) ); 
end

for i = 1:length(led1Work) 
    led1Norm(1,i) = ( ( (led1Work(1,i) ) - mean(led1Work(1,:))) / std2(led1Work(1,:)) ); 
end

for i = 1:length(led4Work) 
    led4Norm(1,i) = ( ( (led4Work(1,i) ) - mean(led4Work(1,:))) / std2(led4Work(1,:)) ); 
end

for i = 1:length(odbaWork) 
    odbaNorm(1,i) = ( ( (odbaWork(1,i) ) - mean(odbaWork(1,:))) / std2(odbaWork(1,:)) ); 
end

% now build first & second derivative estimates from windowed odbaNorm
d1odbaNorm = diff(odbaNorm) / (Fs * 0.1) ; d1odbaNorm(end+1, 1) = 0;
d2odbaNorm = diff(d1odbaNorm) / (Fs * 0.1) ; d2odbaNorm(end+1, 1) = 0;

[pksD1, locsD1, widthsD1, promsD1] = findpeaks(d1odbaNorm(1,:), Fs, 'MinPeakDistance', 0.25);
[pksD2, locsD2, widthsD2, promsD2] = findpeaks(d2odbaNorm(1,:), Fs, 'MinPeakDistance', 0.25);

% Examine the extent of correlation between locsD1 and locsD2
if ( length(locsD2) == length(locsD1) )
    corrD1D2 = corr(locsD1', locsD2');
    fprintf('\tlengths for locsD1 & locsD2 match. corr = %d.\n', corrD1D2);
else
    fprintf('\tlocsD1 != locsD2. Skipping corr calc.\n');
end


% make some pretty pictures of all this...

figure;
plot(led2Norm, 'b-'); 
hold on; 
plot(led3Norm, 'r-');
plot(led1Norm, 'k--');
plot(led4Norm, 'g--');
plot(odbaNorm, 'k-');
plot(d2odbaNorm(:,1)', 'g-.');
plot(locsD1 * 25, pksD1, 'mv', 'MarkerFaceColor', 'm'); 
plot(locsD2 * 25, pksD2, 'gv', 'MarkerFaceColor', 'g'); 
hold off;
xlabel('Samples, 25 Hz');
ylabel('Windowed Norms');
title('Norm LED2 LED3 w/ Norm ODBA & D1 + D2 HR markers');
grid;
% legend('LED2','LED3', 'ODBA','d1ODBA','d2ODBA');


figure;
plot( (abs(led2Norm) - abs(led1Norm) ), 'b-'); 
hold on; 
plot( (abs(led3Norm) - abs(led1Norm) ), 'r-');
% plot(led1Norm, 'k--');
% plot(led4Norm, 'g--');
%plot(odbaNorm, 'k-');
plot(d2odbaNorm(1,:)', 'k-');
plot(locsD1 * 25, pksD1, 'mv', 'MarkerFaceColor', 'm'); 
plot(locsD2 * 25, pksD2, 'gv', 'MarkerFaceColor', 'g'); 
hold off;
xlabel('Samples, 25 Hz');
ylabel('Windowed Norms');
title('Norm LED2 LED3 w/ Norm ODBA & D1 + D2 HR markers');
grid;
% legend('LED2','LED3', 'ODBA','d1ODBA','d2ODBA');

%% do a wsst of normOdba_25Hz

[wsstOdbaNorm, wsstOdbaNorm_F] = wsst(odbaNorm',Fs);

penalty = 10;

wsstOdbaNormRidge = wsstridge(wsstOdbaNorm, penalty, wsstOdbaNorm_F, ...
    'NumRidges',2);

figure;
%contour(waveTime_s,wsstOdbaNorm_F,abs(wsstOdbaNorm));
contour(abs(wsstOdbaNorm));
hold on;
plot(wsstOdbaNormRidge, 'b.');
hold off;
xlabel('Samples, 25 Hz'); ylabel('Frequency, Hz');
grid on;
ylim([0 4]);
caxis([0 5]);
title('Synchrosqueezed Transform - Raw NormOdba PPG');