%% Optical gaussian filtering and peak finder scratch space

%   Written by Dave Haas between 17-26 June 2021
%   
%   Once this works, integrate into the main FaunaDataAnalysis routine.

clc;
clear;

%% load some data in...

global TAG_PATHS;

tag = 'tt21_134a';

rawFileName = sprintf('%s/%sraw.mat', TAG_PATHS.RAW, tag);
prhFileName = sprintf('%s/%sprh.mat', TAG_PATHS.PRH, tag);

if (exist(rawFileName,'file'))
    fprintf('Specified RAW file was found. Loading!\n');
    load(rawFileName);
else
    fprintf('Specified RAW file was not found. Check it out and try again!\n');
    return;
end

if (exist(prhFileName,'file'))
    fprintf('Specified PRH file was found. Loading!\n');
    load(prhFileName);
else
    fprintf('Specified PRH file was not found. Check it out and try again!\n');
    return;
end

if ( exist('TIME','var') )      % if TIME struct is present, use 
   time = TIME.time;            % the TIME.time and TIME.time_s
   time_s = TIME.time_s;        % as default time for x-axis plotting
end

%% gaussian filter experimenter

% compare LED2 with gaussian-filtered LED2

% first perform the reflecance inversion on the data...

% rLedX = reflectance correction via inversion

rLed1 = led1 * -1;      
rLed2 = led2 * -1;
rLed3 = led3 * -1;
rLed4 = led4 * -1;

% get min values of each for DC offset of rLedX

min_rLed1 = min(rLed1);
min_rLed2 = min(rLed2);
min_rLed3 = min(rLed3);
min_rLed4 = min(rLed4);

rdcoLed1 = rLed1 - min_rLed1;
rdcoLed2 = rLed2 - min_rLed2;
rdcoLed3 = rLed3 - min_rLed3;
rdcoLed4 = rLed4 - min_rLed4;

rdcoLed1_25Hz = decimate(rdcoLed1, 10);
rdcoLed2_25Hz = decimate(rdcoLed2, 10);
rdcoLed3_25Hz = decimate(rdcoLed3, 10);
rdcoLed4_25Hz = decimate(rdcoLed4, 10);


%% Do a bandpass filter on rdco 250 Hz & 25 Hz optical data (for HR freqs)

fprintf('\tApplying a band-pass filter with the following properties:\n');

% Bandpass example is a 500 to 560 Hz bandpass, 
% [A,B,C,D] = butter(10,[500 560]/750); 
% d = designfilt('bandpassiir','FilterOrder',20, ...
%    'HalfPowerFrequency1',500,'HalfPowerFrequency2',560, ...
%    'SampleRate',1500);
%    //  where 10 = 1/2 filter order?
%    //       750 = 1/2 sampling frequency   

% so for raw optical data, where opticsFs = 250 & 20th order BW filt...
bwFiltOrder = 6;        % 20th order butterworth filter

% experiment with some low and high-pass frequencies
% my initial guesses for some good band-pass freqs
lowCut      = 0.3333;   % 0.3333 Hz = ~20 beats per minute
% lowCut      = 0.5;      % 0.5 Hz    = ~30 beats per minute
highCut     = 3.5;      % 3.5 Hz    = 210 beats per minute

% numbers suggested by the fieldtriptoolbox folks - good for O2 sat/desat
% lowCut      = 0.01;       % 0.01 Hz   = fluctuations over 10 seconds
% highCut     = 0.1;        % 0.1 Hz    = fluctuations over 100 seconds

% try some different values for nyquistFs...
% here's an explainer of some prelim results / looks at different values
testFs = 8;             % 8 produces big pulses around odba HR in bpLed2
% testFs = 9;             % 9 produces big pulses around odba HR in bpLed2
% testFs = 25;            % 25 produces some overlaps in bpLed3 but washes
                          % out the bpLed2 signal. Maybe this suggests
                          % using separate filters for LED2 and LED3?
% testFs = 250;           % the default value                          
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

% Convert the state-space representation to second-order sections. 
% Visualize the frequency responses using fvtool.

sos = ss2sos(A,B,C,D);
fvt = fvtool(sos, d, 'Fs', nyquistFs);
legend(fvt,'butter','designfilt');

%% Apply the filter to reflectance LEDs with coarse DC offset

% bpRdcoLed1 = filtfilt(d, rdcoLed1);
% bpRdcoLed2 = filtfilt(d, rdcoLed2);
% bpRdcoLed3 = filtfilt(d, rdcoLed3);
% bpRdcoLed4 = filtfilt(d, rdcoLed4);

bpRdcoLed1 = filtfilt(d, OPTICS.rLed1);
bpRdcoLed2 = filtfilt(d, OPTICS.rLed2);
bpRdcoLed3 = filtfilt(d, OPTICS.rLed3);
bpRdcoLed4 = filtfilt(d, OPTICS.rLed4);

% bpRdcoLed1_25Hz = filtfilt(d, rdcoLed1_25Hz);
% bpRdcoLed2_25Hz = filtfilt(d, rdcoLed2_25Hz);
% bpRdcoLed3_25Hz = filtfilt(d, rdcoLed3_25Hz);
% bpRdcoLed4_25Hz = filtfilt(d, rdcoLed4_25Hz');

bpRdcoLed1_25Hz = filtfilt(d, O.rLed1);
bpRdcoLed2_25Hz = filtfilt(d, O.rLed2);
bpRdcoLed3_25Hz = filtfilt(d, O.rLed3);
bpRdcoLed4_25Hz = filtfilt(d, O.rLed4);

bpLed1 = filtfilt(d, O.led1);
bpLed2 = filtfilt(d, O.led2);
bpLed3 = filtfilt(d, O.led3);
bpLed4 = filtfilt(d, O.led4);

%% Plot measured LED reflectance with and without DC offset

figure;

r1 = subplot(211);
plot(OPTICS.time_s, OPTICS.rLed1, 'k-');
hold on;
plot(OPTICS.time_s, OPTICS.rLed2, 'b-');
plot(OPTICS.time_s, OPTICS.rLed3, 'r-');
plot(OPTICS.time_s, OPTICS.rLed4, 'g-');
hold off;
xlabel('Time (local)');
ylabel('Intensity');
title('Reflectance with DC Offset');
grid;

r2 = subplot(212);
plot(OPTICS.time_s, OPTICS.led1, 'k-');
hold on; 
plot(OPTICS.time_s, OPTICS.led2, 'b-.');
plot(OPTICS.time_s, OPTICS.led3, 'r-.');
plot(OPTICS.time_s, OPTICS.led4, 'g-.');
hold off; 
xlabel('Time (local)');
ylabel('Intensity');
title('Raw intensities of LEDs 1-4');
grid;

linkaxes([r1,r2],'x');

figure;
plot(time, bpLed1, 'k-');
hold on;
plot(time, bpLed2, 'b-');
plot(time, bpLed3, 'r-');
plot(time, bpLed4, 'g-');
hold off;
xlabel('Time, seconds');
ylabel('Band-passed intensity');
title('Band-passed raw intensities - all optical channels');
grid;

%% Build and apply the Gaussian filter for OD analysis

% set a 100 ms windowSize to use for gaussian filtering
optGaussianWindowSize = opticsFs / 10;
kinGaussianWindowSize = movementFs / 4;

optB = ones(1, optGaussianWindowSize) / optGaussianWindowSize;
kinB = ones(1, kinGaussianWindowSize) / kinGaussianWindowSize;

A = [1];

gfLed1 = filtfilt(optB, A, OPTICS.led1);
gfLed2 = filtfilt(optB, A, OPTICS.led2);
gfLed3 = filtfilt(optB, A, OPTICS.led3);
gfLed4 = filtfilt(optB, A, OPTICS.led4);

gf_ODLed1 = filtfilt(optB, A,(log(max(OPTICS.led1)./OPTICS.led1)));
gf_ODLed2 = filtfilt(optB, A,(log(max(OPTICS.led2)./OPTICS.led2)));
gf_ODLed3 = filtfilt(optB, A,(log(max(OPTICS.led3)./OPTICS.led3)));
gf_ODLed4 = filtfilt(optB, A,(log(max(OPTICS.led4)./OPTICS.led4)));

gfAx = filtfilt(kinB, A, KINEMATICS.ax);
gfAy = filtfilt(kinB, A, KINEMATICS.ay);
gfAz = filtfilt(kinB, A, KINEMATICS.az);
gfGx = filtfilt(kinB, A, KINEMATICS.gx);
gfGy = filtfilt(kinB, A, KINEMATICS.gy);
gfGz = filtfilt(kinB, A, KINEMATICS.gz);

sgoGfAz = gfAz - (sgolayfilt(KINEMATICS.az,30,401)); 

%% Plot to have a look at gaussian filtered signals

figure;

g1 = subplot(211);
plot(OPTICS.time_s, gfLed1, 'k');
hold on;
plot(OPTICS.time_s, gfLed2, 'b');
plot(OPTICS.time_s, gfLed3, 'r');
plot(OPTICS.time_s, gfLed4, 'g');
hold off;
xlabel('Time, seconds)');
ylabel('Intensity');
title('Raw optical intensities smoothed with 100 ms Gaussian filtering');
legend('ambient1','1050 nm','1200 nm','ambient2');
grid;

g2 = subplot(212);
plot(OPTICS.time_s, gf_ODLed1, 'k');
hold on;
plot(OPTICS.time_s, gf_ODLed2, 'b');
plot(OPTICS.time_s, gf_ODLed3, 'r');
plot(OPTICS.time_s, gf_ODLed4, 'g');
hold off;
xlabel('Time, seconds');
ylabel('OD');
title('Optical Density (OD) smoothed with 100 ms Gaussian filtering');
legend('ambient1','1050 nm','1200 nm','ambient2');
grid;

linkaxes([g1,g2],'x');

%% experiment with some mean- and standard deviation-based smoothing


g1_normRdcoLed2 = normalize(g1_rdcoLed2);
g2_normRdcoLed2 = normalize(g2_rdcoLed2);
g5_normRdcoLed2 = normalize(g5_rdcoLed2);

figure;
plot(OPTICS.time_s, OPTICS.rLed2, 'b-');
hold on;
plot(OPTICS.time_s, g1_rdcoLed2, 'b--');
plot(OPTICS.time_s, g2_rdcoLed2, 'b-.');
plot(OPTICS.time_s, g5_rdcoLed2, 'b:');
hold off;
title('Gaussian filtered LED2');
grid;

figure;
plot(OPTICS.time_s, OPTICS.nrLed2, 'b-');
hold on;
plot(OPTICS.time_s, g1_normRdcoLed2, 'b--');
plot(OPTICS.time_s, g2_normRdcoLed2, 'b-.');
plot(OPTICS.time_s, g5_normRdcoLed2, 'b:');
hold off;
title('Normalized gaussian filtered LED2');
grid;


%% examine effectiveness of sgolayfilt for cleaning up HR LED & ODBA data

% do some work with non-decimated optical data...

led2sgoFilt = sgolayfilt(OPTICS.rLed2, 6, 375);
led3sgoFilt = sgolayfilt(OPTICS.rLed3, 6, 375);

figure;
plot(OPTICS.time_s, OPTICS.rLed2-led2sgoFilt,'b-');
hold on;
plot(OPTICS.time_s, OPTICS.rLed3-led3sgoFilt,'r-');
hold off;
grid;

%% examine effectiveness of sgolayfilt for cleaning up HR LED & ODBA data

% do some work with decimated optical data...

filterOrder = 4;
frameLength = 75;

% d2bpOdba_25Hz = diff(bpOdba_25Hz, 2);
% d2bpOdba_25Hz(1,end+1:end+2) = 0;     % tack on end zeros for equal length
K.d1odba = diff(K.odba, 1);
K.d1odba(1,end+1) = 0;     % tack on end zeros for equal length

K.d2odba = diff(K.odba, 2);
K.d2odba(1,end+1:end+2) = 0;     % tack on end zeros for equal length

led1sgoFilt = sgolayfilt(O.rLed1, filterOrder, frameLength);
led2sgoFilt = sgolayfilt(O.rLed2, filterOrder, frameLength);
led3sgoFilt = sgolayfilt(O.rLed3, filterOrder, frameLength);
led4sgoFilt = sgolayfilt(O.rLed4, filterOrder, frameLength);
odbasgoFilt = sgolayfilt(K.odba, filterOrder, frameLength);
d1odbasgoFilt = sgolayfilt(K.d1odba, filterOrder, frameLength);
d2odbasgoFilt = sgolayfilt(K.d2odba, filterOrder, frameLength);

sgo_rLed1= O.rLed1 - led1sgoFilt;
sgo_rLed2 = O.rLed2 - led2sgoFilt;
sgo_rLed3 = O.rLed3 - led3sgoFilt;
sgo_rLed4 = O.rLed4 - led4sgoFilt;
sgo_odba = K.odba - odbasgoFilt;
sgo_d1odba = K.d1odba - d1odbasgoFilt;
sgo_d2odba = K.d1odba - d2odbasgoFilt;

figure;

f1 = subplot(211);
plot(time_s, sgo_rLed1, 'k--');
hold on;
plot(time_s, sgo_rLed2, 'b--');
plot(time_s, sgo_rLed3, 'r--');
plot(time_s, sgo_rLed4, 'g--');
plot(time_s, O.hrBwLed2, 'b-');
plot(time_s, O.hrBwLed3, 'r-');
hold off;
xlabel('Time (local)');
ylabel('Intensity');
figTitle = sprintf('Savitzky-Golay optics · 25 Hz, order: %d, frame: %d', ... 
    filterOrder, frameLength);
title(figTitle);
grid;
legend('Ambient','sgo 1050 nm','sgo 1200 nm','AltAmbient','bp 1050 nm','bp 1200 nm','Location','southeast');

f2 = subplot(212);
plot(time_s, K.odba, 'k-');
hold on;
plot(time_s, sgo_odba, 'b-');
plot(time_s, sgo_d1odba, 'r-');
plot(time_s, sgo_d2odba, 'g-');
hold off;
xlabel('Time (local)');
ylabel('Movement');
figTitle = sprintf('Savitzky-Golay ODBA · 25 Hz, order: %d, frame: %d', ... 
    filterOrder, frameLength);
title(figTitle);
grid;
legend('bpOdba','sgoOdba','Location','southeast');

linkaxes([f1,f2],'x');

%%

windowSize = 2; % in seconds
windowSamples = windowSize * Fs;

sampleStart = 13675;
sampleEnd = sampleStart + windowSamples;

d1_bpOdba_25Hz = diff(bpOdba_25Hz,1);
d1_bpOdba_25Hz(end+1) = 0;  % same length as bpLed2 & bpOdba
d1_bpOdba_25Hz = d1_bpOdba_25Hz';

Fn = Fs / 2;    %nyquist                             
L = numel(bpRdcoLed2_25Hz(sampleStart:sampleEnd));
timeScale = linspace(0, L, L)/Fs;
N = 2^nextpow2(L);
bpRdcoLed2FTs = fft(bpRdcoLed2_25Hz(sampleStart:sampleEnd),N)/L;
bpOdbaFTs = fft(bpOdba_25Hz(sampleStart:sampleEnd),N)/L;
d1_bpOdbaFTs = fft(d1_bpOdba_25Hz(sampleStart:sampleEnd),N)/L;
Fv = linspace(0, 1, N/2+1)*Fn;
Iv = 1:numel(Fv);


figure;

figA1 = subplot(321);
plot(time(sampleStart:sampleEnd), bpRdcoLed2_25Hz(sampleStart:sampleEnd), 'b-');
grid;

figA2 = subplot(323);
plot(time(sampleStart:sampleEnd), bpOdba_25Hz(sampleStart:sampleEnd), 'k-');
grid;

figA3 = subplot(325);
plot(time(sampleStart:sampleEnd), d1_bpOdba_25Hz(sampleStart:sampleEnd), 'r-');
grid;

figA4 = subplot(322);
plot(Fv, abs(bpRdcoLed2FTs(Iv))*2, 'b-');
grid;

figA5 = subplot(324);
plot(Fv, abs(bpOdbaFTs(Iv))*2, 'k-');
grid;

figA6 = subplot(326);
plot(Fv, abs(d1_bpOdbaFTs(Iv))*2, 'r-');
grid;

linkaxes([figA1 figA2 figA3],'x');
linkaxes([figA4 figA5 figA6], 'x');

%% do peak detections on bpLed2, bpOdba and d1_bpOdba

figure; 

% do bpLed2_25Hz peak detection

sampleMean_bpRdcoLed2       = mean(bpRdcoLed2_25Hz(1,sampleStart:sampleEnd));
sampleSd_bpRdcoLed2         = std(bpRdcoLed2_25Hz(1,sampleStart:sampleEnd));
peakThreshold_bpRdcoLed2	= 0.75 * sampleSd_bpRdcoLed2;

subplot(311);

findpeaks(bpRdcoLed2_25Hz(sampleStart:sampleEnd), ...
    time_s(sampleStart:sampleEnd), ...
    'MinPeakProminence', peakThreshold_bpLed2, 'Annotate', 'extents');
title('Band-passed RDCO LED2 during breath-hold recovery');

% (1,sampleStart:sampleEnd)

[pks_bpRdcoLed2, locs_bpRdcoLed2, widths_bpRdcoLed2, proms_bpRdcoLed2] = ...
    findpeaks(bpRdcoLed2_25Hz(sampleStart:sampleEnd), ...
    time_s(sampleStart:sampleEnd), ...
    'MinPeakProminence', peakThreshold_bpLed2, 'Annotate', 'extents');

% do bpOdba_25Hz peak detection

sampleMean_bpOdba       = mean(bpOdba_25Hz(sampleStart:sampleEnd,1));
sampleSd_bpOdba         = std(bpOdba_25Hz(sampleStart:sampleEnd,1));
peakThreshold_bpOdba    = 0.9 * sd_bpOdba_25Hz;

subplot(312);

findpeaks(bpOdba_25Hz(sampleStart:sampleEnd,1), ...
    time_s(1,sampleStart:sampleEnd), ...
    'MinPeakProminence',peakThreshold_bpOdba, 'Annotate','extents');
title('Band-passed ODBA during breath-hold recovery');

[pks_bpOdba, locs_bpOdba, widths_bpOdba, proms_bpOdba] = ...
    findpeaks(bpOdba_25Hz(sampleStart:sampleEnd,1), ...
    time_s(1,sampleStart:sampleEnd), ...
    'MinPeakProminence',peakThreshold_bpOdba, 'Annotate','extents');

% do d1 bpOdba peak detection

sampleMean_d1bpOdba     = mean(d1_bpOdba_25Hz(sampleStart:sampleEnd));
sampleSd_d1bpOdba       = std(d1_bpOdba_25Hz(sampleStart:sampleEnd));
peakThreshold_d1bpOdba  = sampleMean_d1bpOdba + 2 * sampleSd_d1bpOdba;

subplot(313);

findpeaks(d1_bpOdba_25Hz(sampleStart:sampleEnd), ...
    time_s(sampleStart:sampleEnd), ...
    'MinPeakProminence',threshold,'Annotate','extents');
title('Band-passed 1st Derivative ODBA during breath-hold recovery');

[pks_d1bpOdba, locs_d1bpOdba, widths_d1bpOdba, prom_d1bpOdbas] = ...
    findpeaks(d1_bpOdba_25Hz(sampleStart:sampleEnd), ...
    time_s(sampleStart:sampleEnd),...
    'MinPeakProminence',threshold,'Annotate','extents');