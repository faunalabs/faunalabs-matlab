%% FaunaTag - tagOS csv data import
%       written by Dave Haas, 17 June 2020
%       uses csvproc from dtagtools by Mark Johnson, St. Andrews University

clear;

%% do the work...

% specify the path and filename here:




% specify pathName - was previously done statically by next line
% pathName = '/Users/dave/Documents/FaunaData/DH/';

useGoogle = true;

if (useGoogle)
    pathName = '/Volumes/GoogleDrive/Shared drives/FaunaData/';
else
    pathName = '/Users/dave/Documents/FaunaData/';
end

% specify subpath depending on whether its my data or animal data

animalData = true;

if (animalData)
    subPath = 'DQO/';
else
    subPath = 'DH/';
end

% FaunaTag 110 (reborn!) with Dave's pulse
% fileName = 'FaunaTag110/20210504/ft110_20210504_1551xx-10.txt';

% FaunaTag 111 with Hua's measurements
fileName = 'FaunaTag111/20210503/ft111_20210503_142825-10.txt';

% now, build a complete targetFile for open...
targetFile = sprintf('%s%s%s', pathName, subPath, fileName);

% load targetFile into S
S = csvproc(targetFile, []);

%% parse S into usable data components

% optics (relative insensity)           @ 250 Hz
led2 = str2double(S(1,:));
led3 = str2double(S(2,:));
led1 = str2double(S(3,:));
led4 = str2double(S(4,:));

num250HzSamples = length(led2);
numSecondsData = num250HzSamples / 250;

% accel (meters per second squared)     @ 100 Hz
% ax = str2double(S(5,:));
% ay = str2double(S(6,:));
% az = str2double(S(7,:));
ax = str2double(S(5,1:numSecondsData * 100));
ay = str2double(S(6,1:numSecondsData * 100));
az = str2double(S(7,1:numSecondsData * 100));


% gyro (degrees per second)             @ 100 Hz
% gx = str2double(S(8,:));
% gy = str2double(S(9,:));
% gz = str2double(S(10,:));
gx = str2double(S(8,1:numSecondsData * 100));
gy = str2double(S(9,1:numSecondsData * 100));
gz = str2double(S(10,1:numSecondsData * 100));

% prh (degree)                          @ 100 Hz
% pitch   = str2double(S(11,:));
% roll    = str2double(S(12,:));
% heading = str2double(S(13,:));
pitch   = str2double(S(11,1:numSecondsData * 100));
roll    = str2double(S(12,1:numSecondsData * 100));
heading = str2double(S(13,1:numSecondsData * 100));

% depth (meters)                        @ 100 Hz
% depth = str2double(S(14,:));
depth = str2double(S(14,1:numSecondsData * 100));

% temp (degrees Celsius)                @ 100 Hz
% temperature = str2double(S(15,:));
temperature = str2double(S(15,1:numSecondsData * 100));

% mag (degree)                                  @ 20 Hz
% mx = str2double(S(16,:));
% my = str2double(S(17,:));
% mz = str2double(S(18,:));
mx = str2double(S(16,1:numSecondsData * 20));
my = str2double(S(17,1:numSecondsData * 20));
mz = str2double(S(18,1:numSecondsData * 20));

% timestamp                                     @ 2 Hz
% sampletime = str2double(S(19,:));
sampletime = str2double(S(19,1:numSecondsData * 2));

% power measurements                            @ 2 Hz
% voltage             = str2double(S(20,:));
% current             = str2double(S(21,:));
% power               = str2double(S(22,:));
% chargeState         = str2double(S(23,:));
% remainingCapacity   = str2double(S(24,:));
voltage             = str2double(S(20,1:numSecondsData * 2));
current             = str2double(S(21,1:numSecondsData * 2));
power               = str2double(S(22,1:numSecondsData * 2));
chargeState         = str2double(S(23,1:numSecondsData * 2));
remainingCapacity   = str2double(S(24,1:numSecondsData * 2));


%% Diagnostic optical plots

figure;

led1Plot = subplot(411);
plot(led1, 'k-');
title('LED1');
ylabel('Intensity');

led2Plot = subplot(412);
plot(led2, 'b-');
title('LED2');
ylabel('Intensity');

led3Plot = subplot(413);
plot(led3, 'g-');
title('LED3');
ylabel('Intensity');

led4Plot = subplot(414);
plot(led4, 'r-');
title('LED4');
ylabel('Intensity');

linkaxes([led1Plot led2Plot led3Plot led4Plot], 'x');


%% Diagnostic kinematics plots

figure;

accelPlot = subplot(411);
plot(ax, 'b-');
hold on;
plot(ay, 'r-');
plot(az, 'g-');
hold off;
title('Accelerometer');
ylabel('m/s\^2');

gyroPlot = subplot(412);
plot(gx, 'b-');
hold on;
plot(gy, 'r-');
plot(gz, 'g-');
hold off;
title('Gyroscope');
ylabel('DPS');

magPlot = subplot(413);
plot(mx, 'b-');
hold on;
plot(my, 'r-');
plot(mz, 'g-');
hold off
title('Magnetometer');
ylabel('Degrees');

prhPlot = subplot(414);
plot(pitch, 'b-');
hold on;
plot(roll, 'r-');
plot(heading, 'g-');
hold off
title('Pitch Roll and Heading');
ylabel('Degrees');

linkaxes([accelPlot gyroPlot prhPlot], 'x');

%% Diagnostic depth and temperature plots

figure; 

depthPlot = subplot(211);
plot(depth, 'b-');
title('Depth');
ylabel('Meters');

tempPlot = subplot(212);
plot(temperature, 'r-');
title('Temperature');
ylabel('Deg Celsius');

linkaxes([depthPlot tempPlot], 'x');

%% Diagnostic power measurement plots

figure;

voltagePlot = subplot(511);
plot(voltage, 'b-');
title('Voltage');
ylabel('mV');

currentPlot = subplot(512);
plot(current, 'g-');
title('Current');
ylabel('mA');

powerPlot = subplot(513);
plot(power, 'r-');
title('Power Consumption');
ylabel('mW');

chargePlot = subplot(514);
plot(chargeState, 'm-');
title('State of Charge');
ylabel('%');

capacityPlot = subplot(515);
plot(remainingCapacity, 'k-');
title('Remaining Capacity');
ylabel('mAh');

linkaxes([voltagePlot currentPlot powerPlot chargePlot capacityPlot], 'x');

%% Decimate LED1-4 (250 Hz) & gyroscope data (100 Hz) to 25 Hz

led1_25Hz = decimate(led1, 10);         % decimate 250 Hz by a factor of 10
led2_25Hz = decimate(led2, 10);
led3_25Hz = decimate(led3, 10);
led4_25Hz = decimate(led4, 10);

ax_25Hz = decimate(ax, 4);              % decimate 100 Hz by a factor of 4
ay_25Hz = decimate(ay, 4);
az_25Hz = decimate(az, 4);

gx_25Hz = decimate(gx, 4);
gy_25Hz = decimate(gy, 4);
gz_25Hz = decimate(gz, 4);

pitch_25Hz   = decimate(pitch, 4);
roll_25Hz    = decimate(roll, 4);
heading_25Hz = decimate(heading, 4);

depth_25Hz          = decimate(depth, 4);
temperature_25Hz    = decimate(temperature, 4);

%% try looking at some normalized data for signal:noise comparisons...

led1_d_norm = normalize(led1_25Hz);
led2_d_norm = normalize(led2_25Hz);
led3_d_norm = normalize(led3_25Hz);
led4_d_norm = normalize(led4_25Hz);

figure;
plot(led1_d_norm, 'k.');
hold on;
plot(led2_d_norm, 'r-');
plot(led3_d_norm, 'g-');
plot(led4_d_norm, 'b-');
hold off;
title('Overlaid normalized optical signals');

%% create an optics vector and compute a mean of all sampled wavelengths

opticsVector = [led1_d_norm; led2_d_norm; led3_d_norm; led4_d_norm];
opticalMean  = mean(opticsVector);

figure;
plot(opticalMean, 'b-');
title('Mean optical signal, all wavelengths');
grid;

%% normalize decimated accel & gyro data

ax_d_norm   = normalize(ax_25Hz);
ay_d_norm   = normalize(ay_25Hz);
az_d_norm   = normalize(az_25Hz);

gx_d_norm   = normalize(gx_25Hz);
gy_d_norm   = normalize(gy_25Hz);
gz_d_norm   = normalize(gz_25Hz);

figure;
plot(ax_d_norm, 'b.');
hold on;
plot(ay_d_norm, 'r.');
plot(az_d_norm, 'g.');
hold off;
title('Norm accel');
grid;

figure;
plot(gx_d_norm, 'b.');
hold on;
plot(gy_d_norm, 'r.');
plot(gz_d_norm, 'g.');
hold off;
title('Norm gyro');
grid;

%% create an accel and gyro vector and compute the mean using all axes

accelVector = [ax_d_norm; ay_d_norm; az_d_norm];
gyroVector  = [gx_d_norm; gy_d_norm; gz_d_norm];

accelMean   = mean(accelVector);
gyroMean    = mean(gyroVector);

figure;
subplot(211);
plot(accelMean, 'r-');
grid;
subplot(212);
plot(gyroMean, 'm-');
grid;


%% put it all together...

%  plot overalls of norm'd optical mean, with norm'd accel & gyro mean

figure;
plot(opticalMean, 'b-', 'LineWidth', 2);
hold on;
plot(accelMean, 'r-');
plot(gyroMean, 'm-');
hold off;
grid;
title('Mean normalized decimated optical, accelerometer, and gyroscope signal');

%% AC/DC and FFT with decimated optical signals

% use decimated LED1 as signal length
sigLength = length(led1_25Hz);

% deprecate Ts, since we always know Fs
% Ts     = mean(diff(time_s));          % average time sample

% declare Fs explicitly, since we always decimate to 25 Hz here
% Fs     = 1 / Ts;                      % effective sampling rate
Fs = 25;

% now try some FFT magic...
nFFT   = 2 ^ nextpow2(sigLength);     % next power of 2 from sig length
ambFFT = fft(led1_25Hz, nFFT) / sigLength;      % LED1 FFT
led2FFT = fft(led2_25Hz, nFFT) / sigLength;     % LED2 FFT
led3FFT = fft(led3_25Hz, nFFT) / sigLength;     % LED3 FFT
led4FFT = fft(led4_25Hz, nFFT) / sigLength;     % LED4 FFT

f      = Fs/2 * linspace(0,1,nFFT/2+1);
freq   = -Fs/2 + Fs/nFFT:Fs/nFFT:Fs/2;


ambCenter  = fftshift(ambFFT);
led2Center = fftshift(led2FFT);
led3Center = fftshift(led3FFT);
led4Center = fftshift(led4FFT);

ambAC      = fftshift(ambFFT(2:end));
led2AC     = fftshift(led2FFT(2:end));
led3AC     = fftshift(led3FFT(2:end));
led4AC     = fftshift(led4FFT(2:end));

ambDC      = ambFFT(1);
led2DC     = led2FFT(1);
led3DC     = led3FFT(1);
led4DC     = led4FFT(1);

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
cutoffValue       = 2000;  % was 4000, but trying with 3000 = 1 BPM          
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

%% Build an x-scale in seconds based on samples and decimated frequency

numSamples = length(led1_25Hz);
samplingFreqHz = Fs;
samplingIncrement = 1/samplingFreqHz;

time_s = 0:samplingIncrement:(numSamples/samplingFreqHz)-samplingIncrement;
time = time_s';

%% Plot overlaid filtered ppgWaveforms 

figure;

plot(time_s, filtered_led1_25Hz,'b-');
xlabel('Time, s');
ylabel('LED1-4 Reflectance');
tempTitle = sprintf('Filtered G/R/NIR PPG (%d° BW w %0.3g Hz Cutoff)', ...
    butterOrder, cutoffFrequency);
title(tempTitle);
hold on;
plot(time_s, filtered_led2_25Hz, 'b-');
plot(time_s, filtered_led3_25Hz, 'r-');
plot(time_s, filtered_led4_25Hz, 'g-');
hold off;
grid;

%% Plot the filtered ppgWaveforms in individual panels

figure;

fig3a = subplot(411);
plot(time_s, filtered_led1_25Hz, 'k-');
ylabel('Ambient Reflectance');
% ylim([-200 200]);
tempTitle = sprintf('Filtered ambient ((%d° BW w %0.3g Hz Cutoff)', ...
    butterOrder, cutoffFrequency);
title(tempTitle);
grid;

fig3b = subplot(412);
plot(time_s, filtered_led2_25Hz,'b-');
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

%% Make some spectral analysis plots

welchNFFT     = 2^7;        % try segment lengths of 2^8 = 256
winSize  = hanning(welchNFFT);   % set hanning window shape
nOverlap = welchNFFT / 2;        % set 50% overlap between segments

[Pled2, Fled2] = periodogram(led2_25Hz, [], welchNFFT, Fs, 'power');
figure; 
plot(Fled2,10*log10(Pled2),'b-');
xlabel('Frequency, Hz'); ylabel('Power spectrum (dBW) LED2 1050 nm');

[Pled2Power,Fled2Power] = pwelch(led2_25Hz,ones(welchNFFT,1),0,welchNFFT,Fs,'power');
figure; 
plot(Fled2Power,10*log10(Pled2Power),'b-');
xlabel('Frequency, Hz'); ylabel('Power spectrum (dBW)');

[psd_led1, f_led1] = pwelch(filtered_led1_25Hz, winSize, nOverlap, welchNFFT);
[psd_led2, f_led2] = pwelch(filtered_led2_25Hz, winSize, nOverlap, welchNFFT);
[psd_led3, f_led3] = pwelch(filtered_led3_25Hz, winSize, nOverlap, welchNFFT);
[psd_led4, f_led4] = pwelch(filtered_led4_25Hz, winSize, nOverlap, welchNFFT);


%% Try to make some spectrograms and periodograms of pulse ox data

[S, F, T, P] = spectrogram(led2_25Hz,winSize,nOverlap,nFFT,Fs);

% make surface3d plot
figure;
colormap(parula(5));
% surface with edgecolor
% surfc(T,F,10*log10(abs(P)));
% surface with no edgecolor
surfc(T,F,10*log10(abs(P)), 'EdgeColor','none');
xlabel('Time,s'); ylabel('Frequency, Hz'); zlabel('Magnitude, dB');
axis tight;
view(-45,45);

%% make contour3d plot

figure;
colormap(parula(5));
% surfc(T,F,10*log10(abs(P)), 'EdgeColor','none');
contour3(T,F,10*log10(abs(P)));
xlabel('Time,s'); ylabel('Frequency, Hz'); 
c = colorbar;
c.Label.String = 'Intensity, dB';
axis tight;

%% Try a static spectrogram plot

colormap(parula(5));
colorLimit = [40 90] ;  % color axis limits in dB for specgram
figure; 
imagesc(T, F, 10*log10(abs(P)));
% contourf(10*log10(abs(P)));
xlabel('Time,s'); ylabel('Frequency, Hz'); 
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

    figPSD3 = subplot(413);
    loglog(f_led3, psd_led3, '-', 'Color', 'red');
    xlabel ('Frequency (Hz)');
    ylabel ('PSD');
    hold on;
    loglog( [-log10(cutoffFrequency) -log10(cutoffFrequency)], [10^2 10^8], 'k:');
    hold off;
    tempTitle = sprintf('Filtered LED3 PPG PSD');
    title(tempTitle);
    
    figPSD4 = subplot(414);
    loglog(f_led4, psd_led4, '-', 'Color', 'green');
    xlabel ('Frequency (Hz)');
    ylabel ('PSD');
    hold on;
    loglog( [-log10(cutoffFrequency) -log10(cutoffFrequency)], [10^2 10^8], 'k:');
    hold off;    
    tempTitle = sprintf('Filtered LED4 PPG PSD');
    title(tempTitle);
    
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
ylabel('LED1 Intensity');
% ylim([-200, 200]);
tempTitle = sprintf('Filtered LED1 PPG w/ 1° and 2nd Derivatives');
title(tempTitle);

fig4b = subplot(312);
plot(time_s, d1led1, 'k-');
ylabel('1° LED1 Intensity');

fig4c = subplot(313);
plot(time_s, d2led1,'k-');
xlabel('Time, s');
ylabel('2° Red Intensity');

linkaxes( ([fig4a, fig4b, fig4c]), 'x');

%% plot 1st and 2nd ppg derivatives for led2

figure; 

fig4a = subplot(311);

plot(time_s, filtered_led2_25Hz, 'b-');
% [maxtab, mintab] = peakdet(filteredPPGred, 0.5);
%hold on; plot(mintab(:,1), mintab(:,2), 'g*');
%plot(maxtab(:,1), maxtab(:,2), 'r*');
ylabel('LED2 Intensity');
% ylim([-200, 200]);
tempTitle = sprintf('Filtered LED2 PPG w/ 1° and 2nd Derivatives');
title(tempTitle);

fig4b = subplot(312);
plot(time_s, d1led2, 'b-');
ylabel('1° LED2 Intensity');

fig4c = subplot(313);
plot(time_s, d2led2,'b-');
xlabel('Time, s');
ylabel('2° Red Intensity');

linkaxes( ([fig4a, fig4b, fig4c]), 'x');

%% plot 1st and 2nd ppg derivatives for led3

figure; 

fig4a = subplot(311);

plot(time_s, filtered_led3_25Hz, 'r-');
% [maxtab, mintab] = peakdet(filteredPPGred, 0.5);
%hold on; plot(mintab(:,1), mintab(:,2), 'g*');
%plot(maxtab(:,1), maxtab(:,2), 'r*');
ylabel('LED3 Intensity');
% ylim([-200, 200]);
tempTitle = sprintf('Filtered LED3 PPG w/ 1° and 2nd Derivatives');
title(tempTitle);

fig4b = subplot(312);
plot(time_s, d1led3, 'r-');
ylabel('1° LED3 Intensity');

fig4c = subplot(313);
plot(time_s, d2led3,'r-');
xlabel('Time, s');
ylabel('2° LED3 Intensity');

linkaxes( ([fig4a, fig4b, fig4c]), 'x');

%% plot 1st and 2nd ppg derivatives for led4

figure; 

fig4a = subplot(311);

plot(time_s, filtered_led4_25Hz, 'g-');
% [maxtab, mintab] = peakdet(filteredPPGred, 0.5);
%hold on; plot(mintab(:,1), mintab(:,2), 'g*');
%plot(maxtab(:,1), maxtab(:,2), 'r*');
ylabel('LED4 Intensity');
% ylim([-200, 200]);
tempTitle = sprintf('Filtered LED4 PPG w/ 1° and 2nd Derivatives');
title(tempTitle);

fig4b = subplot(312);
plot(time_s, d1led4, 'g-');
ylabel('1° LED4 Intensity');

fig4c = subplot(313);
plot(time_s, d2led4,'g-');
xlabel('Time, s');
ylabel('2° LED4 Intensity');

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
fig9b = subplot(212);
plot(F, phaseLed2FFT, 'b-');
xlabel('Frequency, Hz');
ylabel('Phase, radians');
tempTitle = sprintf('LED2 Phase FFT');
title(tempTitle);

%% do some wavelet stuff
   
doCWT = 1;

if doCWT == 1

    fprintf('\tWorking on cwt...\n');
    
    % do cwt for raw led2 and led3
    figure;
    subplot(211);
    cwt(led2_25Hz,Fs);
    caxis([0 200]);
    title('cwt(led2)');
    subplot(212);
    caxis([0 500]);
    cwt(led3_25Hz,Fs);
    title('cwt(nir)');

    % do cwt for filtered red and nir
    figure;
    subplot(211);
    cwt(filtered_led2_25Hz,Fs);
    caxis([0 200]);
    title('cwt(filteredLed2)');
    subplot(212);
    cwt(filtered_led3_25Hz,Fs);
    caxis([0 500]);
    title('cwt(filteredLed3)');

end

% do wsst for for raw led2 and led3

doWSST = 1;

if doWSST == 1

    fprintf('\tWorking on wsst...\n');

    % do wsst for for raw led2 and led3

    figure;
    subplot(211);
    wsst(led2_25Hz,Fs);
    ylim([0 5]);
    caxis([0 10]);
    title('wsst(led2)');
    subplot(212);
    wsst(led3_25Hz,Fs);
    ylim([0 5]);
    caxis([0 10]);
    title('wsst(led3)');

    % do wsst for for filtered led2 and led3

    figure;
    figWSSTred = subplot(211);
    wsst(filtered_led2_25Hz,Fs);
    ylim([0 5]);
    caxis([0 5]);
    title('wsst(filteredLed2)');
    figWSSTnir = subplot(212);
    wsst(filtered_led3_25Hz,Fs);
    ylim([0 5]);
    caxis([0 10]);
    title('wsst(filteredLed3)');

    % do some manual wsst stuff for plotting contours of least
    % penalized WSST results

    doManualWSST = 1;

    if doManualWSST == 1

        % build wsst for signals of interest

        [wsstLed2, wsstLed2F]         = wsst(led2_25Hz,Fs);
        [wsstLed3, wsstLed3F]         = wsst(led3_25Hz,Fs);
        
        [wsstFiltLed2, wsstFiltLed2F] = wsst(filtered_led2_25Hz,Fs);
        [wsstFiltLed3, wsstFiltLed3F] = wsst(filtered_led2_25Hz,Fs);

        % build wsst ridges for signals of interest, using 'NumRidges'
        % set to 2, to get primary signal and 1st harmonic
        % and run a loop (for now) to explore the effect of 'penalty'

        doMultiPenaltyTest = 1;

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

            figure;
            contour(time_s,wsstLed2F,abs(wsstLed2));
            hold on;
            plot(time_s, wsstLed2Ridge, 'b.');
            hold off;
            xlabel('Time, s'); ylabel('Frequency, Hz');
            grid on;
            ylim([0 4]);
            caxis([0 5]);
            title('Synchrosqueezed Transform - Raw LED2 PPG');

            figure;
            contour(time_s,wsstFiltLed2F,abs(wsstFiltLed2));
            hold on;
            plot(time_s, wsstFiltLed2Ridge, 'b.');
            hold off;
            xlabel('Time, s'); ylabel('Frequency, Hz');
            grid on;
            ylim([0 4]);
            caxis([0 5]);
            title('Synchrosqueezed Transform for Filtered LED2 PPG');            

            figure;
            contour(time_s,wsstLed3F,abs(wsstLed3));
            hold on;
            plot(time_s, wsstLed3Ridge, 'r.');
            hold off;
            xlabel('Time, s'); ylabel('Frequency, Hz');
            grid on;
            ylim([0 4]);
            caxis([0 10]);
            title('Synchrosqueezed Transform for Raw LED3 PPG');  

            figure;
            contour(time_s,wsstFiltLed3F,abs(wsstFiltLed3));
            hold on;
            plot(time_s, wsstFiltLed3Ridge, 'r.');
            hold off;
            xlabel('Time, s'); ylabel('Frequency, Hz');
            grid on;
            ylim([0 4]);
            caxis([0 10]);
            title('Synchrosqueezed Transform for Filtered LED3 PPG'); 

            % plot 1st wsstridge for filtered red and nir as -
            % plot 2nd wsstridge for filtered red and nir as -.

            figure;
            plot(time_s, wsstFiltLed2Ridge(:,1), 'b-');
            hold on; 
            plot(time_s, wsstFiltLed2Ridge(:,2), 'b--');
            plot(time_s, wsstFiltLed3Ridge(:,1), 'r-');
            plot(time_s, wsstFiltLed3Ridge(:,2), 'r--');
            hold off;
            xlabel('Time, s'); ylabel('Frequency, Hz');
            ylim([0 2.5]);
            title('Synchrosqueezed Ridges - Raw and Filtered PPGs');
            legend('LED2 Pulse','LED3 Pulse','LED2 1st Harmonic', ...
                'LED3 1st Harmonic');

            figure;
            plot(time_s, led2PulseBPM,'b-');
            hold on;
            plot(time_s, led3PulseBPM,'r.');
            hold off;
            xlabel('Time, s'); ylabel('Estimated Heart Rate, BPM');
            title('Synchrosqueeze Ridges - Estimated Heart Rate');
            grid on;

            figure;
            plot(time_s,led2PulseBPMh1,'b-');
            hold on;
            plot(time_s,led3PulseBPMh1,'r.');
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
                contour(time_s,wsstLed2F,abs(wsstLed2));
                hold on;
                plot(time_s, wsstLed2Ridge, 'b.');
                hold off;
                xlabel('Time, s'); ylabel('Frequency, Hz');
                grid on;
                ylim([0 4]);
                caxis([0 5]);
                title('Synchrosqueezed Transform - Raw LED2 PPG');

                figure;
                contour(time_s,wsstFiltLed2F,abs(wsstFiltLed2));
                hold on;
                plot(time_s, wsstFiltLed2Ridge, 'b.');
                hold off;
                xlabel('Time, s'); ylabel('Frequency, Hz');
                grid on;
                ylim([0 4]);
                caxis([0 5]);
                title('Synchrosqueezed Transform for Filtered LED2 PPG');            

                figure;
                contour(time_s,wsstLed3F,abs(wsstLed3));
                hold on;
                plot(time_s, wsstLed3Ridge, 'r.');
                hold off;
                xlabel('Time, s'); ylabel('Frequency, Hz');
                grid on;
                ylim([0 4]);
                caxis([0 10]);
                title('Synchrosqueezed Transform for Raw LED3 PPG');  

                figure;
                contour(time_s,wsstFiltLed3F,abs(wsstFiltLed3));
                hold on;
                plot(time_s, wsstFiltLed3Ridge, 'r.');
                hold off;
                xlabel('Time, s'); ylabel('Frequency, Hz');
                grid on;
                ylim([0 4]);
                caxis([0 10]);
                title('Synchrosqueezed Transform for Filtered LED3 PPG'); 

                % plot 1st wsstridge for filtered LED2 and LED3 as -
                % plot 2nd wsstridge for filtered LED2 and LED3 as -.

                figure;
                plot(time_s, wsstFiltLed2Ridge(:,1), 'b-');
                hold on; 
                plot(time_s, wsstFiltLed2Ridge(:,2), 'b--');
                plot(time_s, wsstFiltLed3Ridge(:,1), 'r-');
                plot(time_s, wsstFiltLed3Ridge(:,2), 'r--');
                hold off;
                xlabel('Time, s'); ylabel('Frequency, Hz');
                ylim([0 2.5]);
                title('Synchrosqueezed Ridges - Raw and Filtered PPGs');
                legend('LED2 Pulse','LED3 Pulse','LED2 1st Harmonic', ...
                    'LED3 1st Harmonic');

                figure;
                plot(time_s,led2PulseBPM,'b-');
                hold on;
                plot(time_s,led3PulseBPM,'r--');
                hold off;
                xlabel('Time, s'); ylabel('Estimated Heart Rate, BPM');
                title('Synchrosqueeze Ridges - Estimated Heart Rate');
                grid on;

                figure;
                plot(time_s,led2PulseBPMh1,'b-');
                hold on;
                plot(time_s,led3PulseBPMh1,'r--');
                hold off;
                xlabel('Time, s'); ylabel('Estimated 1H Heart Rate, BPM');
                title('Synchrosqueezed Ridges - Estimated HR 1st Harmonic');
                grid on;
        end

        end

    end

end


%% Experiment with ambient light subtraction...

led2_minusAmb_25Hz = led2_25Hz - led1_25Hz;
led3_minusAmb_25Hz = led3_25Hz - led1_25Hz;
led4_minusAmb_25Hz = led4_25Hz - led1_25Hz;

figure;
led2lessAmb = subplot(311);
plot(time_s, led2_minusAmb_25Hz, 'b-');
grid;
led3lessAmb = subplot(312);
plot(time_s, led3_minusAmb_25Hz, 'r-');
grid;
led4lessAmb = subplot(313);
plot(time_s, led4_minusAmb_25Hz, 'g-');
hold on;
plot(time_s, accelMean, 'k--');
hold off;
grid;

linkaxes([led2lessAmb led3lessAmb led4lessAmb], 'x');

%% Experiment with subtracting accelMean from normalized ledX_minusAmb_25Hz

led2_minusAmb_25Hz_norm = normalize(led2_minusAmb_25Hz);
led3_minusAmb_25Hz_norm = normalize(led3_minusAmb_25Hz);
led4_minusAmb_25Hz_norm = normalize(led4_minusAmb_25Hz);

led2_noAmb_lessAccel = led2_minusAmb_25Hz_norm - accelMean;
led3_noAmb_lessAccel = led3_minusAmb_25Hz_norm - accelMean;
led4_noAmb_lessAccel = led4_minusAmb_25Hz_norm - accelMean;

figure;
led2lessAmbAccel = subplot(411);
plot(time_s, led2_noAmb_lessAccel, 'b-');
grid;
led3lessAmbAccel = subplot(412);
plot(time_s, led3_noAmb_lessAccel, 'r-');
grid;
led4lessAmbAccel = subplot(413);
plot(time_s, led4_noAmb_lessAccel, 'g-');
grid;
accelMeanPlot = subplot(414);
plot(time_s, accelMean, 'k-');
grid;

linkaxes([led2lessAmbAccel led3lessAmbAccel led4lessAmbAccel accelMeanPlot], 'x');

%% try dddtree cleanup on a segment of unfiltered red ppg

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


