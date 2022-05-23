%% ftEcgFilter.m
%
%       Created and written by Dave Haas, 16-17 May 2022

clc;
clear;

%%  Load global TAG_PATHS and verify they're where they should be...

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

%% define some standard FaunaLabs plotting colors

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

%%

ecgLowCut           = 0.5;
ecgHighCut          = 2;
powerLineFreq       = 60;
ecgBpFilterOrder    = 5;
ecgFs               = opticsFs / 10;


[ecgAbp, ecgBbp, ecgCbp, ecgDbp] = butter(ecgBpFilterOrder, ...
    [ecgLowCut ecgHighCut] / (ecgFs) );

[sosBpEcg,gBpEcg] = ss2sos(ecgAbp,ecgBbp,ecgCbp,ecgDbp);

if (mod(ecgBpFilterOrder,2) == 0)
    ecgDesignFiltOrder = ecgBpFilterOrder;
else
    ecgDesignFiltOrder = ecgBpFilterOrder + 1;
end

ecgBandPassFilter = designfilt('bandpassiir', ...
    'FilterOrder', ecgDesignFiltOrder, ...
    'HalfPowerFrequency1', ecgLowCut, ...
    'HalfPowerFrequency2', ecgHighCut, ...
	'SampleRate', ecgFs);

makePlots = true;

if (makePlots)
    fvtO = fvtool(sosBpEcg, ecgBandPassFilter, 'Fs', ecgFs);
    legend(fvtO,'butter','designfilt');
end

ECG = struct;
ECG.BANDPASS_FILTER.filterOrder = ecgBpFilterOrder;
ECG.BANDPASS_FILTER.lowCut = ecgLowCut;
ECG.BANDPASS_FILTER.highCut = ecgHighCut;
ECG.BANDPASS_FILTER.sampleRate = ecgFs;
ECG.BANDPASS_FILTER.A = ecgAbp;
ECG.BANDPASS_FILTER.B = ecgBbp;
ECG.BANDPASS_FILTER.C = ecgCbp;
ECG.BANDPASS_FILTER.D = ecgDbp;
ECG.BANDPASS_FILTER.SOS = sosBpEcg;
ECG.BANDPASS_FILTER.G = gBpEcg;

fprintf('Band-pass butterworth filter for ECG signal:\n');
fprintf('\tFilter order: %d\n', ecgBpFilterOrder);
fprintf('\tLow cut frequency: %2.2f Hz\n', ecgLowCut);
fprintf('\tHigh cut frequency: %2.2f Hz\n', ecgHighCut);
fprintf('\tECG BPF sample rate (2Hz = normalized frequencies): %d Hz\n', ...
    ecgFs);

filtEcg = filtfilt( sosBpEcg, gBpEcg, optics.led4);


%% Low-pass experiment for sections of ECG without RLD attached

%% Next up, a high-pass filter for DC offset removal

dcoFiltOrder        = 3;                % FieldTrip suggested default
dcoPassbandFreqHz   = 0.01;             % DC-removal (cut below 0.01 Hz)
ecgNyquistFs        = ecgFs / 2;

fprintf('DC-removal (high-pass) filter characteristics:\n');
fprintf('\tFilter order: %d\n', dcoFiltOrder);
fprintf('\tPassband frequency: %d Hz\n', dcoPassbandFreqHz);
fprintf('\tNyquist frequency (ECG): %d\n', ecgNyquistFs);


% filter coefficients for optics high-pass (DC-removal)
[AecgDCHP, BecgDCHP, CecgDCHP, DecgDCHP] = ...
    butter(dcoFiltOrder, ( dcoPassbandFreqHz / (ecgNyquistFs) ), 'high');

[ecgSOS, ecgG] = ss2sos(AecgDCHP, BecgDCHP, CecgDCHP, DecgDCHP);


% design the high-pass filter. NOTE: using samplingRate = 2, which is the
% designfilt default, results in normalized frequencies

dcoEcgFilter = designfilt('highpassiir','FilterOrder',dcoFiltOrder, ...
    'PassbandFrequency', dcoPassbandFreqHz, ... 
    'PassbandRipple', 0.2, ...
    'SampleRate', 2);

if (makePlots)
    
    fvt = fvtool(ecgSOS, dcoEcgFilter, 'Fs', ecgFs);
    legend(fvt,'butter','designfilt');
    
    
end

dcoEcg = filtfilt(dcoEcgFilter', optics.led4);

%% Use low-pass filter to screen out 50/60 Hz noise

lpFiltOrder       = 6;        
ecgPassbandFreqHz = ecgFs / 25;

fprintf('Low-pass filter for removal of high frequency noise signals.\n');
fprintf('\tFilter order: %d\n', lpFiltOrder);
fprintf('\tECG passband frequency: %d Hz\n', ecgPassbandFreqHz);


% filter coefficients for high-pass

[Aolp, Bolp, Colp, Dolp] = butter(lpFiltOrder, ( ecgPassbandFreqHz / (ecgNyquistFs) ) );

% design the high-pass filter for optics
lowPassEcgFilter = designfilt('lowpassiir','FilterOrder', lpFiltOrder, ...
    'PassbandFrequency', ecgPassbandFreqHz, ... 
    'PassbandRipple', 0.2, ...
    'SampleRate', ecgNyquistFs);


lpEcg = filtfilt(lowPassEcgFilter', dcoEcg);

%% Use this filter to make raw and filtered ECG signal plots

f1 = figure('Color', white);

s1 = subplot(311);
plot(optics.Time, optics.led4, 'Color', goldenrod, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Raw ECG');
grid;
s2 = subplot(312);
plot(optics.Time, filtEcg, 'Color', red, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Filtered ECG');
title('Band-pass filter');
grid;
s3 = subplot(313);
plot(optics.Time, dcoEcg, 'Color', red, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Filtered ECG');
title('DCO filter');
grid;
linkaxes([s1 s2 s3],'x');

%% Use this filter to make raw and filtered ECG signal plots

f2 = figure('Color', white);

f2a = subplot(211);

plot(optics.Time, optics.led4, 'Color', goldenrod, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Raw ECG');
grid;

f2b = subplot(212);
plot(optics.Time, optics.led3, 'Color', red, 'LineWidth', 1.5);
xlabel('Time');
ylabel('A.U.');
title('Reflectance, 1200 nm');
grid;

linkaxes([f2a f2b],'x');
