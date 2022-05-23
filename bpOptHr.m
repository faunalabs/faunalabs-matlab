function [bpOptics, BPSTRUCT ] = optHrBp( opticalData, decimate )

bandPassFiltOrder  = 3;      % 6th order butterworth filter
oLowCut            = 0.05;    % ~30 beats per minute
% oHighCut           = 3.5;   % ~210 beats per minute <- 3.5 test 12.5
oHighCut           = 5.5;   % ~210 beats per minute <- 3.5 test 12.5
opticalTestFs      = oHighCut * 2;     % 8 produces good pulses on optics apnea        
    
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

bpOptics = filtfilt(sosBpo, gBpo, opticalData);
BPSTRUCT = OPTICS_BANDPASS_FILTER;

