%% ftKinetoCardiographyEnergetics.m
%
%   Written by Dave Haas, 15 August to 1 September 2021
%
%   This is used to compute kinetic energy, work and power from
%   kinetocardiographic (i.e.: accelerometer and gyroscope) signals in
%   Faunatag attachments.

clc;
clear;

%%  Declare some commonly used variables

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

% set some default LED* colors...

led1Color = black;
led2Color = blue;
led3Color = red;
led4Color = goldenrod;

% ... and some default kinematic colors...

kxColor = blue;
kyColor = red;
kzColor = goldenrod;

%% Make introductions

fprintf('Welcome to the FaunaTag cardiac energetics estimator.\n');


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
                confText = 'Epoch analysis results found! Continue? (y/n): ';
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
                        validEpoch = true;                
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
                        validEpoch = true;                
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
                        validEpoch = true;                
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
                        validEpoch = true;                
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

%% Do some regional selections with wsst and wsstridge finding

screenSize = get(0,'screensize');

figAnalysisRegion = figure('Color', white);
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

%% Make band-pass filtered variables for signal calculations 

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

brickWallHighPassFilter = designfilt('highpassiir','FilterOrder', 3, ...
    'PassbandFrequency', 50, 'PassbandRipple', 0.2, ...
    'SampleRate', 100);

bpAx   = filtfilt(bandPassHrKinematicFilter', KINEMATICS.ax(analysisRange));
bpAy   = filtfilt(bandPassHrKinematicFilter', KINEMATICS.ay(analysisRange));
bpAz   = filtfilt(bandPassHrKinematicFilter', KINEMATICS.az(analysisRange));
bpGx   = filtfilt(bandPassHrKinematicFilter', KINEMATICS.gx(analysisRange));
bpGy   = filtfilt(bandPassHrKinematicFilter', KINEMATICS.gy(analysisRange));
bpGz   = filtfilt(bandPassHrKinematicFilter', KINEMATICS.gz(analysisRange));

bpOdba = sqrt( bpAx.^2 + bpAy.^2 + bpAz.^2);
bpOdav = sqrt( bpGx.^2 + bpGy.^2 + bpGz.^2);

%% try this for the whole of bpA*

if ( mod(numel(bpAx),2) == 0 )      % even
    sampleRange = 1:numel(bpAx);
else                                % odd
    sampleRange = 1:numel(bpAx)-1;
end

% assemble Axyz using band-pass filtered accel

Axyz = [bpAx bpAy bpAz];

Gxyz = [bpGx bpGy bpGz];

% use size(Axyz) to make a container for Vxyz

sizeVxyz = size(Axyz);

linearVxyz(sizeVxyz(1),sizeVxyz(2)) = zeros;


% compute instantaneous velocities using Simpson's rule

for i = 1:numel(sampleRange)-1

    linearVxyz(i,:) = simps(Axyz(i:i+1,:));
    
end

figure; 
plot(KINEMATICS.Time(analysisRange), linearVxyz); grid;

%% try to create kinematic energy for this time series

massTag         = 0.403;    % mass (kg) of FaunaTag 111 = 403 grams
radiusTagCenter = 0.085;    % r from suction cup edge to sensor centerpoint
radiusTagZ      = 0.0043;   % r from contact sensor to point at tag top

Ix = mamssTag * radiusTagCenter ^  2;
Iy = massTag * radiusTagCenter ^ 2;
Iz = massTag * radiusTagZ ^ 2;

keLin = 0.5 * massTag * ( linearVxyz(:,1).^2 + linearVxyz(:,2).^2 + linearVxyz(:,3).^2 );

keRot = 0.5 * ( (Ix * bpGx).^2 + (Iy * bpGy).^2 + (Iz * bpGz).^2);

simpsonLinearKE((length(keLin)-1)) = zeros;

for i = 1:length(linearVxyz) - 1

%     tLinKE(i,:) = ...
%         trapz(0.5 * 250 * ( linearVxyz(i:i+1,:).^2 ) );
        
    iLinKE(i,:) = simps(0.5 * massTag * ( linearVxyz(i:i+1,:).^2 ) );
    iRotKE(i,:) = simps(0.5 * ( (Ix * Gxyz(i:i+1,:) ).^2 ));
    
end

figure('Color',white);
sKElin = subplot(211);
plot(KINEMATICS.time_s(analysisRange), keLin, 'Color', red);
xlabel('Time, seconds');
ylabel('J');
title('Linear kinetic energy (accelerometer)');
sKErot = subplot(212);
plot(KINEMATICS.time_s(analysisRange), keRot, 'Color', green);
xlabel('Time, seconds');
ylabel('J');
title('Rotation kinetic energy (gyroscope)');
linkaxes([sKElin,sKErot],'x');

thisT_s = KINEMATICS.time_s(analysisRange);
sampleStart = 1;
sampleEnd = length(thisT_s);
sampleRange = sampleStart:sampleEnd;

figure('Color',white);
pIklin = subplot(211);
plot(thisT_s(1:end-1), iLinKE);
xlabel('Time, seconds');
ylabel('J·s');
title('iK(lin)');
grid;
legend('ax','ay','az');
pIkrot = subplot(212);
plot(thisT_s(1:end-1), iRotKE);
xlabel('Time, seconds');
ylabel('J·s');
title('iK(rot)');
grid;
legend('gx','gy','gz');
linkaxes([pIklin, pIkrot],'x')

%% 

rmsILinKE = rms(iLinKE);
rmsiRotKE = rms(iRotKE);


lenKE = length(iRotKE);
lenTime = length(KINEMATICS.Time);

tOffset = lenTime - lenKE;

figure('Color', white);

pKEt = subplot(511);
plot(KINEMATICS.Time(1:end-tOffset+1), keLin);
xlabel('Time, seconds');
ylabel('mJ');
grid;
legend('KE_Axyz');

pKEs = subplot(512);
plot(KINEMATICS.Time(1:end-tOffset), iLinKE);
xlabel('Time, seconds');
ylabel('mJ');
grid;
legend('ax','ay','az');

pKEx = subplot(513);
plot(KINEMATICS.Time(1:end-tOffset), iLinKE(:,1), ...
    'Color', blue);
xlabel('Time, seconds');
ylabel('mJ');
title('Cardiac energy: x-axis SCG');
grid;
legend('ax');

pKEy = subplot(514);
plot(KINEMATICS.Time(1:end-tOffset), iLinKE(:,2), ...
    'Color', red);
xlabel('Time, seconds');
ylabel('mJ');
title('Cardiac energy: y-axis SCG');
grid;
legend('ay');

pKEz = subplot(515);
plot(KINEMATICS.Time(1:end-tOffset), iLinKE(:,3), ...
    'Color', goldenrod);
xlabel('Time, seconds');
ylabel('mJ');
title('Cardiac energy: z-axis SCG');
grid;
legend('az');

linkaxes([pKEt pKEs pKEx pKEy pKEz], 'x');