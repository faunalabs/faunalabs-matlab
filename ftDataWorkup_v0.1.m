%% FaunaTag Data Workup - 
%       created by Dave Haas, 17 June 2020
%       updated by Dave Haas, 18 June 2021
%       
%       FaunaTag Data Workup fuses FaunaTagDataImporter.m with
%       pulseOx/ppgPlayground.m (also written by Dave Haas in 2017).
%       This synthesis of this work will direcly produce heart rate
%       estimates from bio-optical data collected by the FaunaTag
%
%       Dependency: csvproc from dtagtools by Mark Johnson

clear;

% do the work...
fprintf('\n\n---=== Starting FaunaTagDataWorkup! ===---\n\n');

useGoogle = false;

if (useGoogle)
    pathName = '/Volumes/GoogleDrive/Shared drives/FaunaData';
else
    pathName = '/Users/dave/Documents/FaunaData';
end

% specify subpath depending on whether its my data or animal data

animalData = true;                 %       <-- BE SURE TO SPECIFY THIS

if (animalData)
    subPath = 'DQO';
else
    subPath = 'DH';
end

tagID = 'FaunaTag111';              %       <-- BE SURE TO SPECIFY THIS

dateID = '20210514';                %       <-- BE SURE TO SPECIFY THIS

if(animalData)

    % FaunaTag 110 - all NIR LEDs failed after leaving NC :(
    
    % ----- ----- ----- ----- 3 May 2021 ----- ----- ----- ----- 
    %fileName = 'ft110_20210503_091613-10.txt';     % Kolohe III R
    %totd = 'X';
    %fileName = 'ft110_20210503_121325-10.txt';     % Kolohe II 
    %totd = 'X';
    %fileName = 'ft110_20210503_143350-10.txt';     % Hua I
    %totd = 'X';
    %fileName = 'ft110_20210503_153936-10.txt';     % Liho II
    %totd = 'X';

    % FaunaTag 111 
    
    % ----- ----- ----- ----- 3 May 2021 ----- ----- ----- ----- 
    %fileName = 'ft111_20210503_113626-10.txt';     % Hoku III R
    %totd = 'a';
    %fileName = 'ft111_20210503_1428250-10.txt';    % Hua II
    %totd = 'b';
    %fileName = 'ft111_20210503_153233-10.txt';     % Liho II
    %totd = 'c';
    % ----- ----- ----- ----- 7 May 2021 ----- ----- ----- ----- 
    %fileName = 'ft111_20210507_140742-10.txt';     % Kolohe I II & III L
    %totd = 'abc';
    %fileName = 'ft111_20210507_145342-10.txt';     % Hoku I & II
    %totd = 'de';
    % ----- ----- ----- ----- 8 May 2021 ----- ----- ----- ----- 
    %fileName = 'ft111_20210508_091836-10.txt';     % Hua I & II
    %totd = 'ab';
    %fileName = 'ft111_20210508_114531-10.txt';     % Hoku I II III L & IV  
    %totd = 'cdef';
    %fileName = 'ft111_20210508_143600-10.txt';     % Kolohe I & II
    %totd = 'gh';
    %fileName = 'ft111_20210508_155047-10.txt';     % Liho I
    %totd = 'i';
    % ----- ----- ----- ----- 10 May 2021 ----- ----- ----- ----- 
    %fileName = 'ft111_20210510_134417-10.txt';     % Kolohe III L
    %totd = 'a';
    %fileName = 'ft111_20210510_141230-10.txt';     % Hua I II III L
    %totd = 'bcd';
    %fileName = 'ft111_20210510_151032-10.txt';     % Noa I & II
    %totd = 'ef';
    %fileName = 'ft111_20210510_161759-10.txt';     % Liho II
    %totd = 'g';
    % ----- ----- ----- ----- 11 May 2021 ----- ----- ----- ----- 
    %fileName = 'ft111_20210511_102145-10.txt';     % Noa III L, no recover
    %totd = 'a';
    %fileName = 'ft111_20210511_113352-10.txt';     % Liho III R
    %totd = 'b';
    %fileName = 'ft111_20210511_132316-10.txt';     % DH & ambient laggon 1
    %totd = 'X';
    %fileName = 'ft111_20210511_141052-10.txt';     % Lono I
    %totd = 'c';
    %fileName = 'ft111_20210511_160612-10.txt';     % Lono II & III L
    %totd = 'de';
    
    
    % ----- ----- ----- -----             ----- ----- ----- ----- 
    % ----- ----- ----- vvvvv Foil Trials vvvvv ----- ----- ----- 
    % ----- ----- ----- -----             ----- ----- ----- ----- 
    
    
    % ----- ----- ----- ----- 13 May 2021 ----- ----- ----- ----- 
    
    %fileName = 'ft111_20210513_110326-10.txt';     % Hoku, foil I & II
    %totd = 'a';
    
    %fileName = 'ft111_20210513_114729-10.txt';     % Kolohe, foil I & II
    %totd = 'b';
    
    %fileName = 'ft111_20210513_140811-10.txt';     % Liho, foil I II III R
    %totd = 'c';
    
    %fileName = 'ft111_20210513_151511-10.txt';     % Kolohe w/ foil III R
    %totd = 'd';
    
    % ----- ----- ----- ----- 14 May 2021 ----- ----- ----- ----- 
    
    fileName = 'ft111_20210514_085503-10.txt';     % Liho w/ foil IV
    totd = 'a';
    
    %fileName = 'ft111_20210514_100519-10.txt';     % Kolohe w/ foil IV
    %totd = 'b';

    %fileName = 'ft111_20210514_121738-10.txt';     % Hoku w/ foil III
    %totd = 'c';

    % ----- ----- ----- ----- 21 May 2021 ----- ----- ----- ----- 
    
    %fileName = 'ft111_20210521_091400-10.txt';     % Hua w/ foil I & III
    %totd = 'a';
    
    %fileName = 'ft111_20210521_094909-10.txt';     % Lono w/ foil I & III
    %totd = 'b';
    
    %fileName = 'ft111_20210521_110615-10.txt';     % Hua w/ foil IV
    %totd = 'c';
    
    % ----- ----- ----- ----- 22 May 2021 ----- ----- ----- ----- 
    
    %fileName = 'ft111_20210522_090501-10.txt';     % Noa w/ foil I III & IV
    %totd = 'a';
    
    %fileName = 'ft111_20210522_112033-10.txt';     % Lono w/ foil IV
    %totd = 'b';
        
   
    % ----- ----- ----- ----- 3 May 2021 ----- ----- ----- ----- 
    
    % FaunaTag 112 - all NIR LEDs failed or there's a problem with the PD
    
    %fileName = 'ft112_20210503_142040-10.txt';     % Hua I
    %totd = 'a';
    
    % ^^^^^^^^^^^^^^^^ DOLPHIN DATA ^^^^^^^^^^^^^^^^^^^^
    
else        
    
    % vvvvvvvvvvvvvvvv HUMAN (DH) DATA vvvvvvvvvvvvvvvvv
    
    %fileName = 'ft106_20210506_110258-10.txt';
    %fileName = 'ft106_20210506_140815-10.txt';
    %fileName = 'ft110_20210505_170048-10.txt'; 
    %fileName = 'ft110_20210506_144129-10.txt';
    %fileName = 'ft110_20210506_145450-10.txt';
    %fileName = 'ft110_20210506_150322-10.txt';
    %fileName = 'ft110_20210506_150856-10.txt';
    %fileName = 'ft110_20210506_153231-10.txt';
    %fileName = 'ft110_20210506_190040-10.txt';
    %fileName = 'ft110_20210507_062207-10.txt';
    %fileName = 'ft110_20210507_062957-10.txt';
    %fileName = 'ft110_20210507_070110-10.txt';
    %fileName = 'ft110_20210507_115233-10.txt';
    %fileName = 'ft110_20210512_164025-10.txt';
    %fileName = 'ft110_20210512_165227-10.txt';
    
    %fileName = 'ft111_20210505_165240-10.txt'; 
    %fileName = 'ft111_20210506_161055-10.txt';
    %fileName = 'ft111_20210506_161747-10.txt';
    %fileName = 'ft111_20210506_164114-10.txt';
    %fileName = 'ft111_20210506_171603-10.txt';
    %fileName = 'ft111_20210506_174142-10.txt';
    %fileName = 'ft111_20210506_181726-10.txt'; 
    %fileName = 'ft111_20210507_073000-10.txt';
    %fileName = 'ft111_20210507_101715-10.txt';
    %fileName = 'ft111_20210507_103456-10.txt';
    %fileName = 'ft111_20210512_162301-10.txt';
    fileName = 'ft111_20210512_180612-10.txt';
    
    %fileName = 'ft112_20210507_071747-10.txt';
    %fileName = 'ft112_20210507_093534-10.txt';
    %fileName = 'ft112_20210507_122200-10.txt';
    
end


% now, build a complete targetFile for open...
targetFile = sprintf('%s/%s/%s/%s/%s', pathName, subPath, tagID, dateID, fileName);

fileExists = isfile(targetFile);

if (fileExists)
    fprintf('\tTarget file exists... reading data...\n');
    youShallNotPass = false;
    % load targetFile into S
    S = csvproc(targetFile, []);    
else
    fprintf('\tFile you specified does not seem to exist...\n\n');
    youShallNotPass = true;
    return;
end


%% parse fileName, make date, time, utcOffset, etc., generate prhFileName

fprintf('\tParsing filename into date-time components...\n');

if (strncmp(fileName, 'ft110', 5))
    FaunaTagID = 110;
elseif (strncmp(fileName, 'ft111', 5))
    FaunaTagID = 111;
elseif (strncmp(fileName, 'ft112', 5))
    FaunaTagID = 112;
end


YYYY = str2double(extractBetween(fileName, 7, 10));
YY   = str2double(extractBetween(fileName, 9, 10));
MM   = str2double(extractBetween(fileName, 11,12));
DD   = str2double(extractBetween(fileName, 13,14));

hh   = str2double(extractBetween(fileName, 16,17));
mm   = str2double(extractBetween(fileName, 18,19));
ss   = str2double(extractBetween(fileName, 20,21));

utc  = str2double(extractBetween(fileName, 22,24));

J    = jday([YYYY MM DD]);

trialDate = [YYYY MM DD];

trialDateTimeLocal = datetime(YYYY, MM, DD, hh, mm, ss);
trialDateTimeUTC   = datetime(YYYY, MM, DD, hh + (utc * -1), mm, ss);

if (YYYY == 2021 & MM == 5)  
    
    % these are all DQO dolphins, so auto-fill this info...
    
    gs = 'tt';
    lat = 21.272037;            % decimal degrees at Dolphin Quest Oahu
    long = -157.773116;         % decimal degrees at Dolphin Quest Oahu
    decl = 9.48;                % units: degrees, on 13 May 2021
    magFieldStrength = 34797.5; % field strength, units: nanoTeslas (nT)
else
    
    % ... but enter this manually if not DQO trials ...
    
    validMetaData = false;
    
    while(validMetaData == false)
        
        fprintf('\n----- Enter tag meta data information manually... -----\n\n');
    
        genus = input('Animal genus (e.g.: Tursiops): ','s');
        species = input('Animal species (e.g.: truncatus): ','s');
        gsCode = input('Enter the genusSpecies code for tagID (e.g.: tt): ', 's');
        latitude = input('Enter tagOn latitude in decimal degrees (S is negative): ', 's');
        longitude = input('Enter tagOn longitude in decimal degrees (W is negative): ', 's');
        declAngle = input('Enter declination angle in degrees: ', 's');
        fieldStrength = input('Enter magnetic field strength in nT: ', 's');

        lat = str2double(latitude);
        long = str2double(longitude);
        decl = str2double(declAngle);
        magFieldStrength = str2double(fieldStrength);

        % verify the condition of user entries...

        if ( isempty(genus) || isempty(species) || isempty(gsCode) || ...
                (strlength(gsCode) ~= 2) || ...
                isnan(lat)|| isnan(long) || ... 
                isnan(decl) || isnan(magFieldStrength ) )
            validMetaData = false;
            fprintf('\nTag meta data entered is invalid... Try re-entering it...\n\n');
        else
            validMetaData = true;
            fprintf('Tag meta data seems valid... continuing!\n\n');
        end
        
    end
    
end

rawFileName = sprintf('%s%d_%d%sraw.mat',gs,YY,J,totd);
prhFileName = sprintf('%s%d_%d%sprh.mat',gs,YY,J,totd);
tag = sprintf('%s%d_%d%s', gs,YY,J,totd);

fprintf('\tRAW filename: %s\n', rawFileName);
fprintf('\tPRH filename: %s\n', prhFileName);
fprintf('\tTag: %s\n', tag);

%% parse S into usable data components

if (youShallNotPass)
    fprintf('\tThe data file was missing. Find the right now and try again!\n');
    return;
end

fprintf('\tParsing data into usable components...\n');

% optics (relative insensity)           @ 250 Hz
led2 = str2double(S(1,:));
led3 = str2double(S(2,:));
led1 = str2double(S(3,:));              % this should be ambient
led4 = str2double(S(4,:));

% transform data for reflectance measurements, i.e.: make negative then
% rebase around lower bound values for rough DC offset

rLed1 = (led1 * -1);
rLed2 = (led2 * -1);
rLed3 = (led3 * -1);
rLed4 = (led4 * -1);

rLed1Min = min(rLed1);
rLed2Min = min(rLed2);
rLed3Min = min(rLed3);
rLed4Min = min(rLed4);

rLed1 = rLed1 - rLed1Min;
rLed2 = rLed2 - rLed2Min;
rLed3 = rLed3 - rLed3Min;
rLed4 = rLed4 - rLed4Min;

% normalize these signals for later use
normLed1 = normalize(led1);
normLed2 = normalize(led2);
normLed3 = normalize(led3);
normLed4 = normalize(led4);

nrLed1 = normalize(rLed1);
nrLed2 = normalize(rLed2);
nrLed3 = normalize(rLed3);
nrLed4 = normalize(rLed4);

% build the log10 optical signals
log10Led1 = log10(led1);
log10Led2 = log10(led2);
log10Led3 = log10(led3);
log10Led4 = log10(led4);
log10active = log10Led2 ./ log10Led3;
log10ambient = log10Led4 ./ log10Led1;

% bioOptics (optical intensity)         @ 250 Hz
num250HzSamples = length(led2);
numSecondsData = num250HzSamples / 250;

% accel (meters per second squared)     @ 100 Hz
ax = str2double(S(5,1:numSecondsData * 100));
ay = str2double(S(6,1:numSecondsData * 100));
az = str2double(S(7,1:numSecondsData * 100));

% use raw accel to create high-resolution ODBA and normalized ODBA
odba = sqrt(ax.^2 + ay.^2 + az.^2);
normOdba = normalize(sqrt(ax.^2 + ay.^2 + az.^2));

% gyro (degrees per second)             @ 100 Hz
gx = str2double(S(8,1:numSecondsData * 100));
gy = str2double(S(9,1:numSecondsData * 100));
gz = str2double(S(10,1:numSecondsData * 100));

% prh (degree)                          @ 100 Hz
pitch   = str2double(S(11,1:numSecondsData * 100));
roll    = str2double(S(12,1:numSecondsData * 100));
heading = str2double(S(13,1:numSecondsData * 100));

% depth (meters)                        @ 100 Hz
depth = str2double(S(14,1:numSecondsData * 100));

% temp (degrees Celsius)                @ 100 Hz
temperature = str2double(S(15,1:numSecondsData * 100));

% mag (degree)                                  @ 20 Hz
mx = str2double(S(16,1:numSecondsData * 20));
my = str2double(S(17,1:numSecondsData * 20));
mz = str2double(S(18,1:numSecondsData * 20));

% timestamp                                     @ 2 Hz
sampletime = str2double(S(19,1:numSecondsData * 2));

% power measurements                            @ 2 Hz
voltage             = str2double(S(20,1:numSecondsData * 2));
current             = str2double(S(21,1:numSecondsData * 2));
power_mW            = str2double(S(22,1:numSecondsData * 2));
chargeState         = str2double(S(23,1:numSecondsData * 2));
remainingCapacity   = str2double(S(24,1:numSecondsData * 2));

%% Build sensor-specific times for the FaunaTag's raw data types

fprintf('\tBuilding frequency-appropriate time scales and sensor structures...\n');

opticsFs = 250;
movementFs = 100;
magFs = 20;
powerFs = 2;

numOpticsSamples = length(led1);
opticsSamplingIncrement = 1/opticsFs;
opticsTime_s = 0:opticsSamplingIncrement:(numOpticsSamples/opticsFs)-opticsSamplingIncrement;
opticsTime = trialDateTimeLocal + seconds(opticsTime_s);

numMovementSamples = length(ax);
moveSamplingIncrement = 1/movementFs;
moveTime_s = 0:moveSamplingIncrement:(numMovementSamples/movementFs)-moveSamplingIncrement;
moveTime = trialDateTimeLocal + seconds(moveTime_s);

numMagSamples = length(mx);
magSamplingIncrement = 1/magFs;
magTime_s = 0:magSamplingIncrement:(numMagSamples/magFs)-magSamplingIncrement;
magTime = trialDateTimeLocal + seconds(magTime_s);

numPowerSamples = length(voltage);
powerSamplingIncrement = 1/powerFs;
powerTime_s = 0:powerSamplingIncrement:(numPowerSamples/powerFs)-powerSamplingIncrement;
powerTime = trialDateTimeLocal + seconds(powerTime_s);


%% Decimate optical (250 Hz) + accel, gyro, depth & temp (100 Hz) to 25 Hz

fprintf('\tDecimating optical, movement, depth & temp data to 25Hz...\n');

led1_25Hz = decimate(led1, 10);         % decimate 250 Hz by a factor of 10
led2_25Hz = decimate(led2, 10);
led3_25Hz = decimate(led3, 10);
led4_25Hz = decimate(led4, 10);

normLed1_25Hz = decimate(normLed1, 10);
normLed2_25Hz = decimate(normLed2, 10);
normLed3_25Hz = decimate(normLed3, 10);
normLed4_25Hz = decimate(normLed4, 10);

% build the log10 optical signals
log10Led1_25Hz = decimate(log10Led1, 10);
log10Led2_25Hz = decimate(log10Led2, 10);
log10Led3_25Hz = decimate(log10Led3, 10);
log10Led4_25Hz = decimate(log10Led4, 10);

log10active_25Hz = decimate(log10active, 10);
log10ambient_25Hz = decimate(log10ambient, 10);

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

odba_25Hz = sqrt(ax_25Hz.^2 + ay_25Hz.^2 + az_25Hz.^2);
normOdba_25Hz = normalize(sqrt(ax_25Hz.^2 + ay_25Hz.^2 + az_25Hz.^2));

%% Build an x-scale in seconds based on samples and decimated frequency

Fs = 25;
numSamples = length(led1_25Hz);
samplingFreqHz = Fs;
samplingIncrement = 1/samplingFreqHz;

time_s = 0:samplingIncrement:(numSamples/samplingFreqHz)-samplingIncrement;
% time = time_s';
time = trialDateTimeLocal + seconds(time_s);

%% Do a bandpass filter 250 Hz & 25 Hz optical data (heart rate filter)

fprintf('\tApplying a band-pass filter with the following properties:\n');

bwFiltOrder = 6;        % 20th order butterworth filter

% experiment with some low and high-pass frequencies
% my initial guesses for some good band-pass freqs
lowCut      = 0.3333;   % 0.3333 Hz = ~20 beats per minute
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

% filter optical data, using d' to row-major output
bpLed1 = filtfilt(d', led1);
bpLed2 = filtfilt(d', led2);
bpLed3 = filtfilt(d', led3);
bpLed4 = filtfilt(d', led4);

bpOdba = filtfilt(d', odba);

bpLed1_25Hz = filtfilt(d', led1_25Hz);
bpLed2_25Hz = filtfilt(d', led2_25Hz);
bpLed3_25Hz = filtfilt(d', led3_25Hz);
bpLed4_25Hz = filtfilt(d', led4_25Hz);

bpOdba_25Hz = filtfilt(d', odba_25Hz);

%%  Where are things now?
%   Raw data has been turned into relevant row-major sensor streams
%   Raw data time scale have been created
%   Decimated data at 25 Hz created from raw data for optics & kinematics
%   Decimated data time scale has been built
%   Band-passed raw and decimated optics have been created to screen for HR

%% Create a TAG struct that contains universal meta data for the record

%  Note: this struct contains elements common to both RAW and PRH

TAG = struct;
TAG.id = FaunaTagID;
TAG.startTimeLocal = trialDateTimeLocal;
TAG.startTimeUTC = trialDateTimeUTC;
TAG.utcOffset = utc;
TAG.latitude = lat;
TAG.longitude = long;
TAG.latitudeUnits = 'decimal degrees';
TAG.longitudeUnits = 'decimal degrees';
TAG.declination = decl;
TAG.declinationUnits = 'degrees';
TAG.magneticFieldStrength = magFieldStrength;
TAG.magFieldStrengthUnits = 'nT';
TAG.originalFileName = fileName;


%% Construct raw data sensor structures

% NOTE: raw sensor data structures have full name for sensor streams, 
% e.g.: OPTICS = optics for 250 Hz data streams. In contrast (see below),
% decimated data sensor structures (at a common 25 Hz) use just the
% upper-case letter, e.g.: the letter O for 25 Hz decimated optics 
% OPTICS = optics

fprintf('\tCreating structures for RAW sensor data streams...\n');

OPTICS = struct;
OPTICS.time = opticsTime;
OPTICS.time_s = opticsTime_s;
OPTICS.led1 = led1;
OPTICS.led2 = led2;
OPTICS.led3 = led3;
OPTICS.led4 = led4;
OPTICS.bpLed1 = bpLed1;
OPTICS.bpLed2 = bpLed2;
OPTICS.bpLed3 = bpLed3;
OPTICS.bpLed4 = bpLed4;
OPTICS.normled1 = normLed1;
OPTICS.normled2 = normLed2;
OPTICS.normled3 = normLed3;
OPTICS.normled4 = normLed4;
OPTICS.log10Led1 = log10Led1;
OPTICS.log10Led2 = log10Led2;
OPTICS.log10Led3 = log10Led3;
OPTICS.log10Led4 = log10Led4;
OPTICS.log10Active = log10active;
OPTICS.log10Ambient = log10ambient;
OPTICS.fs = opticsFs;

% K = kinematics (accel, gyro, pitch, roll, heading, odba)

KINEMATICS = struct;
KINEMATICS.time = moveTime;
KINEMATICS.time_s = moveTime_s;
KINEMATICS.ax = ax;
KINEMATICS.ay = ay;
KINEMATICS.az = az;
KINEMATICS.gx = gx;
KINEMATICS.gy = gy;
KINEMATICS.gz = gz;
KINEMATICS.pitch = pitch;
KINEMATICS.roll = roll;
KINEMATICS.heading = heading;
KINEMATICS.odba = odba;
KINEMATICS.normOdba = normOdba;
KINEMATICS.accelerometerUnits = 'g';
KINEMATICS.accelerometerSensitivity = '1 LSB/mg';
KINEMATICS.accelerometerScale = '±2g';
KINEMATICS.gyroscopeUnits = 'degrees per second (dps)';
KINEMATICS.gyroscopeSensitivity = 'LSB/dps';
KINEMATICS.gyroscopeScale = '±500 dps';
KINEMATICS.prhUnits = 'degrees';
KINEMATICS.odbaUnits = 'g';
KINEMATICS.referenceFrame = 'tag';   % no tag2whale conversion yet

PRESSURE = struct;
PRESSURE.depth = depth;
PRESSURE.temperature = temperature;
PRESSURE.fs = movementFs;
PRESSURE.depthUnits = 'meters';
PRESSURE.temperatureUnits = 'degrees celsius';

% MAG = magnetometer

MAG = struct;
MAG.time = magTime;
MAG.time_s = magTime_s;
MAG.mx = mx;
MAG.my = my;
MAG.mz = mz;
MAG.fs = magFs;
MAG.units = 'microTeslas';
MAG.sensitivity = '0.3 microTeslas';

% power

POWER = struct;
POWER.time = powerTime;
POWER.time_s = powerTime_s;
POWER.voltage = voltage;
POWER.current = current;
POWER.power = power_mW;
POWER.chargeState = chargeState;
POWER.remainingCapacity = remainingCapacity;
POWER.fs = powerFs;
POWER.utcOffset = utc;
POWER.voltageUnit = 'mV';
POWER.currentUnit = 'mA';
POWER.powerUnit = 'mW';
POWER.chargeUnit = '%';
POWER.capacityUnit = 'mAh';


%% Build structures for 25 Hz PRH file

fprintf('\tCreating structures for PRH (25 Hz) sensor data streams...\n');

% O = 25 Hz optics

O = struct;

O.led1 = led1_25Hz;
O.led2 = led2_25Hz;
O.led3 = led3_25Hz;
O.led4 = led4_25Hz;
O.bpLed1 = bpLed1_25Hz;
O.bpLed2 = bpLed2_25Hz;
O.bpLed3 = bpLed3_25Hz;
O.bpLed4 = bpLed4_25Hz;
O.normLed1 = normLed1_25Hz;
O.normLed2 = normLed2_25Hz;
O.normLed3 = normLed3_25Hz;
O.normLed4 = normLed4_25Hz;
O.log10Led1 = log10Led1_25Hz;
O.log10Led2 = log10Led2_25Hz;
O.log10Led3 = log10Led3_25Hz;
O.log10Led4 = log10Led4_25Hz;
O.log10active = log10active_25Hz;
O.log10ambient = log10ambient_25Hz;
O.fs = Fs;

% K = 25 Hz kinematics, depth and temperature

K = struct;

K.ax = ax_25Hz;
K.ay = ax_25Hz;
K.az = ax_25Hz;
K.gx = ax_25Hz;
K.gy = ax_25Hz;
K.gz = ax_25Hz;
K.pitch = pitch_25Hz;
K.roll = roll_25Hz;
K.heading = heading_25Hz;
K.odba = odba_25Hz;
K.normOdba = normOdba_25Hz;
K.fs = Fs;
K.accelerometerUnits = 'g';
K.accelerometerSensitivity = '1 LSB/mg';
K.accelerometerScale = '±2g';
K.gyroscopeUnits = 'degrees per second (dps)';
K.gyroscopeSensitivity = 'LSB/dps';
K.gyroscopeScale = '±500 dps';
K.prhUnits = 'degrees';
K.odbaUnits = 'g';
K.referenceFrame = 'tag';   % no tag2whale conversion yet

% P = pressure (depth and temperature)
P.depth = depth_25Hz;
P.temperature = temperature_25Hz;
P.fs = Fs;
P.depthUnits = 'meters';
P.temperatureUnits = 'degrees celsius';

% skipping magnetometer (20 Hz) and power data (2 Hz) data

% create a TIME struct
TIME = struct;
TIME.time = time;       % datetime version
TIME.time_s = time_s;   % time in 25ths of a second since tag activation
TIME.fs = Fs;


%% save all RAW data structs as a somewhat Dtag-compatible raw file

global TAG_PATHS;

rawFile = sprintf('%s/%s', TAG_PATHS.RAW, rawFileName);

if (exist(rawFile, 'file'))
    fprintf('\tRaw file already exists... skipping re-save for now.\n');
else
    fprintf('\tRaw file does not exist... creating...\n');
    save(rawFile, 'OPTICS', 'KINEMATICS', 'PRESSURE', ...
        'MAG', 'POWER', 'TAG');
    if (exist(rawFile, 'file'))
        fprintf('\tLooks like that file was created!\n');
    else
        fprintf('\tThere was some problem saving that raw file. Halting.\n');
    end
end


%% save decimated and aligned data to the PRH file

prhFile = sprintf('%s/%s', TAG_PATHS.PRH, prhFileName);

if (exist(prhFile, 'file'))
    
    fprintf('\tPRH file exists... ');
    
    prompt = 'over-write and re-create? y/n [n]: ';
    str = input(prompt,'s');
    if isempty(str)
        str = 'n';
    end
    
    if (strcmp(str, 'y'))
       fprintf('\n\tOver-writing and re-creating PRH file...\n');
       
    end
    
else
    fprintf('\tPRH file does not exist... creating...\n');
    save(prhFile, 'O', 'K', 'P', 'TIME', 'TAG');
    if (exist(prhFile, 'file'))
        fprintf('\tLooks like that file was created!\n');
    else
        fprintf('\tThere was some problem saving that raw file. Halting.\n');
    end
end

%% Okay! RAW and PRH files should have been created. Now what?
%   
%   1. recommend separation of all of the above steps as a separate script,
%       e.g.: 'FaunaData_PostProcessing.m'
%   
%   2. recommend writing a separate tool to partition multi-trial data
%       sets, call it 'FaunaData_Partioner.m'
%
%   3. recommend separation of:
%           a) all of the below figure plotting
%           b) 'FaunaTagDataSheetInput.m' 
%           c) 'FaunaTag_OpticalGaussianFilter.m' 
%      into 'FaunaData_Analysis.m'

%% Always have a look at raw optical data before doing anything

figure;

rawLed1 = subplot(411);
plot(opticsTime, led1, 'k-');
xlabel('Time (local)');
ylabel('Intensity');
title('LED1 (amb1)');
grid;

rawLed2 = subplot(412);
plot(opticsTime, led2, 'b-');
xlabel('Time, seconds');
ylabel('Intensity');
title('LED2 (\lambda = 1050nm)');
grid;

rawLed3 = subplot(413);
plot(opticsTime, led3, 'r-');
xlabel('Time (local)');
ylabel('LED3 (1200nm) Intensity');
title('LED3 (\lambda = 1200 nm)');
grid;

rawLed4 = subplot(414);
plot(opticsTime, led4, 'g-');
xlabel('Time (local)');
ylabel('Intensity');
title('LED4 (amb2)');
grid;

linkaxes([ rawLed1, rawLed2, rawLed3, rawLed4], 'x');

%% plot some log plots of active LEDs...

figure;

optPlot1 = subplot(311);
plot(optics.time, optics.log10Led2 * -1, 'b-'); 
hold on; 
plot(optics.time, optics.log10Led3 * -1, 'r-'); 
hold off;
xlabel('Time (local)');
ylabel('Log10 Active LEDs')
title('Log10 LED2 (\lambda = 1050 nm) and LED3 (\lambda = 1200 nm)');
grid;

optPlot2 = subplot(312);
plot(optics.time, optics.log10Active * -1, 'k-');
xlabel('Time (local)');
ylabel('Tissue saturation index');
title('Log10 LED2 (\lambda = 1050 nm) / Log10 LED3 (\lambda = 1200 nm)');
grid;

optPlot3 = subplot(313);
plot(optics.time, optics.log10Ambient * -1, 'g-');
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
plot(pitch, 'b-');
xlabel('Time, seconds');
ylabel('Pitch, °');
grid;

figRoll = subplot(312);
plot(roll, 'r-');
xlabel('Time, seconds');
ylabel('Roll, °');
grid;

figHeading = subplot(313);
plot(heading, 'k-');
xlabel('Time, seconds');
ylabel('Heading, °');
grid;

linkaxes([ figPitch, figRoll, figHeading], 'x');

%% overplot undecimated norm'd rawLed2 & rawLed3

figure;

plot(opticsTime_s, normLed1 * -1, 'k-');
hold on;
plot(opticsTime_s, normLed2 * -1, 'b-');
plot(opticsTime_s, normLed3 * -1, 'r-');
plot(opticsTime_s, normLed4 * -1, 'g-');
hold off;
xlabel('Time, seconds');
ylabel('Norm LED Intensities');
grid;







%% Plot some decimated data to know where to manually trim out bad data

figure;

dec1 = subplot(411);
plot(time, led1_25Hz, 'k-');
hold on;
plot(time, led2_25Hz * -1, 'b-');
plot(time, led3_25Hz * -1, 'r-');
plot(time, led4_25Hz, 'g-');
hold off;
xlabel('Time (local)');
ylabel('LED Intensity');
grid;
title('Decimated optical and kinematic data');

dec2 = subplot(412);
plot(time, pitch_25Hz, 'b-');
hold on;
plot(time, roll_25Hz, 'r-');
plot(time, heading_25Hz, 'g-');
xlabel('Time, seconds');
ylabel('Degrees');
title('Pitch · Roll · Heading');
ylim([-200 200]);
grid;
hold off;

dec3 = subplot(413);
plot(time, odba_25Hz, 'k-');
xlabel('Time (local)');
ylabel('Dirty ODBA');
grid;

dec4 = subplot(414);
plot(time, depth_25Hz, 'b-');
ylabel('Depth, meters');
xlabel('Time (local)');
set(gca, 'YDir','reverse');
grid;

linkaxes([dec1 dec2 dec3 dec4], 'x');

%% Resize sensor data to get rid of pre-tag-on + post-tag-off

    

%% try looking at some normalized data for signal:noise comparisons...


figure;

normDecPlot1a = subplot(411);

plot(time_s, led1_d_norm * -1, 'k-');
hold on;
plot(time_s, led2_d_norm * -1, 'b-');
plot(time_s, led3_d_norm * -1, 'r-');
plot(time_s, led4_d_norm * -1, 'g-');
hold off;
xlabel('Time, seconds');
title('Overlaid normalized optical signals');
legend('amb1','1050 nm','1200 nm','amb2','Location','southeast');
grid;

normDecPlot1b = subplot(412);
plot(time_s, normOdba_25Hz, 'k-');
xlabel('Time, seconds');
ylabel('ODBA');
title('Movement proxy');
grid;

normDecPlot1c = subplot(413);
plot(time_s, depth_25Hz, 'b-');
xlabel('Time, seconds');
ylabel('Meters');
title('Depth');
set(gca, 'YDir','reverse');
grid;

normDecPlot1d = subplot(414);
plot(time_s, temperature_25Hz, 'r-');
xlabel('Time, seconds');
ylabel('°C');
title('Temperature');
grid;


linkaxes([ normDecPlot1a normDecPlot1b normDecPlot1c normDecPlot1d], 'x');


%% create an optics vector and compute a mean of all sampled wavelengths

opticsVector        = [led1_d_norm; led2_d_norm; led3_d_norm; led4_d_norm];
opticalMean         = mean(opticsVector);

% one approach: take the mean of both active LED channels
activeOpticalMean   = mean(opticsVector(2:3,:)); 

figure;
plot(time_s, opticalMean, 'k-');
hold on;
plot(time_s, activeOpticalMean, 'b-');
hold off;
title('Mean optical signals - active vs. ambient');
ylabel('Norm decimated optical signals');
xlabel('Time, seconds');
grid;
legend('Ambient channels','Active channels');

% anoother approach: sum the two ambient and two active LED channels
addAmbientOpticalNormD  = opticsVector(1,:) + opticsVector(4,:);
addActiveOpticalNormD   = opticsVector(2,:) + opticsVector(3,:);
activeSubAmbientNormD = addActiveOpticalNormD - addAmbientOpticalNormD;

figure;
plot(time_s, addAmbientOpticalNormD, 'k-');
hold on;
plot(time_s, addActiveOpticalNormD, 'b-');
plot(time_s, activeSubAmbientNormD, 'r-');
hold off;
title('Mean optical signals - active vs. ambient');
ylabel('Norm decimated optical signals');
xlabel('Time, seconds');
grid;
legend('Ambient channels','Active channels');

%% normalize decimated accel & gyro data

ax_d_norm   = normalize(ax_25Hz);
ay_d_norm   = normalize(ay_25Hz);
az_d_norm   = normalize(az_25Hz);

gx_d_norm   = normalize(gx_25Hz);
gy_d_norm   = normalize(gy_25Hz);
gz_d_norm   = normalize(gz_25Hz);

figure;

normKin1 = subplot(211);
plot(time_s, ax_d_norm, 'b-');
hold on;
plot(time_s, ay_d_norm, 'r-');
plot(time_s, az_d_norm, 'g-');
hold off;
xlabel('Time, seconds');
ylabel('Normalized accel');
title('Norm accel');
legend('ax','ay','az');
grid;

normKin2 = subplot(212);
plot(time_s, gx_d_norm, 'b-');
hold on;
plot(time_s, gy_d_norm, 'r-');
plot(time_s, gz_d_norm, 'g-');
hold off;
xlabel('Time, seconds');
ylabel('Normalized gyro');
title('Norm gyro');
legend('gx','gy','gz');
grid;

linkaxes([ normKin1 normKin2 ], 'x');



%% create an accel and gyro vector and compute the mean using all axes

accelVector = [ax_d_norm; ay_d_norm; az_d_norm];
gyroVector  = [gx_d_norm; gy_d_norm; gz_d_norm];

accelMean   = mean(accelVector);
gyroMean    = mean(gyroVector);

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

%  plot overalls of norm'd optical mean, with norm'd accel & gyro mean

figure;
plot(time_s, opticalMean, 'b-', 'LineWidth', 2);
hold on;
plot(time_s, accelMean, 'r-');
plot(time_s, gyroMean, 'm-');
plot(time_s, normOdba_25Hz, 'k-');
hold off;
grid;
xlabel('Time, seconds');
ylabel('Mean optical accel & gyro');
title('Mean normalized decimated sensor data');




%% Diagnostic optical plots

fprintf('\tPlotting raw LED1-4 signals...\n');

figure;

led1Plot = subplot(411);
plot(opticsTime_s, led1, 'g-');
title('LED1');
xlabel('Time, seconds');
ylabel('Intensity');
grid;

led2Plot = subplot(412);
plot(opticsTime_s, led2, 'b-');
title('LED2');
xlabel('Time, seconds');
ylabel('Intensity');
grid;

led3Plot = subplot(413);
plot(opticsTime_s, led3, 'r-');
title('LED3');
xlabel('Time, seconds');
ylabel('Intensity');
grid;

led4Plot = subplot(414);
plot(opticsTime_s, led4, 'g-');
title('LED4');
xlabel('Time, seconds');
ylabel('Intensity');
grid;

linkaxes([led1Plot led2Plot led3Plot led4Plot], 'x');

%% Diagnostic optical plots with theoretical ambient (LED1) subtraction

figure;

led1Plot = subplot(411);
plot(opticsTime_s, led1, 'k-');
title('LED1');
xlabel('Time, seconds');
ylabel('Intensity');
grid;

led2Plot = subplot(412);
plot(opticsTime_s, led2 - led1, 'b-');
title('LED2 - LED1');
xlabel('Time, seconds');
ylabel('Intensity');
grid;

led3Plot = subplot(413);
plot(opticsTime_s, led3 - led1, 'r-');
title('LED3 - LED1');
xlabel('Time, seconds');
ylabel('Intensity');
grid;

led4Plot = subplot(414);
plot(opticsTime_s, led4 - led1, 'g-');
title('LED4');
xlabel('Time, seconds');
ylabel('Intensity');
grid;

linkaxes([led1Plot led2Plot led3Plot led4Plot], 'x');

%% Do a bandpass filter on 250 Hz optical data and decimate a copy to 25 Hz

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


% [A, B, C, D] = butter(bwFiltOrder/2, [lowCut highCut] / (opticsFs/2));
% 
% d = designfilt('bandpassiir','FilterOrder',bwFiltOrder, ...
%     'HalfPowerFrequency1',lowCut,'HalfPowerFrequency2',highCut, ...
% 	'SampleRate',opticsFs);
% 
% sos = ss2sos(A,B,C,D);
% fvt = fvtool(sos,d,'Fs', opticsFs);
% legend(fvt,'butter','designfilt');

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

bpLed1 = filtfilt(d, led1');
bpLed2 = filtfilt(d, led2');
bpLed3 = filtfilt(d, led3');
bpLed4 = filtfilt(d, led4');

bpLed1_25Hz = filtfilt(d, led1_25Hz');
bpLed2_25Hz = filtfilt(d, led2_25Hz');
bpLed3_25Hz = filtfilt(d, led3_25Hz');
bpLed4_25Hz = filtfilt(d, led4_25Hz');

bpOdba_25Hz = filtfilt(d, odba_25Hz');

%% plot band-passed high resolution optical data

fprintf('\tPlotting band-passed optical data for LEDs 1-4\n');

figure;

bpFigA = subplot(411);
plot(opticsTime_s, bpLed1, 'k-');
xlabel('Time, seconds');
ylabel('LED1 ambient');
title('Bandpass filtered 250 Hz optical data');
grid;

bpFigB = subplot(412);
plot(opticsTime_s, bpLed2, 'b-');
xlabel('Time, seconds');
ylabel('LED2');
grid;

bpFigC = subplot(413);
plot(opticsTime_s, bpLed3, 'r-');
xlabel('Time, seconds');
ylabel('LED3');
grid;

bpFigD = subplot(414);
plot(opticsTime_s, bpLed4, 'g-');
xlabel('Time, seconds');
ylabel('LED4');
grid;

linkaxes( [bpFigA bpFigB bpFigC bpFigD], 'x');

% plot all of the bpLeds overlaid

figure;
plot(opticsTime_s, bpLed1, 'k-');
hold on;
plot(opticsTime_s, bpLed2, 'b-');
plot(opticsTime_s, bpLed3, 'r-');
plot(opticsTime_s, bpLed4, 'g-');
hold off;
xlabel('Time, seconds');
ylabel('Intensity, LEDs 1-4');
title('All Bandpassed LED Samples Overlaid');
grid;

%% plot decimated overlaid optical channels and band-passed ODBA

figure;
bpDec1 = subplot(311);
plot(time_s, bpLed1_25Hz, 'k-');
hold on;
plot(time_s, bpLed2_25Hz * -1, 'b-');
plot(time_s, bpLed3_25Hz * -1, 'r-');
plot(time_s, bpLed4_25Hz, 'g-');
hold off;
xlabel('Time, seconds');
ylabel('Band-passed optics');
title('Overlaid band-passed optical channels');
grid;

bpDec2 = subplot(312);
plot(time_s, normOdba_25Hz * -1, 'k-');
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