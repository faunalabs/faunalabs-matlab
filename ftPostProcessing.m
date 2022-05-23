%% FaunaData Post-Processing -
%
%       Created and written by Dave Haas, 18 June 2021
%
%       Based on FaunaDataDataWorkup.m, which was created by Dave Haas
%       between 17 June 2020 - 18 June 2021. This unpacks FaunaTag data
%       from its CSV format, puts it into row-major data for Matlab,
%       partitions the data into RAW structures (default sampling rates)
%       and PRH structures (25 Hz decimated data).
%
%       Dependency: csvproc from dtagtools by Mark Johnson

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

%%  begin FaunaTag post-processing by choosing a tag, date, and trial set

fprintf('\n\n---=== Starting FaunaData Post-processing! ===---\n\n');
fprintf('Select a FaunaTag .bin file for processing...\n');

[fileName,pathName] = uigetfile('*.bin');

% % now, build a complete targetFile for open...
targetFile = sprintf('%s%s', pathName, fileName);

fileExists = isfile(targetFile);

workingDir = sprintf('%s/FaunaLabs/', pwd);

truncdPathName =  pathName(1:end-1);

%% parse fileName, make date, time, utcOffset, etc., generate prhFileName

fprintf('Parsing filename into date-time components...\n');

if (strncmp(fileName, 'ft110', 5))
    FaunaTagID = 110;
elseif (strncmp(fileName, 'ft111', 5))
    FaunaTagID = 111;
elseif (strncmp(fileName, 'ft112', 5))
    FaunaTagID = 112;
elseif (strncmp(fileName, 'ft200', 5))
    FaunaTagID = 200;
elseif (strncmp(fileName, 'ft201', 5))
    FaunaTagID = 201;
elseif (strncmp(fileName, 'ft202', 5))
    FaunaTagID = 202;
elseif (strncmp(fileName, 'ft203', 5))
    FaunaTagID = 203;    
elseif (strncmp(fileName, 'ft204', 5))
    FaunaTagID = 204;
end


YYYY = str2double(extractBetween(fileName, 7, 10));
YY   = str2double(extractBetween(fileName, 9, 10));
MM   = str2double(extractBetween(fileName, 11,12));
DD   = str2double(extractBetween(fileName, 13,14));

hh   = str2double(extractBetween(fileName, 16,17));
mm   = str2double(extractBetween(fileName, 18,19));
ss   = str2double(extractBetween(fileName, 20,21));

utcStr = string(extractBetween(fileName, 22, 24));
utc	   = str2double(extractBetween(fileName, 22,24));

J    = jday([YYYY MM DD]);

trialDate = [YYYY MM DD];

trialDateTimeLocal = datetime([YYYY MM DD hh mm ss], ...
    'TimeZone', utcStr, 'Format','d-MMM-y HH:mm:ss Z');

trialDateTimeUTC   = datetime([YYYY, MM, DD, hh + (utc * -1), mm, ss], ...
    'TimeZone', utcStr, 'Format','d-MMM-y HH:mm:ss Z');

%% Enter animal details from this trial

fprintf('Enter animal details for this trial...\n');

animalName = input('Animal name: ','s');

while (isempty(animalName))
    animalName = input('Animal name: ','s');
end
    
switch(animalName)
    case 'Kolohe'
        animalId = '6JK5';
        animalSex = 'm';
        animalAge = 26;
        animalMass = 0;
        animalGenus = 'Tursiops';
        animalSpecies = 'truncatus';
        gsCode = 'tt';
    case 'Liho'
        animalId = '01L5';
        animalSex = 'm';
        animalAge = 26;
        animalMass = 0;
        animalGenus = 'Tursiops';
        animalSpecies = 'truncatus';
        gsCode = 'tt';
    case 'Hoku'
        animalId = '63HF';
        animalSex = 'm';
        animalAge = 30;
        animalMass = 0;
        animalGenus = 'Tursiops';
        animalSpecies = 'truncatus';
        gsCode = 'tt';
    case 'Hua'
        animalId = '83H1';
        animalSex = 'm';
        animalAge = 13;
        animalMass = 0;
        animalGenus = 'Tursiops';
        animalSpecies = 'truncatus';
        gsCode = 'tt';
    case 'Noa'
        animalId = '9ON6';
        animalSex = 'm';
        animalAge = 22;
        animalMass = 0;
        animalGenus = 'Tursiops';
        animalSpecies = 'truncatus';
        gsCode = 'tt';
    case 'Lono'
        animalId = '9FL3';
        animalSex = 'm';
        animalAge = 36;
        animalMass = 0;
        animalGenus = 'Tursiops';
        animalSpecies = 'truncatus';
        gsCode = 'tt';
    case 'Keto'
        animalId = '3756';
        animalMass = 3756;
        animalSex = 'm';
        animalAge = 28;
        animalGenus = 'Orcinus';
        animalSpecies = 'orca';
        gsCode = 'oo';
    otherwise
        fprintf('Animal details not present for that name, enter manually...\n');
        animalId = input('Animal ID (use managed care ID): ','s');
        animalSex = input('Animal sex (m/f): ','s');
        animalAge = str2double(input('Animal age: ', 's'));
        animalMass = str2double(input('Animal mass in kg (enter 0 if unknown): ','s'));
        animalGenus = input('Animal genus: ','s');
        animalSpecies = input('Animal species: ','s');
        gsCode = sprintf('%s%s', string(lower(extract(animalGenus,1))), ...
            string(lower(extract(animalSpecies,1)) ) );
               
end

fprintf('AnimalId: %s\n', animalId);
fprintf('Animal sex: male\n');
fprintf('Animal age: %d\n', animalAge);
fprintf('Animal taxon: %s %s\n', animalGenus, animalSpecies);
fprintf('Animal taxon code: %s\n', gsCode);

%%

animalFastedStr = input('Animal fasted? (y/n) [n]: ', 's');
    
if (isempty(animalFastedStr))
    animalFasted = 0;
elseif (strcmp ( lower(animalFastedStr), 'n'))
    animalFasted = 0;
elseif (strcmp ( lower(animalFastedStr), 'y'))
    animalFasted = 1;
else
    animalFasted = 0;
end

fprintf('Animal fasted state: %d\n', animalFasted);

%%

trialLoc = str2double(input('Specify trial location (1 = DQO, 2 = LP, 3 = OG, 4 = Vegas, 0 = enter manually): ', 's'));

switch(trialLoc)
    
    case 1   % Dolphin Quest Oahu
       
        lat = 21.272037;            % decimal degrees at Dolphin Quest Oahu
        long = -157.773116;         % decimal degrees at Dolphin Quest Oahu
        decl = 9.48;                % units: degrees, on 13 May 2021
        magFieldStrength = 34797.5; % field strength, units: nanoTeslas (nT)
        validMetaData = true;
    
    case 2  % Loro Parque, Tenerife

        lat = 28.4082;            % decimal degrees at Dolphin Quest Oahu
        long = -16.5643;         % decimal degrees at Dolphin Quest Oahu
        decl = -4.32;                % units: degrees, on 13 May 2021
        magFieldStrength = 38675; % field strength, units: nanoTeslas (nT)
        validMetaData = true;
    
    case 3  % Oceanografic, Valencia
        
        lat = 39.453070;            % decimal degrees at Dolphin Quest Oahu
        long = -0.34721;            % decimal degrees at Dolphin Quest Oahu
        decl = 0.95;                % units: degrees, on 13 May 2021
        magFieldStrength = 44821;   % field strength, units: nanoTeslas (nT)
        validMetaData = true;

        
    case 4  % __, Las Vegas
        
        validMetaData = false;
        
    case 0

    % enter metadata manually if not DQO or LP trials ...
    
    while(validMetaData == false)

        fprintf('\n----- Enter tag meta data information manually... -----\n\n');

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

fprintf('Trial latitude (DD.dddd): %3.4f\n', lat);
fprintf('Trial longitude (DD.dddd): %3.4f\n', long);
fprintf('Trial declination angle (°): %3.4f\n', decl);
fprintf('Trial magnetic field strength (nT): %5.1f\n', magFieldStrength);

%% Build a struct for the trial's focal animal

ANIMAL = struct;
ANIMAL.name = animalName;
ANIMAL.id = animalId;
ANIMAL.age = animalAge;
ANIMAL.sex = animalSex;
ANIMAL.mass = animalMass;
ANIMAL.genus = animalGenus;
ANIMAL.species = animalSpecies;
ANIMAL.gsCode = gsCode;
ANIMAL.fasted = animalFasted;

%% construct file names for writing all this out...

totd = input('Trial of the day (specify a for 1st, b for 2nd, c for 3rd, ...): ','s');

tag = sprintf('%s%d_%03d%s', gsCode,YY,J,totd);

tsvFileName = sprintf('%s%d_%03d%s.tsv',gsCode,YY,J,totd);
csvFileName = sprintf('%s%d_%03d%s.txt',gsCode,YY,J,totd);
rawFileName = sprintf('%s%d_%03d%sraw.mat',gsCode,YY,J,totd);
prhFileName = sprintf('%s%d_%03d%sprh.mat',gsCode,YY,J,totd);
metadataFileName = sprintf('%s%d_%03d%smeta.mat',gsCode,YY,J,totd);

fprintf('\nTag ID: %s\n', tag);
fprintf('TSV filename: %s\n', tsvFileName);
fprintf('CSV filename: %s\n', csvFileName);
fprintf('RAW (undecimated data) filename: %s\n', rawFileName);
fprintf('PRH (50 Hz decimated) filename: %s\n', prhFileName);
fprintf('METADATA filename: %s\n', metadataFileName);


%% Create the TSV and CSV files

execString = sprintf('%stagOS_dataPrep -i ''%s'' -o ''%s%s''', ...
    workingDir, targetFile, pathName, tag);

status = system(execString);

%%

fullCsvFileName = sprintf('%s%s', pathName, csvFileName);

csvFileExists = isfile(fullCsvFileName);

if (csvFileExists)
    fprintf('CSV file exists... reading in data...\n');
    youShallNotPass = false;
    % load targetFile into S
    S = csvproc(fullCsvFileName, []);    
else
    fprintf('File you specified does not seem to exist...\n\n');
    youShallNotPass = true;
    return;
end


%% parse S into usable data components

if (youShallNotPass)
    fprintf('The data file was missing. Find the right now and try again!\n');
    return;
end

fprintf('Parsing data into usable components...\n');

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

% is there an ODBA-like version for gyroscope? there is now...
% call it Overall Dynamic Gyroscopic Signal
odav = sqrt(gx.^2 + gy.^2 + gz.^2);
normOdav = normalize(sqrt(gx.^2 + gy.^2 + gz.^2));

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

fprintf('Building frequency-appropriate time scales and sensor structures...\n');

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

%% Create MATLAB time tables for each of these

fprintf('Constructing timetables of raw FaunaTag data...\n');

optics = timetable(opticsTime', opticsTime_s', led1', led2', led3', ...
    led4', 'VariableNames', {'time_s','led1','led2','led3','led4'} );

kinematics = timetable(moveTime', moveTime_s', ax', ay', az', ...
    gx', gy', gz', pitch', roll', heading', odba', odav', ...
    'VariableNames', {'time_s','ax','ay','az','gx','gy','gz','pitch', ...
    'roll', 'heading', 'odba', 'odav'} );

pressure = timetable(moveTime', moveTime_s', depth', temperature', ...
    'VariableNames', {'time_s','depth','temperature' } );

mag = timetable(magTime', magTime_s', mx', my', mz', 'VariableName', ...
    {'time_s','mx','my','mz'} );

power = timetable(powerTime', powerTime_s', voltage', current', ...
    power_mW', chargeState', remainingCapacity', 'VariableNames', ...
    {'time_s','voltage','current','power','chargeState', ...
    'remainingCapacity' } );

%% Create a TAG struct that contains universal meta data for the record

%  IS THIS EVEN RELEVANT NOW?!?!  ALL THIS DATA IS IN 'TRIAL' META DATA...

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


%% Continue with building out trial meta data...

fprintf('Answer these questions to create metadata file of trial details...\n');
fprintf('Body location of tag attachment, choose from one of these:\n');
fprintf('\t1: post-nuchal\n');
fprintf('\t2: anterio-dorsal\n');
fprintf('\t3R: dorso-lateral (right)\n');
fprintf('\t3L: dorso-lateral (left)\n');
fprintf('\t4: ventral (chest)\n');
fprintf('\t5R: right saddle-patch\n');
fprintf('\t5L: left saddle-patch\n');

attachId = 0;

while (attachId == 0)
    
    attachIdStr = input('Attachment ID: ','s');    
    
    switch(attachIdStr)
        case '1'
            attachId = 1;
            attachSideId = 0;
        case '2'
            attachId = 2;
            attachSideId = 0;
        case '3R'
            attachId = 3;
            attachSideId = 1;
        case '3L'
            attachId = 3;
            attachSideId = 2;
        case '4'
            attachId = 4;
            attachSideId = 0;
        case '5R'
            attachId = 5;
            attachSideId = 1;                        
        case '5L'
            attachId = 5;
            attachSideId = 2;            
        otherwise
            attachId = 0;
            attachSideId = 0;
    end
end



%% tag activate, on, off, deactivate times, needed to create observe offsets

% tag activation datetime

fprintf('Enter activation, tag-on, tag-off, and deactivation times from data sheet...\n');

% tagTimeZoneStr = input('UTC Offset in hours (e.g.: GMT: +00; Hawaii: -10) [-10]: ','s');
% if (isempty(tagTimeZoneStr))
%    tagTimeZoneStr = '-10';
% end
% tz = str2double(tagTimeZoneStr);

cueTimeText = 'Tag activation local time (hh:mm:ss): ';
tagActivationTimeStr = input(cueTimeText,'s');

validCueTime = false;        
while (validCueTime == false)
    hms = str2double(strsplit(tagActivationTimeStr,':'));
    if ( ( numel(hms) == 3 ) && le(hms(1),23) && ge(hms(1),0) && ...
            le(hms(2),59) && ge(hms(2),0) && ... 
            le(hms(3),59) && ge(hms(3),0) )
        % this should escape the while condition...
        validCueTime = true;
    else
        % this should create another while() cycle
        validCueTime = false;
        tagActivationTimeStr = input(cueTimeText, 's');
    end    
end
% validation above makes this next line unnecessary
% [h, m, s] = hms(datetime(tagActivationTimeStr,'InputFormat','HH:mm:ss'));
tagActivateDateTime = datetime([YYYY MM DD hms], ...
    'TimeZone',utcStr,'Format','d-MMM-y HH:mm:ss Z');
fprintf('Tag activation datetime: %s\n', tagActivateDateTime);

% tag on datetime

cueTimeText = 'Tag ON local time (hh:mm:ss): ';
tagOnTimeStr = input(cueTimeText,'s');
validCueTime = false;        
while (validCueTime == false)
    hms = str2double(strsplit(tagOnTimeStr,':'));
    if ( ( numel(hms) == 3 ) && le(hms(1),23) && ge(hms(1),0) && ...
            le(hms(2),59) && ge(hms(2),0) && ... 
            le(hms(3),59) && ge(hms(3),0) )
        % this should escape the while condition...
        validCueTime = true;
    else
        % this should create another while() cycle
        validCueTime = false;
        tagOnTimeStr = input(cueTimeText, 's');
    end    
end
% validation above makes this next line unnecessary
% [h, m, s] = hms(datetime(tagActivationTimeStr,'InputFormat','HH:mm:ss'));
tagOnDateTime = datetime([YYYY MM DD hms], ...
    'TimeZone',utcStr,'Format','d-MMM-y HH:mm:ss Z');
fprintf('Tag on datetime: %s\n', tagOnDateTime);

% tag off datetime
cueTimeText = 'Tag OFF local time (hh:mm:ss): ';
tagOffTimeStr = input(cueTimeText,'s');
validCueTime = false;        
while (validCueTime == false)
    hms = str2double(strsplit(tagOffTimeStr,':'));
    if ( ( numel(hms) == 3 ) && le(hms(1),23) && ge(hms(1),0) && ...
            le(hms(2),59) && ge(hms(2),0) && ... 
            le(hms(3),59) && ge(hms(3),0) )
        % this should escape the while condition...
        validCueTime = true;
    else
        % this should create another while() cycle
        validCueTime = false;
        tagOffTimeStr = input(cueTimeText, 's');
    end    
end
tagOffDateTime = datetime([YYYY MM DD hms], ...
    'TimeZone',utcStr,'Format','d-MMM-y HH:mm:ss Z');
fprintf('Tag off datetime: %s\n', tagOffDateTime);

% tag deactivate datetime

cueTimeText = 'Tag deactivation local time (hh:mm:ss): ';
tagDeactivateTimeStr = input(cueTimeText,'s');
validCueTime = false;        
while (validCueTime == false)
    hms = str2double(strsplit(tagDeactivateTimeStr,':'));
    if ( ( numel(hms) == 3 ) && le(hms(1),23) && ge(hms(1),0) && ...
            le(hms(2),59) && ge(hms(2),0) && ... 
            le(hms(3),59) && ge(hms(3),0) )
        % this should escape the while condition...
        validCueTime = true;
    else
        % this should create another while() cycle
        validCueTime = false;
        tagDeactivateTimeStr = input(cueTimeText, 's');
    end    
end
tagDeactivateDateTime = datetime([YYYY MM DD hms], ...
    'TimeZone',utcStr,'Format','d-MMM-y HH:mm:ss Z');
fprintf('Tag deactivation datetime: %s\n', tagDeactivateDateTime);


%% Enter trial participants and context notes from this trial...

trialPIStr = input('Enter the name or initials of the trial PI [DKH]: ','s');
if (isempty(trialPIStr))
    trialPIStr = 'DKH';
end

trialParticipantsStr = input('Enter name(s) of trial participants [Nicole West]: ','s');
if (isempty(trialParticipantsStr))
    trialParticipantsStr = 'Nicole West';
end

trialNotesStr = input('Enter any notes about this trial: ','s');
if (isempty(trialNotesStr))
    trialNotesStr = 'No context notes on this log sheet';
end

%% Build a TAG structure of tag-specific meta data from this trial

fprintf('Creating TRIAL metadata structure...\n');

TRIAL = struct;
TRIAL.tag = tag;
TRIAL.tagNumber = FaunaTagID;
TRIAL.attachmentLocation = attachId;
TRIAL.attachmentSide = attachSideId;
TRIAL.dataFile = fileName;

TRIAL.tagStartTimeLocal = trialDateTimeLocal;
TRIAL.tagStartTimeUTC = trialDateTimeUTC;

TRIAL.utcOffset = utc;
TRIAL.loggedActivationTime = tagActivateDateTime;
TRIAL.loggedTagOnTime = tagOnDateTime;
TRIAL.loggedTagOffTime = tagOffDateTime;
TRIAL.loggedDeactivationTime = tagDeactivateDateTime;
TRIAL.tagTimelogTimeDifference = TRIAL.tagStartTimeLocal - ...
    TRIAL.loggedActivationTime;

TRIAL.latitude = lat;
TRIAL.longitude = long;
TRIAL.declination = decl;
TRIAL.magFieldStrength = magFieldStrength;

TRIAL.investigatorName = trialPIStr;
TRIAL.participantNames = trialParticipantsStr;
TRIAL.notes = trialNotesStr;

TRIAL.SUBJECT = ANIMAL;

%% Time cue entries (e.g.: for breath interval timing, etc.)

fprintf('Now is the time to enter detailed time cues, for timing of events like\n');
fprintf('transitions between baseline and apnea phases of trials, times of noted\n');
fprintf('ventilation events, etc.\n');
doTimeCueEntry  = lower(input('Any detailed time cues to enter (y/n)? ','s'));

if (strcmp(doTimeCueEntry,'y'))
   
    fprintf('\nEnter detailed time cues for this trial...\n'); 

    timeCueEntryComplete = false;
    
    if ( exist('CUE','var') )   % if CUEs entry accidentally exited...
    
        cueIndex = numel(CUE.time); % pick up where it left off...
        
    else
        
        cueIndex = 0;

        CUE = struct;
        TRIAL.loggedTagOnTime = tagOnDateTime;
        CUE.type = cell('');
        CUE.comment = cell('');
        CUE.tagTimeOffset = TRIAL.tagTimelogTimeDifference
        CUE.investigator = TRIAL.investigatorName;
        CUE.timingDevice = 'watch or stopwatch';
    
    end

    while( ~timeCueEntryComplete )
        
        cueIndex = cueIndex + 1;
        
        fprintf('\n---------------------------------------------------\n');
        
        % enter the time for this cueType at this cueIndex
        cueTimeText = sprintf('Cue #%d time (hh:mm:ss): ', cueIndex);
        cueTime = input(cueTimeText, 's');
        % check cueTime has valid ##:##:## structure; repeat until good
        validCueTime = false;
        while (validCueTime == false)
            hms = str2double(strsplit(cueTime,':'));
            if ( ( numel(hms) == 3 ) && le(hms(1),23) && ge(hms(1),0) && ...
                    le(hms(2),59) && ge(hms(2),0) && ... 
                    le(hms(3),59) && ge(hms(3),0) )
                % this should escape the while condition...
                validCueTime = true;
            else
                % this should create another while() cycle
                validCueTime = false;
                cueTime = input(cueTimeText, 's');
            end    
        end
        
        CUE.time(cueIndex) = datetime([YYYY MM DD hms], ...
            'TimeZone',utcStr,'Format','d-MMM-y HH:mm:ss Z');
        
        
        % enter the cue type code for this cueIndex
        fprintf('Cue #%d code options:\n', cueIndex);
        fprintf('\t1 = start of `chill` phase\n');
        fprintf('\t2 = start of free-breath (pre-breath-hold) phase\n');
        fprintf('\t3 = start of breath-hold phase\n');
        fprintf('\t4 = start of free-breathe recovery (post-breath-hold)\n');
        fprintf('\t* = respiration (visual or acoustic confirmation [enter])\n');
        fprintf('\to = other cue (note in the comment for this cue)\n');
        cueTypeText = sprintf('Cue #%d code: ', cueIndex);
        cueType = input(cueTypeText,'s');
        if (isempty(cueType))
            cueType = '*';
        end
        CUE.type(cueIndex) = {cueType};
        
        
        % deal with any cueTime context comments and/or 'o' coded cues...
        cueCommentText = sprintf('Cue #%d comments [enter to leave blank]: ', cueIndex);
        cueComment = input(cueCommentText,'s');
        CUE.comment(cueIndex) = {cueComment};
        
        % check with user for more time cue entries
        complete = lower(input('Enter another time cue? (y/n) [y]? ', 's'));
        if (strcmp(complete,'y') || isempty(complete) ) % not finished...
            continue;
        else
            timeCueEntryComplete = true;
            fprintf('Time cue entry completed. Moving on...\n');
        end
        
        
    end
   
else
    
   fprintf('\nSkipping cue entry for this trial...\n'); 
   
    CUE = struct;
    CUE.time = TRIAL.tagStartTimeLocal;
    CUE.type = cell('');
    CUE.comment = cell('');
    CUE.tagTimeOffset = duration([0 0 0]);
    CUE.investigator = TRIAL.investigatorName;
    CUE.timingDevice = 'watch or stopwatch';

end


%% test sort for out-of-order CUE.time entries

if (numel(CUE.time) > 1)
    [sortedCueTimes, sortedIndex] = sort(CUE.time);
    sortedCueTypes = CUE.type(sortedIndex);
    sortedCueComments = CUE.comment(sortedIndex);
    CUE.time = sortedCueTimes;
    CUE.type = sortedCueTypes;
    CUE.comment = sortedCueComments;
end

if (numel(CUE.time) > 1) 
    fprintf('The content of the CUE structure...\n');
    for i = 1:numel(CUE.type)
        fprintf('Times: %s | Cues: %s | Comments: %s\n', ... 
            CUE.time(i), string(CUE.type(i)), string(CUE.comment(i)) );
    end
else
    fprintf('The content of the CUE structure is empty (no entries).\n');    
end

fprintf('\n');

    
%% Build the TRIAL struct contain TAG SUBJECT & CUE

TRIAL.CUES = CUE;

disp(TRIAL);
disp(TRIAL.SUBJECT);
disp(TRIAL.CUES);


%% Write this out to metadata file

% global TAG_PATHS;

if ( (exist(TAG_PATHS.METADATA,'dir')) && (contains(TAG_PATHS.METADATA,'metadata')) )
    
    fullFileName = sprintf('%s/%smeta.mat', TAG_PATHS.METADATA, TRIAL.tag);
    
else
    
   fprintf('TAG_PATHS.METADATA is not set. Specify path here...\n');
   metaDataPath = input('Metadata path: ','s');
   while ( ~exist(metaDataPath,'dir') == 7 )
       metaDataPath = input('Metadata path: ','s');
       fullFileName = sprintf('%s/%smeta.mat', metaDataPath, TRIAL.tag);
   end
   
end

if ( exist(fullFileName,'file') )
    
   confirmText = 'Metadata file already exists. Overwrite? (y/n) [n]: ';
   confirmStr = lower(input(confirmText, 's'));
   switch(confirmStr)
        case 'y'
           save(fullFileName, 'TRIAL');
        case 'n'
            fprintf('Skipping save to prevent meta data overwrite.\n');
            fprintf('Recommend manually changing data then saving using:\n');
            fprintf('save(fullFileName, ''TRIAL\'')\n');
        otherwise
            fprintf('Skipping save to prevent meta data overwrite.\n');
            fprintf('Recommend manually changing data then saving using:\n');
            fprintf('save(fullFileName, ''TRIAL\'')\n');
   end
   
else        % file does not exist, so save it for the first time...
    
   save(fullFileName, 'TRIAL');
   
end


if ( exist(fullFileName,'file') )
    fprintf('Metadata file is saved and in place!\n');
else
    fprintf('There was a problem finding the metadata file.\n');
end

%% review before exiting

fprintf('\n');
acceptStr = input('Press [return] after reviewing TRIAL meta data. ','s');



%% check to see if there's a meta data file for this tag...

fprintf('Metadata file is specified as: %s\n', metadataFileName);
goodMetaData = lower(input('Is this correct? (y/n): ','s'));
switch(goodMetaData)
    
    case 'y'
        fprintf('Good. Continuing with metadata file load...\n');
        
    otherwise   % repeat until good metadata file name is specified...

        validMetaDataFileName = false;
        
        while(~validMetaDataFileName)
            
            fprintf('You said this is the wrong metadata file.\n');
            correctTrial = input('Enter the correct one: ','s');
            
            % validate correctMetaData structure
            gsCode = string(extractBetween(correctTrial, 1, 2));
            yyCode = str2double(extractBetween(correctTrial, 3, 4));
            usCode = string(extractBetween(correctTrial, 5, 5));
            doyCode = str2double(extractBetween(correctTrial, 6, 8));
            trialCode = string(extractBetween(correctTrial, 9, 9));

            if (  ( gsCode >= 'aa' & gsCode <= 'zz' ) & ...
                (yyCode >= 00 & yyCode <= 99) & ...
                (strcmp(usCode,'_')) ...
                & (doyCode >= 001 & doyCode <= 366) & ...
                (trialCode >= 'a' & trialCode <= 'z') )
                
            % correctMetaData seems valid-ish
                fprintf('New metadata file seems valid: %s\n', ...
                    correctTrial);
                fprintf('Correcting tag trial, META, RAW & PRH names...\n');
                
                validMetaDataFileName = true;
                tag = sprintf('%s', correctTrial);
                rawFileName = sprintf('%sraw.mat', correctTrial);
                prhFileName = sprintf('%sprh.mat', correctTrial);
                metadataFileName = sprintf('%smeta.mat', correctTrial);
                
                fprintf('\nTag ID: %s\n', tag);
                fprintf('RAW filename: %s\n', rawFileName);
                fprintf('PRH filename: %s\n', prhFileName);
                fprintf('METADATA filename: %s\n', metadataFileName);
                
            else
                fprintf('New metadata file seems INVALID: %s\n', correctTrial);
                fprintf('Check the name and try this again...\n');
                validMetaDataFileName = false;
            end
        end
end


fullMetadataFileName = sprintf('%s/%s', TAG_PATHS.METADATA, ...
    metadataFileName);
if ( exist(fullMetadataFileName, 'file') )
    fprintf('Metadata file located for this trial... loading it!\n');
    load(fullMetadataFileName);
    % confirm presence of TRIAL variable
    if ( exist('TRIAL','var') )
        fprintf('TRIAL variable is present and loaded...\n');
        disp(TRIAL);
        disp(TRIAL.SUBJECT);
        disp(TRIAL.CUES);
        fprintf('\n');
    else
        fprintf('TRIAL variable is not present. Investigate.\n');
        return;
    end
else
    fprintf('Metadata file not located for this trial...\n');    
    loadMetaDataStr = input('Do you want to specify one? (y/n) [y] ','s');
    if (isempty(loadMetaDataStr))
        loadMetaDataStr = 'y';
    end
    switch(loadMetaDataStr)
        case 'y'
            metaDataFileText = 'Specify a tag (e.g.: tt21_134a): ';
            metaDataFileStr = input(metaDataFileText,'s');
            if (~isempty(metaDataFileStr))
                metaDataFile = sprintf('%s/%smeta.mat', ...
                    TAG_PATHS.METADATA, metaDataFileStr);
                load(metaDataFile);
                if (exist('TRIAL','var'))
                    fprintf('--------------------------------------------------------\n');
                    disp(TRIAL);
                    fprintf('Loaded a valid TRIAL.\n');
                else
                end
            else
                fprintf('Skipping metadata. You will need to manage this later.\n');
            end
        otherwise
            fprintf('Skipping metadata. You will need to manage this later.\n');
    end
end

%% Load TRIAL into a timetable

lenTrialCuesTime = length(TRIAL.CUES.time);

TRIAL.CUES.duration(1:lenTrialCuesTime) = duration([ 0 0 0 ]);

if (lenTrialCuesTime == 1)
    
    if (isempty(TRIAL.CUES.type))
        TRIAL.CUES.type = {'test'};
    end

    if (isempty(TRIAL.CUES.comment))
        TRIAL.CUES.comment = {'test'};
    end
   
end

% make a replacement array of TRIAL.CUE.time(:) datetimes
% 
% cueDateTimes(1:lenTrialCuesTime) = datetime('now', 'InputFormat', ...
%     'd-MMM-y HH:mm:ss','TimeZone','' );
% 
% for i = 1:lenTrialCuesTime
%     cueDateTimes(i) = datetime(TRIAL.CUES.time(i), 'Format', ...
%         'd-MMM-y HH:mm:ss Z', 'TimeZone', utcStr);
% end

% get rid of the now unused CUE struct and make a timetable called CUE

clear CUE;

CUE = timetable(TRIAL.CUES.time', TRIAL.CUES.duration', ...
    TRIAL.CUES.type', TRIAL.CUES.comment','VariableNames', ...
    {'duration','type','comment' } );     


%% Ask about doing partitioning now for multi-trial tag deployments...

partitionAskStr = lower(input('Partition trial data (y/n) [n]: ','s'));
switch(partitionAskStr)
    case 'y'
        doPartitioning = true;
        useRaw = true;
    otherwise
        doPartitioning = false;
end

%% If doPartitioning is true, do this stuff next...

while (doPartitioning)

    % plot LED2 and LED3
    
    figPart = figure;
    
    plot(optics.Time, optics.led2, optics.Time, optics.led3);
    grid;
    
    % ask user to do two clicks for defining the start region
    startRegionDefined = false;
    
    % ask user to do two clicks to zoom in on the start region
    while(~startRegionDefined)
        
        fprintf('Click once to set the left-edge of the start region.\n');
        [xStartZoom1, ~, ~] = ginput(1);
        xAxis = gca;
        xStart1 = num2ruler(xStartZoom1, xAxis.XAxis);
        
        fprintf('Click once to set the right-edge of the start region.\n');
        [xStartZoom2, ~, ~] = ginput(1);
        xStart2 = num2ruler(xStartZoom2, xAxis.XAxis);

        startL = retime(optics, xStart1, 'nearest');
        startR = retime(optics, xStart2, 'nearest');
        
        startRange = timerange(startL.Time, startR.Time);
    
        figure(figPart);
        plot(optics.Time(startRange), optics.led2(startRange), ...
            optics.Time(startRange), optics.led3(startRange) );
        grid;
        
        % ask user if this is good enough to select an end or do more zooming
        startGoodTxt = 'Is this zoom level good enough? (y/n): ';
        startRegionGood = lower(input(startGoodTxt,'s'));        
        
        switch(startRegionGood)
            case 'y'
                startRegionDefined = true;
            otherwise
                startRegionDefined = false;
        end
        
    end
    
    % click once for trialStart and mvoe on once the user confirms
    trialStartConfirm = false;
    while (~trialStartConfirm)
        
        fprintf ('Click once for the trial start point.\n');
        [xStart,~,~] = ginput(1);
        xAxis = gca;
        partStart = num2ruler(xStart, xAxis.XAxis);
        
        opticsStart = retime(optics, partStart, 'nearest');
        kinematicsStart = retime(kinematics, partStart, 'nearest');
        pressureStart = retime(pressure, partStart, 'nearest');
        magStart = retime(mag, partStart, 'nearest');
        powerStart = retime(power, partStart, 'nearest');
        
        fprintf('oStart: %.1f | kStart: %.1f\n', opticsStart.time_s, ...
            kinematicsStart.time_s);
        
        trialStartConfirmTxt = 'Is this a good trialStart point (y/n)? ';
        trialStartConfirmStr = lower(input(trialStartConfirmTxt, 's'));
        switch(trialStartConfirmStr)
            case 'y'
                trialStartConfirm = true;
            otherwise
                trialStartConfirm = false;
        end
    end
    
    newPartRange = timerange(opticsStart.Time, optics.Time(end));
    
    % replot LED2 and LED3, starting from the now-defined start point
    figure(figPart);
    plot(optics.Time(newPartRange), optics.led2(newPartRange), ...
        optics.Time(newPartRange), optics.led3(newPartRange) );
    grid;
    
    % ask user to do two clicks to zoom in on the end region
    
    endRegionDefined = false;

    while(~endRegionDefined)
        
        fprintf('Click once to set the left-edge of end region selection.\n');
        [xEndZoom1, ~, ~] = ginput(1);
        
        fprintf('Click once to set the right-edge of end region selection.\n');
        [xEndZoom2, ~, ~] = ginput(1);
        
        xAxis = gca;
        xEnd1 = num2ruler(xEndZoom1, xAxis.XAxis);
        xEnd2 = num2ruler(xEndZoom2, xAxis.XAxis);
        
        newEndRange = timerange(xEnd1, xEnd2);
        
        figure(figPart);
        plot(optics.Time(newEndRange), optics.led2(newEndRange), ...
            optics.Time(newEndRange), optics.led3(newEndRange) );
        grid;
        
        % ask user if this is good enough to select an end or do more zooming
        endGoodTxt = 'Is this zoom level good enough? (y/n): ';
        endGood = lower(input(endGoodTxt,'s'));        
        
        switch(endGood)
            case 'y'
                endRegionDefined = true;
            otherwise
                endRegionDefined = false;
        end
        
    end
    
    % once region is good, click for end
    
    fprintf('Click once to select the trial end point.\n');
    [xEnd, ~, ~] = ginput(1);
    xAxis = gca;
    partEnd = num2ruler(xEnd, xAxis.XAxis);
    
    opticsEnd = retime(optics, partEnd, 'nearest');
    kinematicsEnd = retime(kinematics, partEnd, 'nearest');
    pressureEnd = retime(pressure, partEnd, 'nearest');
    magEnd = retime(mag, partEnd, 'nearest');
    powerEnd = retime(power, partEnd, 'nearest');
    
    partitionRange = timerange(opticsStart.Time, opticsEnd.Time);
    
    oRange = timerange(opticsStart.Time,opticsEnd.Time);
    kRange = timerange(kinematicsStart.Time, kinematicsEnd.Time);
    pRange = timerange(pressureStart.Time, pressureEnd.Time);
    mRange = timerange(magStart.Time, magEnd.Time);
    wRange = timerange(powerStart.Time, powerEnd.Time);
    
    fprintf('Generating a plot of the proposed partition.\n');
    
    figure(figPart);
    
    pOk1 = subplot(311);
    plot(optics.Time(partitionRange), optics.led2(partitionRange), ...
        optics.Time(partitionRange), optics.led3(partitionRange) );
    xlim([opticsStart.Time opticsEnd.Time]);
    xlabel('Time, local');
    ylabel('Raw intensities');
    grid;
    legend('1050 nm', '1200 nm');
    
    pOk2 = subplot(312);
    plot(kinematics.Time(partitionRange), kinematics.ax(partitionRange), ...
        kinematics.Time(partitionRange), kinematics.ay(partitionRange), ...
        kinematics.Time(partitionRange), kinematics.az(partitionRange) );
    xlim([kinematicsStart.Time kinematicsEnd.Time]);
    xlabel('Time, seconds');
    ylabel('Ax · Ay · Az');
    grid;
    legend('ax','ay','az');
    
    pOk3 = subplot(313);
    plot(pressure.Time(partitionRange), pressure.depth(partitionRange), ...
        'Color', blue);
    set(gca, 'YDir','reverse');
    xlim([pressureStart.Time, pressureEnd.Time]);
    xlabel('Time, seconds');
    ylabel('Depth, meters');
    grid;
    
    linkaxes([pOk1 pOk2 pOk3],'x');
    
    allGoodTxt = 'Are you completely satisifed with this partitioning (y/n): ';
    allGoodStr = lower(input(allGoodTxt,'s'));
    switch(allGoodStr)
        case 'y'
            doPartitioning = false;
        otherwise
            fprintf('Feeling skittish eh? OK, starting over!\n');
            doPartitioning = true;
    end
    
end

%% Ready to create a new timetable with only the partitioned data?

partTimeTableTxt = 'Ready to create a timetable using only partition data (y/n)? ';
partCreateStr = lower(input(partTimeTableTxt,'s'));
switch(partCreateStr)
    case 'y'
        
        oRange = timerange(opticsStart.Time, opticsEnd.Time);
        kRange = timerange(kinematicsStart.Time, kinematicsEnd.Time);
        pRange = timerange(pressureStart.Time, pressureEnd.Time);
        mRange = timerange(magStart.Time, magEnd.Time);
        wRange = timerange(powerStart.Time, powerEnd.Time);
        
        oResize = optics.Time(oRange);
        kResize = kinematics.Time(kRange);
        pResize = pressure.Time(pRange);
        mResize = mag.Time(mRange);
        wResize = power.Time(wRange);
        
        OPTICS = retime(optics, oResize, 'Timestep', ...
            optics.Properties.TimeStep );
        OPTICS.Properties.TimeStep = optics.Properties.TimeStep;
        OPTICS.Properties.SampleRate = optics.Properties.SampleRate;
        OPTICS.tagOnTime_s = OPTICS.time_s;
        opticsTimeOffset = OPTICS.time_s(1);
        OPTICS.time_s = OPTICS.time_s - opticsTimeOffset;
        
        KINEMATICS = retime(kinematics, kResize);
        KINEMATICS.Properties.TimeStep = kinematics.Properties.TimeStep;
        KINEMATICS.Properties.SampleRate = kinematics.Properties.SampleRate;
        KINEMATICS.tagOnTime_s = KINEMATICS.time_s;
        kinematicsTimeOffset = KINEMATICS.time_s(1);
        KINEMATICS.time_s = KINEMATICS.time_s - kinematicsTimeOffset;
        
        PRESSURE = retime(pressure, pResize);
        PRESSURE.Properties.TimeStep = pressure.Properties.TimeStep;
        PRESSURE.Properties.SampleRate = pressure.Properties.SampleRate;
        PRESSURE.tagOnTime_s = PRESSURE.time_s;
        pressureTimeOffset = PRESSURE.time_s(1);
        PRESSURE.time_s = PRESSURE.time_s - pressureTimeOffset;
        
        MAG = retime(mag, mResize);
        MAG.Properties.TimeStep = mag.Properties.TimeStep;
        MAG.Properties.SampleRate = mag.Properties.SampleRate;
        MAG.tagOnTime_s = MAG.time_s;
        magTimeOffset = MAG.time_s(1);
        MAG.time_s = MAG.time_s - magTimeOffset;
        
        POWER = retime(power, wResize);
        POWER.Properties.TimeStep = power.Properties.TimeStep;
        POWER.Properties.SampleRate = power.Properties.SampleRate;
        POWER.tagOnTime_s = POWER.time_s;
        powerTimeOffset = POWER.time_s(1);
        POWER.time_s = POWER.time_s - powerTimeOffset;
        
    otherwise
        
        fprintf('Leaving entire timetable data in place.\n');
        fprintf('This is not recommended and may impair post-processing.\n');
        return;
        
end

close(figPart);

%% Plot and align / adjust CUE entries

maxCueIndex = height(CUE);

originalFirstCueTime = CUE.Time(1);

figTrial = figure;

pAlign1 = subplot(211);
p1 = plot(KINEMATICS.Time, KINEMATICS.odba, 'Color', blue);
xlim([KINEMATICS.Time(1) KINEMATICS.Time(end)]);
xlabel('Time, local');
ylabel('m/s^2');
hold on;

for i = 1:maxCueIndex
    switch(string(CUE.type(i)))
        case '1'
            plot(CUE.Time(i), max(p1.YData), ...
                'Color', maroon, 'Marker', 'd', 'MarkerFaceColor', maroon);            
        case '2'
            plot(CUE.Time(i), max(p1.YData), ...
                'Color', green, 'Marker', '^', 'MarkerFaceColor', green);
        case '3'
            plot(CUE.Time(i), max(p1.YData), ...
                'Color', red, 'Marker', 'v', 'MarkerFaceColor', red);
        case '4'
            plot(CUE.Time(i), max(p1.YData), ...
                'Color', blue, 'Marker', '^', 'MarkerFaceColor', blue);
        case '*'
            plot(CUE.Time(i), max(p1.YData), ...
                'Color', goldenrod, 'Marker', '*');
            xline(CUE.Time(i), 'Color', goldenrod);
    end
end
hold off;
xlabel('Time, local');
ylabel('ODBA');
grid;

pAlign2 = subplot(212);
p2 = plot(PRESSURE.Time, PRESSURE.depth, 'Color', blue);
set(gca,'YDir','Reverse');
hold on;
for i = 1:maxCueIndex
    switch(string(CUE.type(i)))
        case '1'
            plot(CUE.Time(i), max(p2.YData), ...
                'Color', maroon, 'Marker', 'd', 'MarkerFaceColor', maroon);            
        case '2'
            plot(CUE.Time(i), max(p2.YData), ...
                'Color', green, 'Marker', '^', 'MarkerFaceColor', green);
        case '3'
            plot(CUE.Time(i), max(p2.YData), ...
                'Color', red, 'Marker', 'v', 'MarkerFaceColor', red);
        case '4'
            plot(CUE.Time(i), max(p2.YData), ...
                'Color', blue, 'Marker', '^', 'MarkerFaceColor', blue);
        case '*'
            plot(CUE.Time(i), max(p2.YData), ...
                'Color', goldenrod, 'Marker', '*');
            xline(CUE.Time(i), 'Color', goldenrod);
    end
end
hold off;
xlim([PRESSURE.Time(1), PRESSURE.Time(end)]);
xlabel('Time, local');
ylabel('Depth, meters');
grid;

linkaxes([pAlign1 pAlign2],'x');

durationDifference = OPTICS.Time(1) - datetime(CUE.Time(1), ...
    'TimeZone', utcStr ,'Format', 'dd-MMM-y HH:mm:ss' );
secondsDifference = seconds(durationDifference);
fprintf('OPTICS.Time(1) vs CUE.Time(1) difference: %d seconds\n', ...
    secondsDifference);

coarseAdjustChoice = lower(input('Perform coarse adjustment (y/n) [y]): ','s'));
if (isempty(coarseAdjustChoice))
    coarseAdjustChoice = 'y';
end
switch(coarseAdjustChoice)
    case 'y'
        for i = 1:maxCueIndex
            CUE.Time(i) = CUE.Time(i) + durationDifference;
        end
        fprintf('Completed adjustments to TRIAL.CUE.time(*).\n');
    otherwise
        fprintf('Making no adjustments to TRIAL.CUE.time(*).\n');
end

figure(figTrial);
clf(figTrial);
pAlign1 = subplot(211);
p1 = plot(KINEMATICS.Time, KINEMATICS.odba, 'Color', blue);
xlim([KINEMATICS.Time(1) KINEMATICS.Time(end)]);
hold on;
for i = 1:maxCueIndex
    switch(string(CUE.type(i)))
        case '1'
            plot(CUE.Time(i), max(p1.YData), ...
                'Color', maroon, 'Marker', 'd', 'MarkerFaceColor', maroon);            
        case '2'
            plot(CUE.Time(i), max(p1.YData), ...
                'Color', green, 'Marker', '^', 'MarkerFaceColor', green);
        case '3'
            plot(CUE.Time(i), max(p1.YData), ...
                'Color', red, 'Marker', 'v', 'MarkerFaceColor', red);
        case '4'
            plot(CUE.Time(i), max(p1.YData), ...
                'Color', blue, 'Marker', '^', 'MarkerFaceColor', blue);
        case '*'
            plot(CUE.Time(i), max(p1.YData), ...
                'Color', goldenrod, 'Marker', '*');
            xline(CUE.Time(i), 'Color', goldenrod);
    end
end
hold off;
xlabel('Time, local');
ylabel('ODBA');
grid;

pAlign2 = subplot(212);
p2 = plot(PRESSURE.Time, PRESSURE.depth, 'Color', blue);
set(gca,'YDir','Reverse');
hold on;
for i = 1:maxCueIndex
    switch(string(CUE.type(i)))
        case '1'
            plot(CUE.Time(i), max(p2.YData), ...
                'Color', maroon, 'Marker', 'd', 'MarkerFaceColor', maroon);            
        case '2'
            plot(CUE.Time(i), max(p2.YData), ...
                'Color', green, 'Marker', '^', 'MarkerFaceColor', green);
        case '3'
            plot(CUE.Time(i), max(p2.YData), ...
                'Color', red, 'Marker', 'v', 'MarkerFaceColor', red);
        case '4'
            plot(CUE.Time(i), max(p2.YData), ...
                'Color', blue, 'Marker', '^', 'MarkerFaceColor', blue);
        case '*'
            plot(CUE.Time(i), max(p2.YData), ...
                'Color', goldenrod, 'Marker', '*');
            xline(CUE.Time(i), 'Color', goldenrod);
    end
end
hold off;
xlim([PRESSURE.Time(1), PRESSURE.Time(end)]);
xlabel('Time, local');
ylabel('Depth, meters');
grid;

linkaxes([pAlign1 pAlign2],'x');

cueAlignComplete = false;
updatePlot = false;

while (~cueAlignComplete)
    
    fprintf('[b]ack 5s · [f]orward 5s · [space] +1s · [m]ove +/- seconds · [g]ood · [q]uit\n');
    
    [~,~,button] = ginput(1);
    
    switch(char(button))
        case 'f'
            for i = 1:maxCueIndex
            	CUE.Time(i) = CUE.Time(i) + duration([0 0 5]);
            end
            cueAlignComplete = false
            updatePlot = true;
            
        case 'b'
            for i = 1:maxCueIndex
            	CUE.Time(i) = CUE.Time(i) - duration([0 0 5]);
            end            
            cueAlignComplete = false
            updatePlot = true;
        
        case ' '
            for i = 1:maxCueIndex
            	CUE.Time(i) = CUE.Time(i) + duration([0 0 1]);
            end            
            cueAlignComplete = false
            updatePlot = true;
            
        case 'm'
            secondsUpdateText = 'Enter +/- seconds adjustment: ';
            secondsUpdate = str2double(input(secondsUpdateText,'s'));
            if (isnumeric(secondsUpdate))
                fprintf('Updating by %d seconds...\n', secondsUpdate);
                newDurationOffset = duration([ 0 0 secondsUpdate]);
                for i = 1:maxCueIndex
                    CUE.Time(i) = CUE.Time(i) + newDurationOffset;
                end            
                cueAlignComplete = false;
                updatePlot = true;
            else
                fprintf('That was an invalid number of seconds. Try again.\n');
            end
            
        case 'g'
            
            alignDurationDiff = originalFirstCueTime - CUE.Time(1);
            fprintf('Alignment difference: %d seconds.\n', ...
                seconds(alignDurationDiff) );
            alignSaveText = 'You said aligment looks good. Want to save? (y/n) [y]: ';
            alignSaveStr = lower(input(alignSaveText, 's'));
            switch(alignSaveStr)
                case 'n'
                    fprintf('Exiting without saving. Re-run to save.\n');
                otherwise
                    TRIAL.tagTimelogTimeDifference = alignDurationDiff;
                    TRIAL.CUES.tagTimeOffset = alignDurationDiff;
                    TRIAL.CUES.timeSyncComplete = 'yes';
                    TRIAL.CUES.timeSyncBy = 'Dave Haas';
                    save(fullMetadataFileName, 'TRIAL');
            end
            
            cueAlignComplete = true;
            
        case 'q'
            fprintf('Quitting cue alignment tool.\n');
            cueAlignComplete = true;
            updatePlot = false;
            
        otherwise
            cueAlignComplete = false;
            updatePlot = false;
            
    end
    
    if (updatePlot)

        if (isempty(figTrial.findobj))
            figTrial = figure;
        else
            
            clf(figTrial);
            figure(figTrial);
            
        end
        
        pAlign1 = subplot(211);
        p1 = plot(KINEMATICS.Time, KINEMATICS.odba, 'Color', blue);
        xlim([KINEMATICS.Time(1) KINEMATICS.Time(end)]);
        hold on;
        for i = 1:maxCueIndex
            switch(string(CUE.type(i)))
                case '1'
                    plot(CUE.Time(i), max(p1.YData), ...
                        'Color', maroon, 'Marker', 'd', 'MarkerFaceColor', maroon);            
                case '2'
                    plot(CUE.Time(i), max(p1.YData), ...
                        'Color', green, 'Marker', '^', 'MarkerFaceColor', green);
                case '3'
                    plot(CUE.Time(i), max(p1.YData), ...
                        'Color', red, 'Marker', 'v', 'MarkerFaceColor', red);
                case '4'
                    plot(CUE.Time(i), max(p1.YData), ...
                        'Color', blue, 'Marker', '^', 'MarkerFaceColor', blue);
                case '*'
                    plot(CUE.Time(i), max(p1.YData), ...
                        'Color', goldenrod, 'Marker', '*');
                    xline(CUE.Time(i), 'Color', goldenrod);
            end
        end         % end maxCue iteration
        hold off;   % make sure to do no more plot updating
        xlabel('Time, local');
        ylabel('ODBA');
        grid;
        
        pAlign2 = subplot(212);
        p2 = plot(PRESSURE.Time, PRESSURE.depth, 'Color', blue);
        set(gca,'YDir','Reverse');
        hold on;
        for i = 1:maxCueIndex
            switch(string(CUE.type(i)))
                case '1'
                    plot(CUE.Time(i), max(p2.YData), ...
                        'Color', maroon, 'Marker', 'd', 'MarkerFaceColor', maroon);            
                case '2'
                    plot(CUE.Time(i), max(p2.YData), ...
                        'Color', green, 'Marker', '^', 'MarkerFaceColor', green);
                case '3'
                    plot(CUE.Time(i), max(p2.YData), ...
                        'Color', red, 'Marker', 'v', 'MarkerFaceColor', red);
                case '4'
                    plot(CUE.Time(i), max(p2.YData), ...
                        'Color', blue, 'Marker', '^', 'MarkerFaceColor', blue);
                case '*'
                    plot(CUE.Time(i), max(p2.YData), ...
                        'Color', goldenrod, 'Marker', '*');
                    xline(CUE.Time(i), 'Color', goldenrod);
            end
        end
        hold off;
        xlim([PRESSURE.Time(1), PRESSURE.Time(end)]);
        xlabel('Time, local');
        ylabel('Depth, meters');
        grid;

        linkaxes([pAlign1 pAlign2],'x');

    end             % end updatePlot
end                 % end cue alignment


%% Confirm that these partition-adjusted time series are good

figAlign = figure;

s1 = subplot(411);
p1 = plot(OPTICS.Time, OPTICS.led1, 'Color', black);
hold on;
p2 = plot(OPTICS.Time, OPTICS.led2, 'Color', blue);
p3 = plot(OPTICS.Time, OPTICS.led3, 'Color', red);
p4 = plot(OPTICS.Time, OPTICS.led4, 'Color', goldenrod);
yPos = max([ max(p1.YData) max(p2.YData) max(p3.YData) max(p4.YData) ]);
numCuePlots = ftPlotCues(CUE.Time(:), CUE.type(:), yPos);
hold off;
xlabel('Time, seconds');
ylabel('Reflectance (A.U.)');
xlim([OPTICS.Time(1) OPTICS.Time(end)]);
grid;
legend('amb1','1050 nm','1200 nm','amb2');

s2 = subplot(412);
p5 = plot(KINEMATICS.Time, KINEMATICS.ax, 'Color', blue);
hold on;
p6 = plot(KINEMATICS.Time, KINEMATICS.ay, 'Color', red);
p7 = plot(KINEMATICS.Time, KINEMATICS.az, 'Color', goldenrod);
yPos = max([ max(p5.YData) max(p6.YData) max(p7.YData) ]);
numCuePlots = ftPlotCues(CUE.Time(:), CUE.type(:), yPos);
hold off;
xlabel('Time, seconds');
ylabel('m/s^2');
xlim([KINEMATICS.Time(1) KINEMATICS.Time(end)]);
grid;

s3 = subplot(413);
p8 = plot(KINEMATICS.Time, KINEMATICS.gx, 'Color', blue);
hold on;
p9 = plot(KINEMATICS.Time, KINEMATICS.gy, 'Color', red);
p10 = plot(KINEMATICS.Time, KINEMATICS.gz, 'Color', goldenrod);
yPos = max([ max(p8.YData) max(p9.YData) max(p10.YData) ]);
numCuePlots = ftPlotCues(CUE.Time(:), CUE.type(:), yPos);
hold off;
xlabel('Time, seconds');
ylabel('°/s');
xlim([KINEMATICS.Time(1) KINEMATICS.Time(end)]);
grid;

s4 = subplot(414);
p11 = plot(PRESSURE.Time, PRESSURE.depth, 'Color', blue);
hold on;
yPos = max(p11.YData);
numCuePlots = ftPlotCues(CUE.Time(:), CUE.type(:), yPos);
hold off;
xlabel('Time, seconds');
ylabel('Depth, meters');
set(gca, 'YDir','reverse');
xlim([PRESSURE.Time(1) PRESSURE.Time(end)]);
grid;

linkaxes([s1 s2 s3 s4],'x');

%% Does all that look good

confirmTxt = 'Does everything in this plot look good? (y/n): ';
confirmStr = lower(input(confirmTxt,'s'));
switch(confirmStr)
    case 'y'
        okayToProceed = true;
        
    otherwise
        okayToProceed = false;
end

%% Safe to continue? Is post-processing all done?

proceedTxt = 'Continue with decimation? (y/n): ';
proceedStr = lower(input(proceedTxt,'s'));
switch(proceedStr)
    case 'y'
        fprintf('Proceeding with decimation...\n');
        safeToProceed = true;
    otherwise
        fprintf('You said no, so stopping execution here. Rerun when ready.\n');
        return;
end


%% Decimate optical (250 Hz) + accel, gyro, depth & temp (100 Hz) to 25 Hz

fprintf('Decimating optical, movement, depth & temp data to 25Hz...\n');

oT = downsample(OPTICS.Time, 10);
oT_s = downsample(OPTICS.time_s, 10);
led1_25Hz = decimate(OPTICS.led1, 10);         % decimate 250 Hz by a factor of 10
led2_25Hz = decimate(OPTICS.led2, 10);
led3_25Hz = decimate(OPTICS.led3, 10);
led4_25Hz = decimate(OPTICS.led4, 10);

kT = downsample(KINEMATICS.Time, 4);
kT_s = downsample(KINEMATICS.time_s, 4);
ax_25Hz = decimate(KINEMATICS.ax, 4);              % decimate 100 Hz by a factor of 4
ay_25Hz = decimate(KINEMATICS.ay, 4);
az_25Hz = decimate(KINEMATICS.az, 4);

gx_25Hz = decimate(KINEMATICS.gx, 4);
gy_25Hz = decimate(KINEMATICS.gy, 4);
gz_25Hz = decimate(KINEMATICS.gz, 4);

pitch_25Hz   = decimate(KINEMATICS.pitch, 4);
roll_25Hz    = decimate(KINEMATICS.roll, 4);
heading_25Hz = decimate(KINEMATICS.heading, 4);

odba_25Hz = decimate(KINEMATICS.odba, 4);
odav_25Hz = decimate(KINEMATICS.odav, 4);

depth_25Hz          = decimate(PRESSURE.depth, 4);
temperature_25Hz    = decimate(PRESSURE.temperature, 4);

% confirm lengths of optics vs. kinematics

if (length(odba_25Hz) == length(led1_25Hz) )
    
    fprintf('Decimated optics and kinematics have matching lengths.\n');

else
    
    fprintf('Decimated optics and kinematics have differing lengths.\n');
    dOpticsLength = length(led1_25Hz);
    dKinematicsLength = length(odba_25Hz);
    fprintf('Optics length: %d | Kinematics length: %d\n', ...
        dOpticsLength, dKinematicsLength);
    
    figDec = figure;
    plot(oT, normalize(led2_25Hz), 'Color', blue);
    hold on;
    plot(oT, normalize(led3_25Hz), 'Color', red);
    plot(kT, normalize(odba_25Hz), 'Color', goldenrod);
    grid;
    legend('1050 nm','1200 nm', 'ODBA');
    
    fprintf('Exam this figure closely before approving truncation.\n');
    
    if (dOpticsLength > dKinematicsLength)
        
        dOffset = dOpticsLength - dKinematicsLength;
        fprintf('Proposal: truncate the last %d optical sample(s)\n', ...
            dOffset);
        truncateStr = lower(input('Okay to truncate (y/n): ','s'));
        switch(truncateStr)
            case 'y'
            	oT = oT(1:end-dOffset);
                oT_s = oT_s(1:end-dOffset);
                led1_25Hz = led1_25Hz(1:end-dOffset);
                led2_25Hz = led2_25Hz(1:end-dOffset);
                led3_25Hz = led3_25Hz(1:end-dOffset);
                led4_25Hz = led4_25Hz(1:end-dOffset);
            otherwise
                fprintf('Leaving optics with that extra %d sample(s).\n', ...
                    dOffset);
        end
        
        
    else
        
        dOffset = dKinematicsLength - dOpticsLength;

        fprintf('Proposal: truncate the last %d kinematic sample(s)\n', ...
            dOffset);
       truncateStr = lower(input('Okay to truncate (y/n): ','s'));
        switch(truncateStr)
            case 'y'
                kT = kT(1:end-dOffset);
                kT_s = kT_s(1:end-dOffset);
                ax_25Hz = ax_25Hz(1:end-dOffset);
                ay_25Hz = ay_25Hz(1:end-dOffset);
                az_25Hz = az_25Hz(1:end-dOffset);
                gx_25Hz = gx_25Hz(1:end-dOffset);
                gy_25Hz = gy_25Hz(1:end-dOffset);
                gz_25Hz = gz_25Hz(1:end-dOffset);
                pitch_25Hz = pitch_25Hz(1:end-dOffset);
                roll_25Hz = roll_25Hz(1:end-dOffset);
                heading_25Hz = heading_25Hz(1:end-dOffset);
                odba_25Hz = odba_25Hz(1:end-dOffset);
                odav_25Hz = odav_25Hz(1:end-dOffset);
                depth_25Hz = depth_25Hz(1:end-dOffset);
                temperature_25Hz = temperature_25Hz(1:end-dOffset);
            otherwise
                fprintf('Leaving optics with that extra %d sample(s).\n', ...
                    dOffset);
        end
        
    end         % end differing optics length from kinematics length
    
    if (exist('figDec','var'))
        close(figDec);
    end
    
end             % end decimated time series difference tests



%% Build timetable for 25 Hz PRH file

O = timetable(oT, oT_s, led1_25Hz, led2_25Hz, led3_25Hz, led4_25Hz, ...
    'VariableNames', {'time_s','led1','led2','led3','led4'} );

K = timetable(kT, kT_s, ax_25Hz, ay_25Hz, az_25Hz, gx_25Hz, ...
    gy_25Hz, gz_25Hz, pitch_25Hz, roll_25Hz, heading_25Hz, ...
    odba_25Hz, odav_25Hz, 'VariableNames', ...
    {'time_s','ax','ay','az','gx','gy','gz','pitch', ...
    'roll', 'heading', 'odba', 'odav'} );

P = timetable(kT, kT_s, depth_25Hz, temperature_25Hz, ...
    'VariableNames', {'time_s','depth','temperature' } );


%% Make sure that 'tag' and TRIAL.tag are the same, 
%   and have the user fix if not

if (strcmp(tag, TRIAL.tag))     % tag and TRIAL.tag match
    fprintf('tag and TRIAL.tag match. This is good news!\n');
% else                            % tag and TRIAL.tag do not match
%     fprintf('CUE.Props...tag and tag variable do not match. Resolve now...\n');
%     fprintf('\tCUE.Props...tag = %s\n', CUE.Properties.CustomProperties.tag);
%     fprintf('\ttag = %s\n', tag);
%     tagIdResolveText = 'Do you want to set `tag` equal to CUE.Props...tag? (y/n) [y]: ';
%     tagIdResolveStr = lower(input(tagIdResolveText,'s'));
%     if (isempty(tagIdResolveStr))
%         tagIdResolveStr = 'y';
%     end
%     switch (tagIdResolveStr)
%         case 'y'
%             tag = CUE.Properties.CustomProperties.tag;
%             fprintf('Renamed `tag`: %s\n', tag);
%             rawFileName = sprintf('%sraw.mat', tag);
%             fprintf('Renamed RAW file: %s\n', rawFileName);
%             prhFileName = sprintf('%sprh.mat', tag);
%             fprintf('Renamed PRH file: %s\n', prhFileName);
%         otherwise
%             fprintf('Taking no action and going with aggregated RAW & PRH files.\n');
%     end 
end


%% save all RAW data structs as a somewhat Dtag-compatible raw file

global TAG_PATHS;

rawFile = sprintf('%s/%s', TAG_PATHS.RAW, rawFileName);

fprintf('Creating RAW file: %s\n', rawFile);

if (exist(rawFile, 'file'))
    
    fprintf('RAW file exists... ');
    
    prompt = 'over-write and re-create? y/n [n]: ';
    str = input(prompt,'s');
    if isempty(str)
        str = 'n';
    end
    
    if (strcmp(str, 'y'))
       fprintf('\nOver-writing and re-creating RAW file...\n');
       save(rawFile, 'OPTICS', 'KINEMATICS', 'PRESSURE', ...
            'MAG', 'POWER', 'TRIAL', 'CUE');
    end    
    
else
    
    fprintf('Raw file does not exist... creating...\n');
    save(rawFile, 'OPTICS', 'KINEMATICS', 'PRESSURE', ...
        'MAG', 'POWER', 'TRIAL', 'CUE');
    if (exist(rawFile, 'file'))
        fprintf('Looks like that file was created!\n');
    else
        fprintf('There was some problem saving that raw file. Halting.\n');
    end
    
end


%% save decimated and aligned data to the PRH file

prhFile = sprintf('%s/%s', TAG_PATHS.PRH, prhFileName);

fprintf('Creating PRH file: %s\n', prhFile);

if (exist(prhFile, 'file'))
    
    fprintf('PRH file exists... ');
    
    prompt = 'over-write and re-create? y/n [n]: ';
    str = input(prompt,'s');
    if isempty(str)
        str = 'n';
    end
    
    if (strcmp(str, 'y'))
       fprintf('\nOver-writing and re-creating PRH file...\n');
       save(prhFile, 'O', 'K', 'P', 'TRIAL', 'CUE');
    end
    
else
    fprintf('PRH file does not exist... creating...\n');
    save(prhFile, 'O', 'K', 'P', 'TRIAL', 'CUE');
    if (exist(prhFile, 'file'))
        fprintf('Looks like that file was created!\n');
    else
        fprintf('There was some problem saving that raw file. Halting.\n');
    end
end

%% Confirm RAW and PRH files are where they should be...

if (exist(rawFile, 'file'))
    fprintf('Good news: RAW file is where it should be!\n');
else
    fprintf('Bad news: RAW file is missing from PRH folder. Investigate.\n');
end

if (exist(prhFile, 'file'))
    fprintf('Good news: PRH file is where it should be!\n');
else
    fprintf('Bad news: PRH file is missing from PRH folder. Investigate.\n');
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

fprintf('Basic RAW and PRH files are created.\n');


%% Done!

fprintf('Post-processing complete! Trial is: %s\n', tag);


%% %% Do a bandpass filter of 250 Hz for O2 sat

% %   These are filter options for relative O2 tissue saturation
% 
% filterStr = 'High-pass (h) | Low-pass (l) | Band-pass (b) [b]: ';
% validFilter = false;
% 
% while (~validFilter)
%     
%     o2FilterOptionStr = lower(input(filterStr,'s'));
% 
%     if (isempty(o2FilterOptionStr))
%         o2FilterOptionStr = 'b';
%     end
%     
%     switch(o2FilterOptionStr)
%         
%         case 'l'
%             validFilter = false;
%             fprintf('For now, there is no low-pass filter option.\n');
% %             bwFiltOrder     = 3;        % FieldTrip suggested default
% %             passbandFreqHz  = 0.8;     % FieldTrip suggested default
% %             nyquistFs       = 250;       % should probably experiment with range
% % 
% %             fprintf('Filter order: %d\n', bwFiltOrder);
% %             fprintf('Passband frequency: %d Hz\n', passbandFreqHz);
% %             fprintf('Nyquist frequency: %d\n', nyquistFs);
% % 
% % 
% %             % filter coefficients for low-pass
% %             [Alow, Blow] = butter(bwFiltOrder, ( passbandFreqHz / (nyquistFs / 2) ) );
% % 
% %             ODLed1 = filter(Blow, Alow, odLed1);    
%         
%         case 'h'
%             validFilter = true;
%             bwFiltOrder     = 3;        % FieldTrip suggested default
%             passbandFreqHz  = 0.01;     % FieldTrip suggested default
%             nyquistFs       = 125;       % should probably experiment with range
% 
%             fprintf('Filter order: %d\n', bwFiltOrder);
%             fprintf('Passband frequency: %d Hz\n', passbandFreqHz);
%             fprintf('Nyquist frequency: %d\n', nyquistFs);
% 
%             % filter coefficients for high-pass
%             %[Ao2, Bo2, Co2, Do2] = butter(bwFiltOrder/2, ( passbandFreqHz / (nyquistFs/2) ) );
%             [Ao2, Bo2, Co2, Do2] = butter(bwFiltOrder, ( passbandFreqHz / (nyquistFs) ) );
% 
%             % design the high-pass filter
%             o2BwFilt = designfilt('highpassiir','FilterOrder',bwFiltOrder, ...
%                 'PassbandFrequency',passbandFreqHz, ... 
%                 'PassbandRipple', 0.2, ...
%                 'SampleRate', nyquistFs);
%             
%         case 'b'
%             
%             validFilter = true;
%             bwFiltOrder = 4;        % 20th order butterworth filter
%             lowCut      = 0.01;       % 0.01 Hz   = fluctuations over 10 seconds
%             highCut     = 0.1;        % 0.1 Hz    = fluctuations over 100 seconds
%             nyquistFs   = 25;        
% 
%             fprintf('Applying a band-pass filter with the following properties:\n');
%             fprintf('Filter order: %d\n', bwFiltOrder);
%             fprintf('Low cut frequency: %d Hz\n', lowCut );
%             fprintf('High cut frequency: %d Hz\n', highCut );
%             fprintf('Nyquist frequency: %d\n', nyquistFs);
% 
%             % band-pass
%             [Ao2, Bo2, Co2, Do2] = butter(bwFiltOrder/2, [lowCut highCut] / (nyquistFs/2) );
% 
%             % design the band-pass filter
%             o2BwFilt = designfilt('bandpassiir','FilterOrder',bwFiltOrder, ...
%                 'HalfPowerFrequency1',lowCut,'HalfPowerFrequency2',highCut, ...
%                 'SampleRate', nyquistFs);   
%             
%         otherwise
%             validFilter = false;
%             
%     end
% 
% end
% 
% % Convert the state-space representation to second-order sections. 
% % Visualize the frequency responses using fvtool.
% 
% sos = ss2sos(Ao2,Bo2,Co2,Do2);
% fvt = fvtool(sos, o2BwFilt, 'Fs', nyquistFs);
% legend(fvt,'butter','designfilt');
% 
% approveFilterStr = lower(input('Is this filter good to go? (y/n) [n]: ','s'));
% 
% switch(approveFilterStr)
% 
%     case 'y'
%         
%         % filter optical data, using d' to row-major output
% 
%         fprintf('Using this filter to for TSI / O2 optical data...\n');
%         
%         o2filtLed1 = filtfilt(o2BwFilt', led1);
%         o2filtLed2 = filtfilt(o2BwFilt', led2);
%         o2filtLed3 = filtfilt(o2BwFilt', led3);
%         o2filtLed4 = filtfilt(o2BwFilt', led4);
%         
%         o2FiltLed1 = rescale(o2filtLed1, 0, max(o2filtLed1));
%         o2FiltLed2 = rescale(o2filtLed2, 0, max(o2filtLed2));
%         o2FiltLed3 = rescale(o2filtLed3, 0, max(o2filtLed3));
%         o2FiltLed4 = rescale(o2filtLed4, 0, max(o2filtLed4));
%         
%         odLed1 = log( max(o2FiltLed1) ./ o2FiltLed1 );
%         odLed2 = log( max(o2FiltLed2) ./ o2FiltLed2 );
%         odLed3 = log( max(o2FiltLed3) ./ o2FiltLed3 );
%         odLed4 = log( max(o2FiltLed4) ./ o2FiltLed4 );
% 
%         figure;
%         plot(opticsTime_s, odLed1, 'k-');
%         hold on;
%         plot(opticsTime_s, odLed2, 'b-');
%         plot(opticsTime_s, odLed3, 'r-');
%         plot(opticsTime_s, odLed4, 'g-');
%         hold off;
%         grid;
%         xlabel('Time, seconds');
%         ylabel('Optical density');
%         titleStr = ... 
%             sprintf('OD (0.01 Hz HPF w/ rescale) - RAW %s', TRIAL.tag);
%         title(titleStr);
%         legend('Amb1','1050 nm','1200 nm','amb2');
%         
%         odLed1_25Hz = downsample(odLed1, 10);
%         odLed2_25Hz = downsample(odLed2, 10);
%         odLed3_25Hz = downsample(odLed3, 10);
%         odLed4_25Hz = downsample(odLed4, 10);
%         
%         figure;
%         plot(time_s, odLed1_25Hz, 'k-');
%         hold on;
%         plot(time_s, odLed2_25Hz, 'b-');
%         plot(time_s, odLed3_25Hz, 'r-');
%         plot(time_s, odLed4_25Hz, 'g-');
%         hold off;
%         grid;
%         xlabel('Time, seconds');
%         ylabel('Optical density');
%         titleStr = ... 
%             sprintf('OD (0.01 Hz HPF w/ rescale) - PRH %s', TRIAL.tag);
%         title(titleStr);
%         legend('Amb1','1050 nm','1200 nm','amb2');
% 
%         validFilter = true;
% 
%     otherwise
%         
%         fprintf('Re-run this optical filtering for O2.\n');
%         return;
% 
% end
% 
% 
