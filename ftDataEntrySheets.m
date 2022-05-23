%% FaunaTag Data Sheet Input for Matlab
%
%   Written by Dave Haas, 16 June 2021
%
%   This is an app to manually enter time cues and behavioral observations
%   made during animal field trials involving the FaunaTag.
%   
%   To make this app maximally effective, the time of activation should be
%   recorded on the data sheet to the nearest possible second. In cases
%   where the exact time can't be perfectly known, e.g.: in wild field
%   settings with pole attachments, to the nearest minute is probably
%   sufficient. The time difference between observed tag activation will be
%   offset by the time of activation recorded by the FaunaTag. This offset
%   can then be applied to create time-synchronized audit files of
%   releveant behavioral or observation comment cues.

clear;
clc;

global TAG_PATHS;

if ( ~exist(TAG_PATHS.METADATA, 'var') )
end

%% Enter animal details from this trial

fprintf('Enter trial details...\n');

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

%% 

fprintf('Body location of tag attachment, choose from one of these:\n');
fprintf('\t1: post-nuchal\n');
fprintf('\t2: anterio-dorsal\n');
fprintf('\t3R: dorso-lateral (right)\n');
fprintf('\t3L: dorso-lateral (left)\n');
fprintf('\t4: ventral chest\n');
fprintf('\t5R: right saddle-patch');
fprintf('\t5L: left saddle-patch');

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

%% Tag details

tagIdStr = input('FaunaTag ID (110, 111, or 112) [111]: ', 's');

switch (tagIdStr)
    case ''
        tagId = '111';
    case '110'
        tagId = 110;
    case '111'
        tagId = 111;
    case '112'
        tagId = 112;
    otherwise
        tagId = 111;
end


%% Date and time 

year = 0;

while (year == 0)
    yearStr = input('Trial Year [2021]: ','s');
    if (isempty(yearStr))
        year = 2021;
    else
        year = str2double(yearStr);
        if (isnan(year))
            year = 0;
        end
    end
end

month = 0;

while (month == 0)
    monthStr = input('Trial Month (e.g.: 5 for May) [5]: ','s');
    if (isempty(monthStr))
        month = 5;
    else
        month = str2double(monthStr);
        if (isnan(month))
            month = 0;
        end
    end
end

day = 0;
while (day == 0)
    dayStr = input('Trial Day: ','s');
    if (isempty(dayStr))
        year = 0;
    else
        day = str2double(dayStr);
        if ( isnan(day) || day > 31 || day < 0 )
            day = 0;
        end
    end
end

tagTimeZoneStr = input('UTC Offset in hours (e.g.: Hawaii: -10) [-10]: ','s');
if (isempty(tagTimeZoneStr))
   tagTimeZoneStr = '-10';
end
tz = str2double(tagTimeZoneStr);

tagDate = datetime( [year, month, day], 'InputFormat','yyyy-MM-dd Z', ...
    'TimeZone',tagTimeZoneStr);


% tag activation datetime

cueTimeText = 'Tag activation time (hh:mm:ss): ';
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
tagActivateDateTime = datetime([year month day hms], ...
    'TimeZone',tagTimeZoneStr,'Format','d-MMM-y HH:mm:ss Z');
fprintf('Tag activation datetime: %s\n', tagActivateDateTime);

% tag on datetime

cueTimeText = 'Tag ON time (hh:mm:ss): ';
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
tagOnDateTime = datetime([year month day hms], ...
    'TimeZone',tagTimeZoneStr,'Format','d-MMM-y HH:mm:ss Z');
fprintf('Tag on datetime: %s\n', tagOnDateTime);

% tag off datetime
cueTimeText = 'Tag OFF time (hh:mm:ss): ';
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
tagOffDateTime = datetime([year month day hms], ...
    'TimeZone',tagTimeZoneStr,'Format','d-MMM-y HH:mm:ss Z');
fprintf('Tag off datetime: %s\n', tagOffDateTime);

% tag deactivate datetime

cueTimeText = 'Tag deactivation time (hh:mm:ss): ';
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
tagDeactivateDateTime = datetime([year month day hms], ...
    'TimeZone',tagTimeZoneStr,'Format','d-MMM-y HH:mm:ss Z');
fprintf('Tag deactivation datetime: %s\n', tagDeactivateDateTime);

%% Build additional tag meta data from this trial

tagOrdId = lower(input('Trial of the day letter (e.g.: 2nd: b; 5th: e): ','s'));

yearId = tagOnDateTime.Year-2000;
doyId = jday([tagOnDateTime.Year tagOnDateTime.Month tagOnDateTime.Day]);

tag = sprintf('%s%d_%d%s', gsCode, yearId, doyId, tagOrdId);

fprintf('Tag ID is %s\n', tag);

%% get data fileName, do some comparisons, and improve metadata

fileName = '';
while (isempty(fileName))
    fileName = input('FaunaData file (do not leave blank): ', 's');
    validNameTest1 = extractBetween(fileName, 6, 6);
    validNameTest2 = extractBetween(fileName, 15, 15);
    if ( (~strcmp(validNameTest1, '_')) || (~strcmp(validNameTest2, '_')) )
        fprintf('Invalid file name. Check it and re-enter it.\n');
        fileName = '';
    end
    
end

fprintf('Parsing fileName and comparing with entered datetimes...\n');

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
trialDateTimeLocal = datetime(trialDateTimeLocal, ...
    'Format','yyyy-MM-dd HH:mm:ss Z','TimeZone','-10');

trialDateTimeUTC   = datetime(YYYY, MM, DD, hh + (utc * -1), mm, ss);
trialDateTimeUTC = datetime(trialDateTimeUTC, ...
    'Format','yyyy-MM-dd HH:mm:ss Z','TimeZone','+00');

if (J == doyId && utc == tz && ...
        ( datenum(trialDate) == datenum(tagDate) ) )  

    if (YYYY == 2021 && MM == 5)

        % these are all at DQO, so auto-fill this info...
        
        lat = 21.272037;            % decimal degrees at Dolphin Quest Oahu
        long = -157.773116;         % decimal degrees at Dolphin Quest Oahu
        decl = 9.48;                % units: degrees, on 13 May 2021
        magFieldStrength = 34797.5; % field strength, units: nanoTeslas (nT)

    else

        % ... but enter this manually if not DQO trials ...

        validMetaData = false;

        while(validMetaData == false)

            fprintf('\n----- Enter additional tag meta data manually... -----\n\n');

            latitude = input('Enter tagOn latitude in decimal degrees (S is negative): ', 's');
            longitude = input('Enter tagOn longitude in decimal degrees (W is negative): ', 's');
            declAngle = input('Enter declination angle in degrees: ', 's');
            fieldStrength = input('Enter magnetic field strength in nT: ', 's');

            lat = str2double(latitude);
            long = str2double(longitude);
            decl = str2double(declAngle);
            magFieldStrength = str2double(fieldStrength);

            % verify the condition of user entries...

            if ( 	isnan(lat)|| isnan(long) || ... 
                    isnan(decl) || isnan(magFieldStrength ) )
                validMetaData = false;
                fprintf('\nTag meta data entered is invalid... Try re-entering it...\n\n');
            else
                validMetaData = true;
                fprintf('Tag meta data seems valid... continuing!\n\n');
            end

        end
    end
end 

%% Enter trial participants and context notes from this trial...

trialPIStr = input('Enter the name of the trial PI [Dave Haas]: ','s');
if (isempty(trialPIStr))
    trialPIStr = 'Dave Haas';
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

TRIAL = struct;
TRIAL.tag = tag;
TRIAL.tagNumber = tagId;
TRIAL.attachmentLocation = attachId;
TRIAL.attachmentSide = attachSideId;
TRIAL.dataFile = fileName;

TRIAL.tagStartTimeLocal = trialDateTimeLocal;
TRIAL.tagStartTimeUTC = trialDateTimeUTC;
TRIAL.utcOffset = tz;
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

doTimeCueEntry  = lower(input('Any detailed time cues to enter (y/n)? ','s'));

if (strcmp(doTimeCueEntry,'y'))
   
    fprintf('\nEnter detailed time cues for this trial...\n'); 

    timeCueEntryComplete = false;
    
    if ( exist('CUE','var') )   % if CUEs entry accidentally exited...
    
        cueIndex = numel(CUE.time); % pick up where it left off...
        
    else
        
        cueIndex = 0;

        CUE = struct;
        CUE.time = datetime();
        CUE.type = cell('');
        CUE.comment = cell('');
        CUE.tagTimeOffset = duration([0 0 0]);
        CUE.investigator = 'Dave Haas';
        CUE.timingDevice = 'Garmin Fenix 5x';
    
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
        CUE.time(cueIndex) = datetime([year month day]) + duration(hms);
        
        
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
    CUE.investigator = 'Dave Haas';
    CUE.timingDevice = 'Garmin Fenix 5x';

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


%% 

global TAG_PATHS;

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

