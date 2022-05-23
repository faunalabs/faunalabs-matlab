%% ftCardioEnergyStats.m
%
%   Written by Dave Haas on 18 September

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
grey = [0.6 0.6 0.6];               % '#999999'
black = [0 0 0];                    % '#000000
white = [1 1 1];                    % '#ffffff'

% set some default LED* colors...

led1Color = black;
led2Color = blue;
led3Color = red;
led4Color = goldenrod;

% ... and some default kinematic colors...

kxColor = blue;
kyColor = red;
kzColor = goldenrod;


%% 

fullDataFileName = '/Users/dave/Documents/MATLAB/tursiopsTimetables.mat';

if (exist(fullDataFileName, 'file'))
    
    getOk = lower(input('File exists, overwrite (y/n) [y]: ','s'));
    switch(getOk)
        case 'n'
            fprintf('Okay. Bailing!\n');
            goAhead = false;
            return;
        otherwise
            fprintf('Cool. Processing!\n');
            load(fullDataFileName);
            goAhead = true;
    end
    
else
    
    goAhead = true;     % file doesn't exist, so proceed!

end
   
if (goAhead)
    
    % aggregate all the vent data into a single timetable
    
    % HOKU
    
    load('/Users/dave/Desktop/DTAG/raw/tt21_128fraw.mat');

    HOKU_E1_TIME = E1_KINEMATICS.Time;
    HOKU_E1_WSST_HR = E1_KINEMATICS{:,22:24};
    HOKU_E1_EPOCH(1:height(E1_KINEMATICS),1) = 1;
    HOKU_E1_MASS(1:height(E1_KINEMATICS),1) = 179.6;
    HOKU_E1_ID(1:height(E1_KINEMATICS),1) = {'63HF'};
    HOKU_E1_ALTID(1:height(E1_KINEMATICS),1) = 1;
    HOKU_E1_VENT_DATA  = E1_VENT_DATA;
    HOKU_E1_CC_DATA    = E1_CC_ENERGY;
    
    HOKU_E3_TIME = E3_KINEMATICS.Time;
    HOKU_E3_WSST_HR = E3_KINEMATICS{:,22:24};
    HOKU_E3_EPOCH(1:height(E3_KINEMATICS),1) = 3;
    HOKU_E3_MASS(1:height(E3_KINEMATICS),1) = 179.6;
    HOKU_E3_ID(1:height(E3_KINEMATICS),1) = {'63HF'};
    HOKU_E3_ALTID(1:height(E3_KINEMATICS),1) = 1;
    HOKU_E3_VENT_DATA  = E3_VENT_DATA;
    HOKU_E3_CC_DATA    = E3_CC_ENERGY;
    
    % no E5 - tag came off during fast swim transition
    
    % LIHO
    
    load('/Users/dave/Desktop/DTAG/raw/tt21_134araw.mat');
    
    LIHO_E1_TIME = E1_KINEMATICS.Time;
    LIHO_E1_WSST_HR = E1_KINEMATICS{:,22:24};
    LIHO_E1_EPOCH(1:height(E1_KINEMATICS),1) = 1;
    LIHO_E1_MASS(1:height(E1_KINEMATICS),1) = 163.7;
    LIHO_E1_ID(1:height(E1_KINEMATICS),1) = {'01L5'};
    LIHO_E1_ALTID(1:height(E1_KINEMATICS),1) = 2;
    LIHO_E1_VENT_DATA  = E1_VENT_DATA;
    LIHO_E1_CC_DATA    = E1_CC_ENERGY;

    LIHO_E3_TIME = E3_KINEMATICS.Time;
    LIHO_E3_WSST_HR = E3_KINEMATICS{:,22:24};
    LIHO_E3_EPOCH(1:height(E3_KINEMATICS),1) = 3;
    LIHO_E3_MASS(1:height(E3_KINEMATICS),1) = 163.7;
    LIHO_E3_ID(1:height(E3_KINEMATICS),1) = {'01L5'};
    LIHO_E3_ALTID(1:height(E3_KINEMATICS),1) = 2;
    LIHO_E3_VENT_DATA  = E3_VENT_DATA;
    LIHO_E3_CC_DATA    = E3_CC_ENERGY;
    
    LIHO_E5_TIME = E5_KINEMATICS.Time;
    LIHO_E5_WSST_HR = E5_KINEMATICS{:,22:24};
    LIHO_E5_EPOCH(1:height(E5_KINEMATICS),1) = 5;
    LIHO_E5_MASS(1:height(E5_KINEMATICS),1) = 163.7;
    LIHO_E5_ID(1:height(E5_KINEMATICS),1) = {'01L5'};
    LIHO_E5_ALTID(1:height(E5_KINEMATICS),1) = 2;
    LIHO_E5_VENT_DATA  = E5_VENT_DATA;
    LIHO_E5_CC_DATA    = E5_CC_ENERGY;
    
    % KOLOHE
    
    load('/Users/dave/Desktop/DTAG/raw/tt21_134braw.mat');
    
    KOLOHE_E1_TIME = E1_KINEMATICS.Time;
    KOLOHE_E1_WSST_HR = E1_KINEMATICS{:,22:24};
    KOLOHE_E1_EPOCH(1:height(E1_KINEMATICS),1) = 1;
    KOLOHE_E1_MASS(1:height(E1_KINEMATICS),1) = 200.9;
    KOLOHE_E1_ID(1:height(E1_KINEMATICS),1) = {'6JK5'};
    KOLOHE_E1_ALTID(1:height(E1_KINEMATICS),1) = 3;
    KOLOHE_E1_VENT_DATA  = E1_VENT_DATA;
    KOLOHE_E1_CC_DATA    = E1_CC_ENERGY;

    KOLOHE_E3_TIME = E3_KINEMATICS.Time;
    KOLOHE_E3_WSST_HR = E3_KINEMATICS{:,22:24};
    KOLOHE_E3_EPOCH(1:height(E3_KINEMATICS),1) = 3;
    KOLOHE_E3_MASS(1:height(E3_KINEMATICS),1) = 200.9;
    KOLOHE_E3_ID(1:height(E3_KINEMATICS),1) = {'6JK5'};
    KOLOHE_E3_ALTID(1:height(E3_KINEMATICS),1) = 3;
    KOLOHE_E3_VENT_DATA  = E3_VENT_DATA;
    KOLOHE_E3_CC_DATA    = E3_CC_ENERGY;
    
    
    KOLOHE_E5_TIME = E5_KINEMATICS.Time;
    KOLOHE_E5_WSST_HR = E5_KINEMATICS{:,22:24};
    KOLOHE_E5_EPOCH(1:height(E5_KINEMATICS),1) = 5;
    KOLOHE_E5_MASS(1:height(E5_KINEMATICS),1) = 200.9;
    KOLOHE_E5_ID(1:height(E5_KINEMATICS),1) = {'6JK5'};
    KOLOHE_E5_ALTID(1:height(E5_KINEMATICS),1) = 3;
    KOLOHE_E5_VENT_DATA  = E5_VENT_DATA;
    KOLOHE_E5_CC_DATA    = E5_CC_ENERGY;
    
    % HUA
    
    load('/Users/dave/Desktop/DTAG/raw/tt21_141eraw.mat');

    HUA_E1_TIME = E1_KINEMATICS.Time;
    HUA_E1_WSST_HR = E1_KINEMATICS{:,22:24};
    HUA_E1_EPOCH(1:height(E1_KINEMATICS),1) = 1;
    HUA_E1_MASS(1:height(E1_KINEMATICS),1) = 147.0;
    HUA_E1_ID(1:height(E1_KINEMATICS),1) = {'83H1'};
    HUA_E1_ALTID(1:height(E1_KINEMATICS),1) = 4;
    HUA_E1_VENT_DATA  = E1_VENT_DATA;
    HUA_E1_CC_DATA    = E1_CC_ENERGY;

    HUA_E3_TIME = E3_KINEMATICS.Time;
    HUA_E3_WSST_HR = E3_KINEMATICS{:,22:24};
    HUA_E3_EPOCH(1:height(E3_KINEMATICS),1) = 3;
    HUA_E3_MASS(1:height(E3_KINEMATICS),1) = 147.0;
    HUA_E3_ID(1:height(E3_KINEMATICS),1) = {'83H1'};
    HUA_E3_ALTID(1:height(E3_KINEMATICS),1) = 4;
    HUA_E3_VENT_DATA  = E3_VENT_DATA;
    HUA_E3_CC_DATA    = E3_CC_ENERGY;
    
    HUA_E5_TIME = E5_KINEMATICS.Time;
    HUA_E5_WSST_HR = E5_KINEMATICS{:,22:24};
    HUA_E5_EPOCH(1:height(E5_KINEMATICS),1) = 5;
    HUA_E5_MASS(1:height(E5_KINEMATICS),1) = 147.0;
    HUA_E5_ID(1:height(E5_KINEMATICS),1) = {'83H1'};
    HUA_E5_ALTID(1:height(E5_KINEMATICS),1) = 4;
    HUA_E5_VENT_DATA  = E5_VENT_DATA;
    HUA_E5_CC_DATA    = E5_CC_ENERGY;    

    % NOA
    
    load('/Users/dave/Desktop/DTAG/raw/tt21_142craw.mat');
    
    NOA_E1_TIME = E1_KINEMATICS.Time;
    NOA_E1_WSST_HR = E1_KINEMATICS{:,22:24};
    NOA_E1_EPOCH(1:height(E1_KINEMATICS),1) = 1;
    NOA_E1_MASS(1:height(E1_KINEMATICS),1) = 192.8;
    NOA_E1_ID(1:height(E1_KINEMATICS),1) = {'9ON6'};
    NOA_E1_ALTID(1:height(E1_KINEMATICS),1) = 5;
    NOA_E1_VENT_DATA  = E1_VENT_DATA;
    NOA_E1_CC_DATA    = E1_CC_ENERGY;

    NOA_E3_TIME = E3_KINEMATICS.Time;
    NOA_E3_WSST_HR = E3_KINEMATICS{:,22:24};
    NOA_E3_EPOCH(1:height(E3_KINEMATICS),1) = 3;
    NOA_E3_MASS(1:height(E3_KINEMATICS),1) = 192.8;
    NOA_E3_ID(1:height(E3_KINEMATICS),1) = {'9ON6'};
    NOA_E3_ALTID(1:height(E3_KINEMATICS),1) = 5;
    NOA_E3_VENT_DATA  = E3_VENT_DATA;
    NOA_E3_CC_DATA    = E3_CC_ENERGY;
    
    NOA_E5_TIME = E5_KINEMATICS.Time;
    NOA_E5_WSST_HR = E5_KINEMATICS{:,22:24};
    NOA_E5_EPOCH(1:height(E5_KINEMATICS),1) = 5;
    NOA_E5_MASS(1:height(E5_KINEMATICS),1) = 192.8;
    NOA_E5_ID(1:height(E5_KINEMATICS),1) = {'9ON6'};
    NOA_E5_ALTID(1:height(E5_KINEMATICS),1) = 5;
    NOA_E5_VENT_DATA  = E5_VENT_DATA;
    NOA_E5_CC_DATA    = E5_CC_ENERGY;   

    % LONO
    
    load('/Users/dave/Desktop/DTAG/raw/tt21_142draw.mat');
    
    LONO_E1_TIME = E1_KINEMATICS.Time;
    LONO_E1_WSST_HR = E1_KINEMATICS{:,22:24};
    LONO_E1_EPOCH(1:height(E1_KINEMATICS),1) = 1;
    LONO_E1_MASS(1:height(E1_KINEMATICS),1) = 251.7;
    LONO_E1_ID(1:height(E1_KINEMATICS),1) = {'9FL3'};
    LONO_E1_ALTID(1:height(E1_KINEMATICS),1) = 6;
    LONO_E1_VENT_DATA  = E1_VENT_DATA;
    LONO_E1_CC_DATA    = E1_CC_ENERGY;

    LONO_E3_TIME = E3_KINEMATICS.Time;
    LONO_E3_WSST_HR = E3_KINEMATICS{:,22:24};
    LONO_E3_EPOCH(1:height(E3_KINEMATICS),1) = 3;
    LONO_E3_MASS(1:height(E3_KINEMATICS),1) = 251.7;
    LONO_E3_ID(1:height(E3_KINEMATICS),1) = {'9FL3'};
    LONO_E3_ALTID(1:height(E3_KINEMATICS),1) = 6;
    LONO_E3_VENT_DATA  = E3_VENT_DATA;
    LONO_E3_CC_DATA    = E3_CC_ENERGY;
    
    % Now aggregate everything...
    
    Time = [ HOKU_E1_TIME; HOKU_E3_TIME; ...
        LIHO_E1_TIME; LIHO_E3_TIME; LIHO_E5_TIME; ...
        KOLOHE_E1_TIME; KOLOHE_E3_TIME; KOLOHE_E5_TIME; ...
        HUA_E1_TIME; HUA_E3_TIME; HUA_E5_TIME; ...
        NOA_E1_TIME; NOA_E3_TIME; NOA_E5_TIME; ...
        LONO_E1_TIME; LONO_E3_TIME ];

    wsstHR = [ HOKU_E1_WSST_HR; HOKU_E3_WSST_HR; ...
        LIHO_E1_WSST_HR; LIHO_E3_WSST_HR; LIHO_E5_WSST_HR; ...
        KOLOHE_E1_WSST_HR; KOLOHE_E3_WSST_HR; KOLOHE_E5_WSST_HR; ...
        HUA_E1_WSST_HR; HUA_E3_WSST_HR; HUA_E5_WSST_HR; ...
        NOA_E1_WSST_HR; NOA_E3_WSST_HR; NOA_E5_WSST_HR; ...
        LONO_E1_WSST_HR; LONO_E3_WSST_HR ];
    
    epoch = [ HOKU_E1_EPOCH; HOKU_E3_EPOCH; ...
        LIHO_E1_EPOCH; LIHO_E3_EPOCH; LIHO_E5_EPOCH; ...
        KOLOHE_E1_EPOCH; KOLOHE_E3_EPOCH; KOLOHE_E5_EPOCH; ...
        HUA_E1_EPOCH; HUA_E3_EPOCH; HUA_E5_EPOCH; ...
        NOA_E1_EPOCH; NOA_E3_EPOCH; NOA_E5_EPOCH; ...
        LONO_E1_EPOCH; LONO_E3_EPOCH ];
    
    mass = [ HOKU_E1_MASS; HOKU_E3_MASS; ...
        LIHO_E1_MASS; LIHO_E3_MASS; LIHO_E5_MASS; ...
        KOLOHE_E1_MASS; KOLOHE_E3_MASS; KOLOHE_E5_MASS; ...
        HUA_E1_MASS; HUA_E3_MASS; HUA_E5_MASS; ...
        NOA_E1_MASS; NOA_E3_MASS; NOA_E5_MASS; ...
        LONO_E1_MASS; LONO_E3_MASS ];

    id = [ HOKU_E1_ID; HOKU_E3_ID; ...
        LIHO_E1_ID; LIHO_E3_ID; LIHO_E5_ID; ...
        KOLOHE_E1_ID; KOLOHE_E3_ID; KOLOHE_E5_ID; ...
        HUA_E1_ID; HUA_E3_ID; HUA_E5_ID; ...
        NOA_E1_ID; NOA_E3_ID; NOA_E5_ID; ...
        LONO_E1_ID; LONO_E3_ID ];
    
    altId = [ HOKU_E1_ALTID; HOKU_E3_ALTID; ...
        LIHO_E1_ALTID; LIHO_E3_ALTID; LIHO_E5_ALTID; ...
        KOLOHE_E1_ALTID; KOLOHE_E3_ALTID; KOLOHE_E5_ALTID; ...
        HUA_E1_ALTID; HUA_E3_ALTID; HUA_E5_ALTID; ...
        NOA_E1_ALTID; NOA_E3_ALTID; NOA_E5_ALTID; ...
        LONO_E1_ALTID; LONO_E3_ALTID ];    
    
    WSST_HR_DATA = timetable(Time, wsstHR(:,1), wsstHR(:,2), wsstHR(:,3), ...
        epoch, mass, id, altId, 'VariableNames', ...
        {'linWsstHr','rotWsstHr','ensWsstHr','epoch','mass','id','altId' } );
    
    WSST_HR_DATA.Time.Second = round(WSST_HR_DATA.Time.Second, 2);
    
    VENT_DATA = [ HOKU_E1_VENT_DATA; HOKU_E3_VENT_DATA; ...
        LIHO_E1_VENT_DATA; LIHO_E3_VENT_DATA; LIHO_E5_VENT_DATA; ...
        KOLOHE_E1_VENT_DATA; KOLOHE_E3_VENT_DATA; KOLOHE_E5_VENT_DATA; ...
        HUA_E1_VENT_DATA; HUA_E3_VENT_DATA; HUA_E5_VENT_DATA; ...
        NOA_E1_VENT_DATA; NOA_E3_VENT_DATA; NOA_E5_VENT_DATA; ...
        LONO_E1_VENT_DATA; LONO_E3_VENT_DATA];
        
    VENT_DATA.Time.Second = round(VENT_DATA.Time.Second, 2);
    
    CC_DATA = [ HOKU_E1_CC_DATA; HOKU_E3_CC_DATA; ...
        LIHO_E1_CC_DATA; LIHO_E3_CC_DATA; LIHO_E5_CC_DATA; ...
        KOLOHE_E1_CC_DATA; KOLOHE_E3_CC_DATA; KOLOHE_E5_CC_DATA; ...
        HUA_E1_CC_DATA; HUA_E3_CC_DATA; HUA_E5_CC_DATA; ...
        NOA_E1_CC_DATA; NOA_E3_CC_DATA; NOA_E5_CC_DATA; ...
        LONO_E1_CC_DATA; LONO_E3_CC_DATA ];
    
    CC_DATA.Time.Second = round(CC_DATA.Time.Second, 2);

    % make an aggregated version of CC_DATA + HR_DATA
    
    fprintf('Making CC_HR_DATA...\n');
    
    tempTT = synchronize(CC_DATA, HR_DATA);

    overlapRows = ( ~isnan(tempTT.iHR));

    index = 0;

    for zz = 1:length(overlapRows)

        if ( overlapRows(zz) == 1 )

            index = index + 1;
            getRow(index) = zz;

        end

    end

    CC_HR_DATA = tempTT(getRow,:); 
    
    CC_HR_DATA = removevars(CC_HR_DATA, { 'ID', 'mass_CC_DATA', ...
        'epoch_CC_DATA'} );
    
    CC_HR_DATA = splitvars(CC_HR_DATA, ...
        {'iLinKE_xyz','iRotKE_xyz', 'iTotalKE_xyz'}, 'NewVariableNames', ...
        { ...
        {'iLinKE_x','iLinKE_y','iLinKE_z'}, ...
        {'iRotKE_x','iRotKE_y','iRotKE_z'}, ...
        {'iTotalKE_x','iTotalKE_y','iTotalKE_z'} ...
        } );
    
    CC_HR_DATA = renamevars(CC_HR_DATA, ...
        {'iHR','epoch_HR_DATA','mass_HR_DATA' }, ...
        {'iHr','epoch','oldMass' } );
    
    
    %
    % ---=== HOKU ===---
    %


    HOKU        = struct;
    HOKU.mass   = 188;
    HOKU.length = 2.54;         % 254 cm
    HOKU.girth  = 1.3081;       % 130.81 cm

    %
    % ---=== LIHO ===---
    %

    LIHO        = struct;
    LIHO.mass   = 165;
    LIHO.length = 2.4003;         % 254 cm
    LIHO.girth  = 1.2426;       % 130.81 cm

    %
    % ---=== KOLOHE ===---
    %

    KOLOHE        = struct;
    KOLOHE.mass   = 211;
    KOLOHE.length = 2.5908;         % 254 cm
    KOLOHE.girth  = 1.4478;       % 130.81 cm

    %
    % ---=== HUA ===---
    %

    HUA        = struct;
    HUA.mass   = 154;
    HUA.length = 2.4130;         % 254 cm
    HUA.girth  = 1.2573;       % 130.81 cm

    %
    % ---=== NOA ===---
    %

    NOA        = struct;
    NOA.mass   = 187;
    NOA.length = 2.4765;         % 254 cm
    NOA.girth  = 1.3716;       % 130.81 cm

    %
    % ---=== LONO ===---
    %

    LONO        = struct;
    LONO.mass   = 228;
    LONO.length = 2.7432;         % 254 cm
    LONO.girth  = 1.4542;       % 130.81 cm    

    % build a reduced and updated 'mass' variable
    
    clear mass;
    
    mass(2517,1) = zeros;
    mass(1:185,1) = HOKU.mass;
    mass(186:720,1) = LIHO.mass;
    mass(721:1303,1) = KOLOHE.mass;
    mass(1304:1668,1) = HUA.mass;
    mass(1669:2200,1) = NOA.mass;
    mass(2201:2517,1) = LONO.mass;    
    
    len(2517,1) = zeros;
    len(1:185,1) = HOKU.length;
    len(186:720,1) = LIHO.length;
    len(721:1303,1) = KOLOHE.length;
    len(1304:1668,1) = HUA.length;
    len(1669:2200,1) = NOA.length;
    len(2201:2517,1) = LONO.length;        
    
    girth(2517,1) = zeros;
    girth(1:185,1) = HOKU.girth;
    girth(186:720,1) = LIHO.girth;
    girth(721:1303,1) = KOLOHE.girth;
    girth(1304:1668,1) = HUA.girth;
    girth(1669:2200,1) = NOA.girth;
    girth(2201:2517,1) = LONO.girth;    
    
    CC_HR_DATA = addvars(CC_HR_DATA, mass, len, girth, ...
        'NewVariableNames',{'mass','length','girth'}, ...
        'After','epoch' );    
    
    clear tempTT:
    clear overlapRows;
    clear getRow;
    
    % make an aggregated version of VENT_DATA + HR_DATA
    
    fprintf('Making VENT_HR_DATA...\n');
    
    tempTT = synchronize(VENT_DATA, HR_DATA);
    
    overlapRows = ( ~isnan(tempTT.ibi));

    index = 0;

    for zz = 1:length(overlapRows)

        if ( overlapRows(zz) == 1 )

            index = index + 1;
            getRow(index) = zz;

        end

    end

    VENT_HR_DATA = tempTT(getRow,:); 
    
    VENT_HR_DATA = removevars(VENT_HR_DATA, { 'id_VENT_DATA', ...
        'mass_VENT_DATA', 'epoch_VENT_DATA' } );
    
    VENT_HR_DATA = renamevars(VENT_HR_DATA, ...
        {'epoch_HR_DATA','mass_HR_DATA','id_HR_DATA' }, ...
        {'epoch', 'mass', 'id' } );
    
        
else
    
	fprintf('File already appears to exist... Loading it!\n');
    load(fullDataFileName);

end


%% Ask about saving or just save the stuff...


fprintf('Saving VENT_DATA, CC_DATA, HR_DATA, CC_HR_DATA & VENT_HR_DATA to the following file:\n');
fprintf('\t%s\n', fullDataFileName);

save(fullDataFileName, ...
    'VENT_DATA','CC_DATA', 'HR_DATA', 'CC_HR_DATA', 'VENT_HR_DATA', ...
    'HOKU_E1_CC_DATA', 'HOKU_E1_VENT_DATA', 'HOKU_E3_CC_DATA', ...
    'HOKU_E3_VENT_DATA', 'LIHO_E1_CC_DATA', 'LIHO_E1_VENT_DATA', ...
    'LIHO_E3_CC_DATA', 'LIHO_E3_VENT_DATA', 'LIHO_E5_CC_DATA', ...
    'LIHO_E5_VENT_DATA', 'KOLOHE_E1_CC_DATA', 'KOLOHE_E1_VENT_DATA', ...
    'KOLOHE_E3_CC_DATA', 'KOLOHE_E3_VENT_DATA', 'KOLOHE_E5_CC_DATA', ...
    'KOLOHE_E5_VENT_DATA', 'HUA_E1_CC_DATA', 'HUA_E1_VENT_DATA', ...
    'HUA_E3_CC_DATA', 'HUA_E3_VENT_DATA', 'HUA_E5_CC_DATA', ...
    'HUA_E5_VENT_DATA', 'NOA_E1_CC_DATA', 'NOA_E1_VENT_DATA', ...
    'NOA_E3_CC_DATA', 'NOA_E3_VENT_DATA', 'NOA_E5_CC_DATA', ...
    'NOA_E5_VENT_DATA', 'LONO_E1_CC_DATA', 'LONO_E1_VENT_DATA', ...
    'LONO_E3_CC_DATA', 'LONO_E3_VENT_DATA', 'WSST_HR_DATA'  );

fprintf('Writing HR_DATA to wsstHrData.csv\n');
writetimetable(HR_DATA, ...
    '/Users/dave/Documents/MATLAB/wsstHrData.csv');

fprintf('Writing VENT_DATA to ventData.csv\n');
writetimetable(VENT_DATA, ...
    '/Users/dave/Documents/MATLAB/ventData.csv');

fprintf('Writing CC_DATA to ccData.csv\n');
writetimetable(CC_DATA, ...
    '/Users/dave/Documents/MATLAB/ccData.csv');

fprintf('Writing CC_HR_DATA to ccHrData.csv\n');
writetimetable(CC_HR_DATA, ...
    '/Users/dave/Documents/MATLAB/ccHrData.csv');

fprintf('Writing VENT_HR_DATA to ventHrData.csv\n');
writetimetable(VENT_HR_DATA, ...
    '/Users/dave/Documents/MATLAB/ventHrData.csv');

stats = input('Press [enter] for summary stats...','s');

clc;
summary(WSST_HR_DATA);
summary(HR_DATA);
summary(CC_DATA);
summary(VENT_DATA);
summary(CC_HR_DATA);
summary(VENT_HR_DATA);


%% 
 
HOKU = '63HF';
LIHO = '01L5';
KOLOHE = '6JK5';
HUA = '83H1';
NOA = '9ON6';
LONO = '9FL3';

HokuStahlHr = ...
    unique(241 * HR_DATA.mass(string(HR_DATA.id) == HOKU) .^ -0.25);
LihoStahlHr = ...
    unique(241 * HR_DATA.mass(string(HR_DATA.id) == LIHO) .^ -0.25);
KoloheStahlHr = ...
    unique(241 * HR_DATA.mass(string(HR_DATA.id) == KOLOHE) .^ -0.25);
HuaStahlHr = ...
    unique(241 * HR_DATA.mass(string(HR_DATA.id) == HUA) .^ -0.25);
NoaStahlHr = ...
    unique(241 * HR_DATA.mass(string(HR_DATA.id) == NOA) .^ -0.25);
LonoStahlHr = ...
    unique(241 * HR_DATA.mass(string(HR_DATA.id) == LONO) .^ -0.25);

meanStahlHr = mean( [ HokuStahlHr LihoStahlHr KoloheStahlHr ...
    HuaStahlHr NoaStahlHr LonoStahlHr] );


HOKU_E1_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 1 & contains(HR_DATA.id, {HOKU}) );
HOKU_E3_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 3 & contains(HR_DATA.id, {HOKU}) );

LIHO_E1_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 1 & contains(HR_DATA.id, {LIHO}) );
LIHO_E3_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 3 & contains(HR_DATA.id, {LIHO}) );
LIHO_E5_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 5 & contains(HR_DATA.id, {LIHO}) );

KOLOHE_E1_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 1 & contains(HR_DATA.id, {KOLOHE}) );
KOLOHE_E3_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 3 & contains(HR_DATA.id, {KOLOHE}) );
KOLOHE_E5_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 5 & contains(HR_DATA.id, {KOLOHE}) );   

HUA_E1_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 1 & contains(HR_DATA.id, {HUA}) );
HUA_E3_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 3 & contains(HR_DATA.id, {HUA}) );
HUA_E5_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 5 & contains(HR_DATA.id, {HUA}) );       

NOA_E1_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 1 & contains(HR_DATA.id, {NOA}) );
NOA_E3_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 3 & contains(HR_DATA.id, {NOA}) );
NOA_E5_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 5 & contains(HR_DATA.id, {NOA}) );       

LONO_E1_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 1 & contains(HR_DATA.id, {LONO}) );
LONO_E3_HR = ...
    HR_DATA.ensWsstHr(HR_DATA.epoch == 3 & contains(HR_DATA.id, {LONO}) );

%% Compile stats

% ========= Hoku stats ==========

mean_Hoku_E1_wsstHR = mean(HOKU_E1_HR);
mean_Hoku_E3_wsstHR = mean(HOKU_E3_HR);
mean_Hoku_E1_iHR = mean(HOKU_E1_CC_DATA.iHR);
mean_Hoku_E3_iHR = mean(HOKU_E3_CC_DATA.iHR);

sem_Hoku_E1_wsstHR  = sem(HOKU_E1_HR);
sem_Hoku_E3_wsstHR  = sem(HOKU_E3_HR);
sem_Hoku_E1_iHR  = sem(HOKU_E1_CC_DATA.iHR);
sem_Hoku_E3_iHR  = sem(HOKU_E3_CC_DATA.iHR);

min_Hoku_E1_wsstHR = min(HOKU_E1_HR);
max_Hoku_E1_wsstHR = max(HOKU_E1_HR);
min_Hoku_E3_wsstHR = min(HOKU_E3_HR);
max_Hoku_E3_wsstHR = max(HOKU_E3_HR);

min_Hoku_E1_iHR = min(HOKU_E1_CC_DATA.iHR);
max_Hoku_E1_iHR = max(HOKU_E1_CC_DATA.iHR);
min_Hoku_E3_iHR = min(HOKU_E3_CC_DATA.iHR);
max_Hoku_E3_iHR = max(HOKU_E3_CC_DATA.iHR);

meanIBI_Hoku_E1 = mean(HOKU_E1_VENT_DATA.ibi);
semIBI_Hoku_E1  = sem(HOKU_E1_VENT_DATA.ibi);

apneaDuration_Hoku_E3 = seconds( HOKU_E3_TIME(end) - HOKU_E3_TIME(1) );


% ========= Liho stats ==========

mean_Liho_E1_wsstHR = mean(LIHO_E1_HR);
mean_Liho_E3_wsstHR = mean(LIHO_E3_HR);
mean_Liho_E5_wsstHR = mean(LIHO_E5_HR);
mean_Liho_E1_iHR = mean(LIHO_E1_CC_DATA.iHR);
mean_Liho_E3_iHR = mean(LIHO_E3_CC_DATA.iHR);
mean_Liho_E5_iHR = mean(LIHO_E5_CC_DATA.iHR);

sem_Liho_E1_wsstHR  = sem(LIHO_E1_HR);
sem_Liho_E3_wsstHR  = sem(LIHO_E3_HR);
sem_Liho_E5_wsstHR  = sem(LIHO_E5_HR);
sem_Liho_E1_iHR  = sem(LIHO_E1_CC_DATA.iHR);
sem_Liho_E3_iHR  = sem(LIHO_E3_CC_DATA.iHR);
sem_Liho_E5_iHR  = sem(LIHO_E5_CC_DATA.iHR);

min_Liho_E1_wsstHR = min(LIHO_E1_HR);
max_Liho_E1_wsstHR = max(LIHO_E1_HR);
min_Liho_E3_wsstHR = min(LIHO_E3_HR);
max_Liho_E3_wsstHR = max(LIHO_E3_HR);
min_Liho_E5_wsstHR = min(LIHO_E5_HR);
max_Liho_E5_wsstHR = max(LIHO_E5_HR);

min_Liho_E1_iHR = min(LIHO_E1_CC_DATA.iHR);
max_Liho_E1_iHR = max(LIHO_E1_CC_DATA.iHR);
min_Liho_E3_iHR = min(LIHO_E3_CC_DATA.iHR);
max_Liho_E3_iHR = max(LIHO_E3_CC_DATA.iHR);
min_Liho_E5_iHR = min(LIHO_E5_CC_DATA.iHR);
max_Liho_E5_iHR = max(LIHO_E5_CC_DATA.iHR);

meanIBI_Liho_E1 = mean(LIHO_E1_VENT_DATA.ibi);
meanIBI_Liho_E5 = mean(LIHO_E5_VENT_DATA.ibi);
semIBI_Liho_E1 = sem(LIHO_E1_VENT_DATA.ibi);
semIBI_Liho_E5 = sem(LIHO_E5_VENT_DATA.ibi);

apneaDuration_Liho_E3 = seconds( LIHO_E3_TIME(end) - LIHO_E3_TIME(1) );


% ========= Kolohe stats ==========

mean_Kolohe_E1_HR = mean(KOLOHE_E1_HR);
mean_Kolohe_E3_HR = mean(KOLOHE_E3_HR);
mean_Kolohe_E5_HR = mean(KOLOHE_E5_HR);
mean_Kolohe_E1_iHR = mean(KOLOHE_E1_CC_DATA.iHR);
mean_Kolohe_E3_iHR = mean(KOLOHE_E3_CC_DATA.iHR);
mean_Kolohe_E5_iHR = mean(KOLOHE_E5_CC_DATA.iHR);

sem_Kolohe_E1_wsstHR  = sem(KOLOHE_E1_HR);
sem_Kolohe_E3_wsstHR  = sem(KOLOHE_E3_HR);
sem_Kolohe_E5_wsstHR  = sem(KOLOHE_E5_HR);
sem_Kolohe_E1_iHR  = sem(KOLOHE_E1_CC_DATA.iHR);
sem_Kolohe_E3_iHR  = sem(KOLOHE_E3_CC_DATA.iHR);
sem_Kolohe_E5_iHR  = sem(KOLOHE_E5_CC_DATA.iHR);

min_Kolohe_E1_wsstHR = min(KOLOHE_E1_HR);
max_Kolohe_E1_wsstHR = max(KOLOHE_E1_HR);
min_Kolohe_E3_wsstHR = min(KOLOHE_E3_HR);
max_Kolohe_E3_wsstHR = max(KOLOHE_E3_HR);
min_Kolohe_E5_wsstHR = min(KOLOHE_E5_HR);
max_Kolohe_E5_wsstHR = max(KOLOHE_E5_HR);

min_Kolohe_E1_iHR = min(KOLOHE_E1_CC_DATA.iHR);
max_Kolohe_E1_iHR = max(KOLOHE_E1_CC_DATA.iHR);
min_Kolohe_E3_iHR = min(KOLOHE_E3_CC_DATA.iHR);
max_Kolohe_E3_iHR = max(KOLOHE_E3_CC_DATA.iHR);
min_Kolohe_E5_iHR = min(KOLOHE_E5_CC_DATA.iHR);
max_Kolohe_E5_iHR = max(KOLOHE_E5_CC_DATA.iHR);

meanIBI_Kolohe_E1 = mean(KOLOHE_E1_VENT_DATA.ibi);
meanIBI_Kolohe_E5 = mean(KOLOHE_E5_VENT_DATA.ibi);
semIBI_Kolohe_E1 = sem(KOLOHE_E1_VENT_DATA.ibi);
semIBI_Kolohe_E5 = sem(KOLOHE_E5_VENT_DATA.ibi);

apneaDuration_Kolohe_E3 = seconds( KOLOHE_E3_TIME(end) - KOLOHE_E3_TIME(1) );


% ========= Hua stats ==========

mean_Hua_E1_wsstHR = mean(HUA_E1_HR);
mean_Hua_E3_wsstHR = mean(HUA_E3_HR);
mean_Hua_E5_wsstHR = mean(HUA_E5_HR);
mean_Hua_E1_iHR = mean(HUA_E1_CC_DATA.iHR);
mean_Hua_E3_iHR = mean(HUA_E3_CC_DATA.iHR);
mean_Hua_E5_iHR = mean(HUA_E5_CC_DATA.iHR);

sem_Hua_E1_wsstHR  = sem(HUA_E1_HR);
sem_Hua_E3_wsstHR  = sem(HUA_E3_HR);
sem_Hua_E5_wsstHR  = sem(HUA_E5_HR);
sem_Hua_E1_iHR  = sem(HUA_E1_CC_DATA.iHR);
sem_Hua_E3_iHR  = sem(HUA_E3_CC_DATA.iHR);
sem_Hua_E5_iHR  = sem(HUA_E5_CC_DATA.iHR);

min_Hua_E1_wsstHR = min(HUA_E1_HR);
max_Hua_E1_wsstHR = max(HUA_E1_HR);
min_Hua_E3_wsstHR = min(HUA_E3_HR);
max_Hua_E3_wsstHR = max(HUA_E3_HR);
min_Hua_E5_wsstHR = min(HUA_E5_HR);
max_Hua_E5_wsstHR = max(HUA_E5_HR);

min_Hua_E1_iHR = min(HUA_E1_CC_DATA.iHR);
max_Hua_E1_iHR = max(HUA_E1_CC_DATA.iHR);
min_Hua_E3_iHR = min(HUA_E3_CC_DATA.iHR);
max_Hua_E3_iHR = max(HUA_E3_CC_DATA.iHR);
min_Hua_E5_iHR = min(HUA_E5_CC_DATA.iHR);
max_Hua_E5_iHR = max(HUA_E5_CC_DATA.iHR);

meanIBI_Hua_E1 = mean(HUA_E1_VENT_DATA.ibi);
meanIBI_Hua_E5 = mean(HUA_E5_VENT_DATA.ibi);
semIBI_Hua_E1 = sem(HUA_E1_VENT_DATA.ibi);
semIBI_Hua_E5 = sem(HUA_E5_VENT_DATA.ibi);

apneaDuration_Hua_E3 = seconds( HUA_E3_TIME(end) - HUA_E3_TIME(1) );


% ========= Noa stats ==========

mean_Noa_E1_wsstHR = mean(NOA_E1_HR);
mean_Noa_E3_wsstHR = mean(NOA_E3_HR);
mean_Noa_E5_wsstHR = mean(NOA_E5_HR);
mean_Noa_E1_iHR = mean(NOA_E1_CC_DATA.iHR);
mean_Noa_E3_iHR = mean(NOA_E3_CC_DATA.iHR);
mean_Noa_E5_iHR = mean(NOA_E5_CC_DATA.iHR);

sem_Noa_E1_wsstHR  = sem(NOA_E1_HR);
sem_Noa_E3_wsstHR  = sem(NOA_E3_HR);
sem_Noa_E5_wsstHR  = sem(NOA_E5_HR);
sem_Noa_E1_iHR  = sem(NOA_E1_CC_DATA.iHR);
sem_Noa_E3_iHR  = sem(NOA_E3_CC_DATA.iHR);
sem_Noa_E5_iHR  = sem(NOA_E5_CC_DATA.iHR);

min_Noa_E1_wsstHR = min(NOA_E1_HR);
max_Noa_E1_wsstHR = max(NOA_E1_HR);
min_Noa_E3_wsstHR = min(NOA_E3_HR);
max_Noa_E3_wsstHR = max(NOA_E3_HR);
min_Noa_E5_wsstHR = min(NOA_E5_HR);
max_Noa_E5_wsstHR = max(NOA_E5_HR);

min_Noa_E1_iHR = min(NOA_E1_CC_DATA.iHR);
max_Noa_E1_iHR = max(NOA_E1_CC_DATA.iHR);
min_Noa_E3_iHR = min(NOA_E3_CC_DATA.iHR);
max_Noa_E3_iHR = max(NOA_E3_CC_DATA.iHR);
min_Noa_E5_iHR = min(NOA_E5_CC_DATA.iHR);
max_Noa_E5_iHR = max(NOA_E5_CC_DATA.iHR);

meanIBI_Noa_E1 = mean(NOA_E1_VENT_DATA.ibi);
meanIBI_Noa_E5 = mean(NOA_E5_VENT_DATA.ibi);
semIBI_Noa_E1 = sem(NOA_E1_VENT_DATA.ibi);
semIBI_Noa_E5 = sem(NOA_E5_VENT_DATA.ibi);

apneaDuration_Noa_E3 = seconds( NOA_E3_TIME(end) - NOA_E3_TIME(1) );


% ========= Lono stats ==========

mean_Lono_E1_wsstHR = mean(LONO_E1_HR);
mean_Lono_E3_wsstHR = mean(LONO_E3_HR);
mean_Lono_E1_iHR = mean(LONO_E1_CC_DATA.iHR);
mean_Lono_E3_iHR = mean(LONO_E3_CC_DATA.iHR);

sem_Lono_E1_wsstHR  = sem(LONO_E1_HR);
sem_Lono_E3_wsstHR  = sem(LONO_E3_HR);
sem_Lono_E1_iHR  = sem(LONO_E1_CC_DATA.iHR);
sem_Lono_E3_iHR  = sem(LONO_E3_CC_DATA.iHR);

min_Lono_E1_wsstHR = min(LONO_E1_HR);
max_Lono_E1_wsstHR = max(LONO_E1_HR);
min_Lono_E3_wsstHR = min(LONO_E3_HR);
max_Lono_E3_wsstHR = max(LONO_E3_HR);
min_Lono_E1_iHR = min(LONO_E1_CC_DATA.iHR);
max_Lono_E1_iHR = max(LONO_E1_CC_DATA.iHR);
min_Lono_E3_iHR = min(LONO_E3_CC_DATA.iHR);
max_Lono_E3_iHR = max(LONO_E3_CC_DATA.iHR);

meanIBI_Lono_E1 = mean(LONO_E1_VENT_DATA.ibi);
semIBI_Lono_E1 = sem(LONO_E1_VENT_DATA.ibi);

apneaDuration_Lono_E3 = seconds( LONO_E3_TIME(end) - LONO_E3_TIME(1) );

%% Make a nice table of all these values...
%                                        HOKU  LIHO  KOLOHE  HUA  NOA  LONO
% mean_E1_wsstHr    w/ sem_E1_wsstHr
% mean_E3_wsstHr    w/ sem_E3_wsstHr
% mean_E5_wsstHr    w/ sem_E5_wsstHr
% mean_E1_iHr       w/ sem_E1_iHr
% mean_E3_iHr       w/ sem_E3_iHr
% mean_E5_iHr       w/ sem_E5_iHr

% min_E1_wsstHr
% min_E3_wsstHr
% min_E5_wsstHr
% min_E1_iHr
% min_E3_iHr
% min_E5_iHr

% max_E1_wsstHr
% max_E3_wsstHr
% max_E5_wsstHr
% max_E1_iHr
% max_E3_iHr
% max_E5_iHr

% mean_E1_ibi       w/ sem_E1_ibi  
% mean_E3_ibi       (- no sem  -)
% mean_E5_ibi       w/ sem_E5_ibi



%% IBI, mean ifH, min ifH, and max ifH by animal ID

figure('Color', [1 1 1]);
s1 = subplot(411);
boxplot( VENT_DATA.ibi, VENT_DATA.id, ...
    'Labels', {'63HF','01L5','6JK5','83H1','9ON6','9FL3'} );
xlabel('Animal ID');
ylabel('Time, seconds');
title('Inter-breath interval');
grid;
s2 = subplot(412);
boxplot( VENT_DATA.mean_ifH, VENT_DATA.id);
xlabel('Animal ID');
ylabel('beats · min^-^1');
title('Mean instantaneous heart rate');
grid;
s3 = subplot(413);
boxplot( VENT_DATA.min_ifH, VENT_DATA.id);
xlabel('Animal ID');
ylabel('beats · min^-^1');
title('Minimum instantaneous heart rate');
grid;
s4 = subplot(414);
boxplot( VENT_DATA.max_ifH, VENT_DATA.id);
xlabel('Animal ID');
ylabel('beats · min^-^1');
title('Maximum instantaneous heart rate');
grid;

%% IBI, mean ifH, min ifH, and max ifH by epoch

figure('Color', [1 1 1]);
s1 = subplot(411);
boxplot( VENT_DATA.ibi, VENT_DATA.epoch, ...
    'Labels', {'Baseline','Apnea','Recovery'} );
xlabel('Epoch');
ylabel('Time, seconds');
title('Inter-breath interval');
grid;
s2 = subplot(412);
boxplot( VENT_DATA.mean_ifH, VENT_DATA.epoch, ...
    'Labels', {'Baseline','Apnea','Recovery'} );
xlabel('Epoch');
ylabel('beats · min^-^1');
title('Mean instantaneous heart rate');
grid;
s3 = subplot(413);
boxplot( VENT_DATA.min_ifH, VENT_DATA.epoch, ...
    'Labels', {'Baseline','Apnea','Recovery'} );
xlabel('Epoch');
ylabel('beats · min^-^1');
title('Minimum instantaneous heart rate');
grid;
s4 = subplot(414);
boxplot( VENT_DATA.max_ifH, VENT_DATA.epoch, ...
    'Labels', {'Baseline','Apnea','Recovery'} );
xlabel('Epoch');
ylabel('beats · min^-^1');
title('Maximum instantaneous heart rate');
grid;

%%

figure('Color',[1 1 1]);
subplot(311);
boxplot(CC_HR_DATA.ensWsstHr(CC_HR_DATA.epoch == 1), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 1), 'Notch','on');
xlabel('Animal ID');
ylabel('beats · min^-^1');
title('Baseline Epoch: instantaneous heart rate');
grid;
subplot(312);
boxplot(CC_HR_DATA.ensWsstHr(CC_HR_DATA.epoch == 3), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 3), 'Notch','on');
xlabel('Animal ID');
ylabel('beats · min^-^1');
title('Apnea Epoch: instantaneous heart rate');
grid;
subplot(313);
boxplot(CC_HR_DATA.ensWsstHr(CC_HR_DATA.epoch == 5), ...
    CC_HR_DATA.id(CC_HR_DATA.epoch == 5), 'Notch','on');
xlabel('Animal ID');
ylabel('beats · min^-^1');
title('Recovery Epoch: instantaneous heart rate');
grid;

%% Run some ANOVAs

% a two-way ANOVA with unbalanced design to compare heart rates during
% baseline, apnea, and recovery epochs

y1 = CC_HR_DATA.iHr;
y2 = CC_HR_DATA.ensWsstHr;
g1 = CC_HR_DATA.epoch;
g2 = CC_HR_DATA.altId;

[pAOV1, tblAOV1, statsAOV1, termsAOV1] = ...
    anovan( y2, {g1,g2}, ...
    'model','interaction', ...
    'varnames', {'Epoch','AnimalID'} );

figure('Color', white);
[cAOV1, mAOV1, hAOV1, nmsAOV1] = ...
    multcompare(statsAOV1, 'CType', 'hsd','Dimension',[1 2]);

%

[pAOV2, tblAOV2, statsAOV2, termsAOV2] = ...
    anovan(CC_HR_DATA.linWsstHr, CC_HR_DATA.rotWsstHr);

figure('Color', white);
[cAOV2, mAOV2, hAOV2, nmsAOV2] = ...
    multcompare(statsAOV1, 'CType', 'hsd','Dimension',[1 2]);

%

[pAOV3, tblAOV3, statsAOV3, termsAOV3] = ...
    anovan(CC_HR_DATA.ensWsstHr, CC_HR_DATA.iHr);
figure('Color', white);
[cAOV3, mAOV3, hAOV3, nmsAOV3] = ...
    multcompare(statsAOV1, 'CType', 'hsd','Dimension',[1 2]);


%% 

