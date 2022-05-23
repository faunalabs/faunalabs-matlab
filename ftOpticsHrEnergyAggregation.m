%% ftOpticsHrEnergyAggregation.m
%
%   Written by Dave Haas on 17 October 2021

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
    
    getOk = lower(input('File exists, reprocess? (y/n) [y]: ','s'));
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

    HOKU_E1_DECIMATED  = E1_DECIMATED;
    HOKU_E3_DECIMATED  = E3_DECIMATED;
    
    HOKU_E1_CC_OPTICS = E1_CC_OPTICS;
    HOKU_E3_CC_OPTICS = E3_CC_OPTICS;
    
    % no E5 - tag came off during fast swim transition
    
    % LIHO
    
    load('/Users/dave/Desktop/DTAG/raw/tt21_134araw.mat');

    LIHO_E1_DECIMATED  = E1_DECIMATED;
    LIHO_E3_DECIMATED  = E3_DECIMATED;
    LIHO_E5_DECIMATED  = E5_DECIMATED;
    
    LIHO_E1_CC_OPTICS = E1_CC_OPTICS;
    LIHO_E3_CC_OPTICS = E3_CC_OPTICS;
    LIHO_E5_CC_OPTICS = E5_CC_OPTICS;    
    
    % KOLOHE
    
    load('/Users/dave/Desktop/DTAG/raw/tt21_134braw.mat');
    
    KOLOHE_E1_DECIMATED  = E1_DECIMATED;
    KOLOHE_E3_DECIMATED  = E3_DECIMATED;
    KOLOHE_E5_DECIMATED  = E5_DECIMATED;
    
    KOLOHE_E1_CC_OPTICS = E1_CC_OPTICS;
    KOLOHE_E3_CC_OPTICS = E3_CC_OPTICS;
    KOLOHE_E5_CC_OPTICS = E5_CC_OPTICS;    
    
    % HUA
    
    load('/Users/dave/Desktop/DTAG/raw/tt21_141eraw.mat');

    HUA_E1_DECIMATED  = E1_DECIMATED;
    HUA_E3_DECIMATED  = E3_DECIMATED;
    HUA_E5_DECIMATED  = E5_DECIMATED;
    
    HUA_E1_CC_OPTICS = E1_CC_OPTICS;
    HUA_E3_CC_OPTICS = E3_CC_OPTICS;
    HUA_E5_CC_OPTICS = E5_CC_OPTICS;    

    % NOA
    
    load('/Users/dave/Desktop/DTAG/raw/tt21_142craw.mat');
    
    NOA_E1_DECIMATED  = E1_DECIMATED;
    NOA_E3_DECIMATED  = E3_DECIMATED;
    NOA_E5_DECIMATED  = E5_DECIMATED;
    
    NOA_E1_CC_OPTICS = E1_CC_OPTICS;
    NOA_E3_CC_OPTICS = E3_CC_OPTICS;
    NOA_E5_CC_OPTICS = E5_CC_OPTICS;  

    % LONO
    
    load('/Users/dave/Desktop/DTAG/raw/tt21_142draw.mat');
    
    LONO_E1_DECIMATED  = E1_DECIMATED;
    LONO_E3_DECIMATED  = E3_DECIMATED;
    LONO_E5_DECIMATED  = E5_DECIMATED;
    
    LONO_E1_CC_OPTICS = E1_CC_OPTICS;
    LONO_E3_CC_OPTICS = E3_CC_OPTICS;
    LONO_E5_CC_OPTICS = E5_CC_OPTICS;  
    
    % Now aggregate everything...
    
    % put together the entire decimated data time series
    
    fprintf('Making DEC_DATA (50 Hz decimated full sensor suite)...\n');
    
    DEC_DATA = [ HOKU_E1_DECIMATED; HOKU_E3_DECIMATED; ...
        LIHO_E1_DECIMATED; LIHO_E3_DECIMATED; LIHO_E5_DECIMATED; ...
        KOLOHE_E1_DECIMATED; KOLOHE_E3_DECIMATED; KOLOHE_E5_DECIMATED; ...
        HUA_E1_DECIMATED; HUA_E3_DECIMATED; HUA_E5_DECIMATED; ...
        NOA_E1_DECIMATED; NOA_E3_DECIMATED; NOA_E5_DECIMATED; ...
        LONO_E1_DECIMATED; LONO_E3_DECIMATED ];
        
    DEC_DATA.Time.Second = round(DEC_DATA.Time.Second, 2);
    
    % make an aggregated version of CC_DATA + HR_DATA
    
    fprintf('Making CC_OPTICS_DATA...\n');
    
    CC_OPTICS_DATA = [ HOKU_E1_CC_OPTICS; HOKU_E3_CC_OPTICS; ...
        LIHO_E1_CC_OPTICS; LIHO_E3_CC_OPTICS; LIHO_E5_CC_OPTICS; ...
        KOLOHE_E1_CC_OPTICS; KOLOHE_E3_CC_OPTICS; KOLOHE_E5_CC_OPTICS; ...
        HUA_E1_CC_OPTICS; HUA_E3_CC_OPTICS; HUA_E5_CC_OPTICS; ...
        NOA_E1_CC_OPTICS; NOA_E3_CC_OPTICS; NOA_E5_CC_OPTICS; ...
        LONO_E1_CC_OPTICS; LONO_E3_CC_OPTICS ];
    
    CC_OPTICS_DATA.Time.Second = round(CC_OPTICS_DATA.Time.Second, 2);
    
        
else
    
	fprintf('File already appears to exist... Loading it!\n');
    load(fullDataFileName);

end


%% Because you forgot, build SNR for led2, led3, and led4 w/ led1 as noise

snr_dbFiltLed2 = CC_OPTICS_DATA.dbFiltLed2 - CC_OPTICS_DATA.dbFiltLed1;
snr_dbFiltLed3 = CC_OPTICS_DATA.dbFiltLed3 - CC_OPTICS_DATA.dbFiltLed1;
snr_dbFiltLed4 = CC_OPTICS_DATA.dbFiltLed4 - CC_OPTICS_DATA.dbFiltLed1;

CC_OPTICS_DATA = addvars(CC_OPTICS_DATA, snr_dbFiltLed2, snr_dbFiltLed3, ...
    snr_dbFiltLed4, 'After', 'dbFiltLed4', 'NewVariableNames', ...
    {'snr_dbFiltLed2', 'snr_dbFiltLed3', 'snr_dbFiltLed4'} );


%% Ask about saving or just save the stuff...

dontAsk = false;

if (dontAsk)
   askAns = lower(input('Save everything (y/n) [y]: ','s'));
   switch(askAns)
       case 'n'
           fprintf('Okay. Come back to this section to save.\n');
           okToWrite = false;
       otherwise
           fprintf('Okay, saving it all!\n');
           okToWrite = true;
   end
else
    okToWrite = true;
end


%%

if (okToWrite)
    
    fprintf('Appending DEC_DATA and CC_OPTICS to the following file:\n');
    fprintf('\t%s\n', fullDataFileName);

    save(fullDataFileName, ...
        'DEC_DATA','CC_OPTICS_DATA', ...
        'HOKU_E1_CC_OPTICS', 'HOKU_E3_CC_OPTICS', ...                               
        'HOKU_E1_DECIMATED', 'HOKU_E3_DECIMATED', ...
        'LIHO_E1_CC_OPTICS', 'LIHO_E3_CC_OPTICS', 'LIHO_E5_CC_OPTICS', ...          
        'LIHO_E1_DECIMATED', 'LIHO_E3_DECIMATED', 'LIHO_E5_DECIMATED', ...
        'KOLOHE_E1_CC_OPTICS', 'KOLOHE_E3_CC_OPTICS', 'KOLOHE_E5_CC_OPTICS', ...    
        'KOLOHE_E1_DECIMATED', 'KOLOHE_E3_DECIMATED', 'KOLOHE_E5_DECIMATED', ...
        'HUA_E1_CC_OPTICS', 'HUA_E3_CC_OPTICS', 'HUA_E5_CC_OPTICS', ...             
        'HUA_E1_DECIMATED', 'HUA_E3_DECIMATED', 'HUA_E5_DECIMATED', ...
        'NOA_E1_CC_OPTICS', 'NOA_E3_CC_OPTICS', 'NOA_E5_CC_OPTICS', ...             
        'NOA_E1_DECIMATED', 'NOA_E3_DECIMATED', 'NOA_E5_DECIMATED', ...    
        'LONO_E1_CC_OPTICS', 'LONO_E3_CC_OPTICS', ...                               
        'LONO_E1_DECIMATED', 'LONO_E3_DECIMATED', ...    
        '-append' );


fprintf('Writing CC_OPTICS_DATA to ccOpticsData.csv\n');
writetimetable(CC_OPTICS_DATA, ...
    '/Users/dave/Documents/MATLAB/ccOpticsData.csv');


stats = input('Press [enter] for summary stats...','s');

clc;
summary(CC_OPTICS_DATA);

end

