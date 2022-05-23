%% ftVentStats.m
%
%   Written by Dave Haas on 14 September

clc;
clear;


%% 

fullDataFileName = '/Users/dave/Documents/MATLAB/dqoVentTrialData.mat';

if (~exist(fullDataFileName, 'file'))

    % aggregate all the vent data into a single timetable
    
    load('/Users/dave/Desktop/DTAG/raw/tt21_128fraw.mat');
    HOKU_E1_VENT_DATA = E1_VENT_DATA;
    HOKU_E3_VENT_DATA = E3_VENT_DATA;

    load('/Users/dave/Desktop/DTAG/raw/tt21_134araw.mat');
    LIHO_E1_VENT_DATA = E1_VENT_DATA;
    LIHO_E3_VENT_DATA = E3_VENT_DATA;
    LIHO_E5_VENT_DATA = E5_VENT_DATA;

    load('/Users/dave/Desktop/DTAG/raw/tt21_134braw.mat');
    KOLOHE_E1_VENT_DATA = E1_VENT_DATA;
    KOLOHE_E3_VENT_DATA = E3_VENT_DATA;    
    KOLOHE_E5_VENT_DATA = E5_VENT_DATA;

    load('/Users/dave/Desktop/DTAG/raw/tt21_141eraw.mat');
    HUA_E1_VENT_DATA = E1_VENT_DATA;
    HUA_E3_VENT_DATA = E3_VENT_DATA;
    HUA_E5_VENT_DATA = E5_VENT_DATA;

    load('/Users/dave/Desktop/DTAG/raw/tt21_142craw.mat');
    NOA_E1_VENT_DATA = E1_VENT_DATA;
    NOA_E3_VENT_DATA = E3_VENT_DATA;
    NOA_E5_VENT_DATA = E5_VENT_DATA;

    load('/Users/dave/Desktop/DTAG/raw/tt21_142draw.mat');
    LONO_E1_VENT_DATA = E1_VENT_DATA;
    LONO_E3_VENT_DATA = E3_VENT_DATA;

    VENT_DATA = [ HOKU_E1_VENT_DATA; HOKU_E3_VENT_DATA; ...
        LIHO_E1_VENT_DATA; LIHO_E3_VENT_DATA; LIHO_E5_VENT_DATA; ...
        KOLOHE_E1_VENT_DATA; KOLOHE_E3_VENT_DATA; KOLOHE_E5_VENT_DATA; ...
        HUA_E1_VENT_DATA; HUA_E3_VENT_DATA; HUA_E5_VENT_DATA; ...
        NOA_E1_VENT_DATA; NOA_E3_VENT_DATA; NOA_E5_VENT_DATA; ...
        LONO_E1_VENT_DATA; LONO_E3_VENT_DATA];
    
	VENT_DATA.Time.Second = round(VENT_DATA.Time.Second, 2);
    
    save(fullDataFileName, 'VENT_DATA');
    
    writetimetable(VENT_DATA, '/Users/dave/Documents/MATLAB/ventData.csv');
    
else
    
	fprintf('File already appears to exist... Loading it!\n');
    load(fullDataFileName);

end


%% IBI, mean ifH, min ifH, and max ifH by animal ID

figure('Color', [1 1 1]);
s1 = subplot(411);
boxplot( VENT_DATA.ibi, VENT_DATA.id);
xlabel('IBI');
grid;
s2 = subplot(412);
boxplot( VENT_DATA.mean_ifH, VENT_DATA.id);
xlabel('Mean ifH');
grid;
s3 = subplot(413);
boxplot( VENT_DATA.min_ifH, VENT_DATA.id);
xlabel('Min ifH');
grid;
s4 = subplot(414);
boxplot( VENT_DATA.max_ifH, VENT_DATA.id);
xlabel('Max ifH');
grid;

%% IBI, mean ifH, min ifH, and max ifH by epoch

figure('Color', [1 1 1]);
s1 = subplot(411);
boxplot( VENT_DATA.ibi, VENT_DATA.epoch);
xlabel('IBI');
grid;
s2 = subplot(412);
boxplot( VENT_DATA.mean_ifH, VENT_DATA.epoch);
xlabel('Mean ifH');
grid;
s3 = subplot(413);
boxplot( VENT_DATA.min_ifH, VENT_DATA.epoch);
xlabel('Min ifH');
grid;
s4 = subplot(414);
boxplot( VENT_DATA.max_ifH, VENT_DATA.epoch);
xlabel('Max ifH');
grid;




