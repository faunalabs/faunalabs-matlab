
% E1

e1_meanOdba = mean(LONO.E1_KINEMATICS.hrTimeConsensusOdba) ;
e1_stdOdba  = std(LONO.E1_KINEMATICS.hrTimeConsensusOdba) ;

e1_meanOdav = mean(LONO.E1_KINEMATICS.hrTimeConsensusOdav) ;
e1_stdOdav  = std(LONO.E1_KINEMATICS.hrTimeConsensusOdav) ;

e1_meanHr = mean(LONO.E1_KINEMATICS.hrTimeConsensus) ;
e1_stdHr  = std(LONO.E1_KINEMATICS.hrTimeConsensus) ;

% E3

e3_meanOdba = mean(LONO.E3_KINEMATICS.hrTimeConsensusOdba) ;
e3_stdOdba  = std(LONO.E3_KINEMATICS.hrTimeConsensusOdba) ;

e3_meanOdav = mean(LONO.E3_KINEMATICS.hrTimeConsensusOdav) ;
e3_stdOdav  = std(LONO.E3_KINEMATICS.hrTimeConsensusOdav) ;


e3_meanHr = mean(LONO.E3_KINEMATICS.hrTimeConsensus) ;
e3_stdHr  = std(LONO.E3_KINEMATICS.hrTimeConsensus) ;

% E5
% 
% e5_meanOdba = mean(LONO.E5_KINEMATICS.hrTimeConsensusOdba) ;
% e5_stdOdba  = std(LONO.E5_KINEMATICS.hrTimeConsensusOdba) ;
% 
% e5_meanOdav = mean(LONO.E5_KINEMATICS.hrTimeConsensusOdav) ;
% e5_stdOdav  = std(LONO.E5_KINEMATICS.hrTimeConsensusOdav) ;
% 
% 
% e5_meanHr = mean(LONO.E5_KINEMATICS.hrTimeConsensus) ;
% e5_stdHr  = std(LONO.E5_KINEMATICS.hrTimeConsensus) ;

%%

clc;

fprintf('E1 odba mean (sd): %3.1f (±%3.1f)\n', e1_meanOdba, e1_stdOdba);
fprintf('E1 odav mean (sd): %3.1f (±%3.1f)\n', e1_meanOdav, e1_stdOdav);
fprintf('E1 consensus mean (sd): %3.1f (±%3.1f)\n', e1_meanHr, e1_stdHr);
fprintf('\n');

fprintf('E3 odba mean (sd): %3.1f (±%3.1f)\n', e3_meanOdba, e3_stdOdba);
fprintf('E3 odav mean (sd): %3.1f (±%3.1f)\n', e3_meanOdav, e3_stdOdav);
fprintf('E3 consensus mean (sd): %3.1f (±%3.1f)\n', e3_meanHr, e3_stdHr);
fprintf('\n');

% fprintf('E5 odba mean (sd): %3.1f (±%3.1f)\n', e5_meanOdba, e5_stdOdba);
% fprintf('E5 odav mean (sd): %3.1f (±%3.1f)\n', e5_meanOdav, e5_stdOdav);
% fprintf('E5 consensus mean (sd): %3.1f (±%3.1f)\n', e5_meanHr, e5_stdHr);
% fprintf('\n');

%%

% Step-by-step through ANOVA testing of each triplet

% Kolohe E1

[p, tbl, stats] = anova1( ...
    [ KOLOHE.E1_KINEMATICS.consensusHrOdba KOLOHE.E1_KINEMATICS.consensusHrOdav ] );
figure; [results, means ] = multcompare(stats, 'CType','hsd')
 
[p, tbl, stats] = anova1( [ KOLOHE.E1_KINEMATICS.hrTimeOdbaRidges ] );
figure; [results, means ] = multcompare(stats, 'CType','hsd')

[p, tbl, stats] = anova1( [ KOLOHE.E1_KINEMATICS.hrTimeOdavRidges ] );
figure; [results, means ] = multcompare(stats, 'CType','hsd')

sizeRidges = size(KOLOHE.E1_KINEMATICS.hrTimeOdbaRidges);
if (exist('ridges','var'))
    clear ridges;
end
ridges(sizeRidges(1),sizeRidges(2)*2) = zeros;
ridges(:,1:5) = KOLOHE.E1_KINEMATICS.hrTimeOdbaRidges;
ridges(:,6:10) = KOLOHE.E1_KINEMATICS.hrTimeOdavRidges;

[p, tbl, stats] = anova1( ridges  );
figure; [results, means ] = multcompare(stats, 'CType','hsd')


figure; [p, tbl, stats] = kruskalwallis( ridges );

figure; [results, means ] = multcompare(stats, 'CType','hsd')

% Kolohe E3

[p, tbl, stats] = anova1( [ KOLOHE.E3_KINEMATICS.consensusHrOdba KOLOHE.E3_KINEMATICS.bpOdav ] )
figure; [results, means ] = multcompare(stats, 'CType','hsd')

[p, tbl, stats] = anova1( [ KOLOHE.E3_KINEMATICS.hrTimeOdbaRidges ] );
figure; [results, means ] = multcompare(stats, 'CType','hsd')


if (exist('ridges','var'))
    clear ridges;
end
sizeRidges = size(KOLOHE.E3_KINEMATICS.hrTimeOdbaRidges);
ridges(sizeRidges(1),sizeRidges(2)*2) = zeros;
ridges(:,1:5) = KOLOHE.E3_KINEMATICS.hrTimeOdbaRidges;
ridges(:,6:10) = KOLOHE.E3_KINEMATICS.hrTimeOdavRidges;

[p, tbl, stats] = anova2( ridges  );
figure; [results, means ] = multcompare(stats, 'CType','hsd')

% Kolohe E5

[p, tbl, stats] = anova1( [ KOLOHE.E5_KINEMATICS.consensusHrOdba KOLOHE.E5_KINEMATICS.bpOdav ] )
figure; [results, means ] = multcompare(stats, 'CType','hsd')

[p, tbl, stats] = anova1( [ KOLOHE.E5_KINEMATICS.hrTimeOdbaRidges ] );
figure; [results, means ] = multcompare(stats, 'CType','hsd')

if (exist('ridges','var'))
    clear ridges;
end
sizeRidges = size(KOLOHE.E5_KINEMATICS.hrTimeOdbaRidges);
ridges(sizeRidges(1),sizeRidges(2)*2) = zeros;
ridges(:,1:5) = KOLOHE.E5_KINEMATICS.hrTimeOdbaRidges;
ridges(:,6:10) = KOLOHE.E5_KINEMATICS.hrTimeOdavRidges;

[p, tbl, stats] = anova2( ridges  );
figure; [results, means ] = multcompare(stats, 'CType','hsd')

%% Run EPOCHal ANOVA tests

% one-way ANOVA for E1 epoch

animalData = [ ones(size(HOKU.E1_KINEMATICS.hrTimeConsensus)); ...
    2 * ones(size(LIHO.E1_KINEMATICS.hrTimeConsensus)); ...
    3 * ones(size(KOLOHE.E1_KINEMATICS.hrTimeConsensus)); ...
    4 * ones(size(HUA.E1_KINEMATICS.hrTimeConsensus)); ...
    5 * ones(size(NOA.E1_KINEMATICS.hrTimeConsensus)); ...
    6 * ones(size(LONO.E1_KINEMATICS.hrTimeConsensus)) ];

[p, tbl, stats] = anova1( [ HOKU.E1_KINEMATICS.hrTimeConsensus ; ...
    LIHO.E1_KINEMATICS.hrTimeConsensus ; ...
    KOLOHE.E1_KINEMATICS.hrTimeConsensus ; ...
    HUA.E1_KINEMATICS.hrTimeConsensus ; ...
    NOA.E1_KINEMATICS.hrTimeConsensus ; ...
    LONO.E1_KINEMATICS.hrTimeConsensus  ], animalData );

figure; [results, means ] = multcompare(stats, 'CType','hsd')

% one-way ANOVA for E3 epoch

animalData = [ ones(size(HOKU.E3_KINEMATICS.hrTimeConsensus)); ...
    2 * ones(size(LIHO.E3_KINEMATICS.hrTimeConsensus)); ...
    3 * ones(size(KOLOHE.E3_KINEMATICS.hrTimeConsensus)); ...
    4 * ones(size(HUA.E3_KINEMATICS.hrTimeConsensus)); ...
    5 * ones(size(NOA.E3_KINEMATICS.hrTimeConsensus)); ...
    6 * ones(size(LONO.E3_KINEMATICS.hrTimeConsensus)) ];

[p, tbl, stats] = anova1( [ HOKU.E3_KINEMATICS.hrTimeConsensus ; ...
    LIHO.E3_KINEMATICS.hrTimeConsensus ; ...
    KOLOHE.E3_KINEMATICS.hrTimeConsensus ; ...
    HUA.E3_KINEMATICS.hrTimeConsensus ; ...
    NOA.E3_KINEMATICS.hrTimeConsensus ; ...
    LONO.E3_KINEMATICS.hrTimeConsensus  ], animalData );

figure; [results, means ] = multcompare(stats, 'CType','hsd')

% one-way ANOVA for E5 epoch (exclude Hoku and Lono, tags were off for E5)

animalData = [ ones(size(LIHO.E5_KINEMATICS.hrTimeConsensus)); ...
    2 * ones(size(KOLOHE.E5_KINEMATICS.hrTimeConsensus)); ...
    3 * ones(size(HUA.E5_KINEMATICS.hrTimeConsensus)); ...
    4 * ones(size(NOA.E5_KINEMATICS.hrTimeConsensus)) ];

[p, tbl, stats] = anova1( [ LIHO.E5_KINEMATICS.hrTimeConsensus ; ...
    KOLOHE.E5_KINEMATICS.hrTimeConsensus ; ...
    HUA.E5_KINEMATICS.hrTimeConsensus ; ...
    NOA.E5_KINEMATICS.hrTimeConsensus   ], animalData );

figure; [results, means ] = multcompare(stats, 'CType','hsd')

%% linear mixed-effects model

makeThis = size( [ LIHO.E5_KINEMATICS.hrTimeConsensus ; ...
    KOLOHE.E5_KINEMATICS.hrTimeConsensus ; ...
    HUA.E5_KINEMATICS.hrTimeConsensus ; ...
    NOA.E5_KINEMATICS.hrTimeConsensus ] );



lme = fitlme([ LIHO.E5_KINEMATICS.hrTimeConsensus ; ...
    KOLOHE.E5_KINEMATICS.hrTimeConsensus ; ...
    HUA.E5_KINEMATICS.hrTimeConsensus ; ...
    NOA.E5_KINEMATICS.hrTimeConsensus ], 



%%

figure('Color', white);
s1a = subplot(411);
plot(E1_KINEMATICS.Time, E1_KINEMATICS.sumImfBpOdba); 
hold on; 
plot(E1_KINEMATICS.Time,E1_KINEMATICS.bestImfBpOdba); 
plot(E2_KINEMATICS.Time, E2_KINEMATICS.sumImfBpOdba); 
plot(E2_KINEMATICS.Time, E2_KINEMATICS.bestImfBpOdba); 
plot(E3_KINEMATICS.Time, E3_KINEMATICS.sumImfBpOdba); 
plot(E3_KINEMATICS.Time, E3_KINEMATICS.bestImfBpOdba); 
plot(E4_KINEMATICS.Time, E4_KINEMATICS.sumImfBpOdba); 
plot(E4_KINEMATICS.Time, E4_KINEMATICS.bestImfBpOdba); 
plot(E5_KINEMATICS.Time, E5_KINEMATICS.sumImfBpOdba); 
plot(E5_KINEMATICS.Time, E5_KINEMATICS.bestImfBpOdba); 
hold off; 
grid;
s1b = subplot(412);
plot(E1_KINEMATICS.Time,E1_KINEMATICS.bestImfBpOdba);
hold on; 
plot(E2_KINEMATICS.Time, E2_KINEMATICS.bestImfBpOdba); 
plot(E3_KINEMATICS.Time, E3_KINEMATICS.bestImfBpOdba); 
plot(E4_KINEMATICS.Time, E4_KINEMATICS.bestImfBpOdba); 
plot(E5_KINEMATICS.Time, E5_KINEMATICS.bestImfBpOdba); 
hold off; 
grid;
s1c = subplot(413);
plot(E1_KINEMATICS.Time, E1_KINEMATICS.sumImfBpOdav); 
hold on; 
plot(E1_KINEMATICS.Time,E1_KINEMATICS.bestImfBpOdav); 
plot(E2_KINEMATICS.Time, E2_KINEMATICS.sumImfBpOdav); 
plot(E2_KINEMATICS.Time, E2_KINEMATICS.bestImfBpOdav); 
plot(E3_KINEMATICS.Time, E3_KINEMATICS.sumImfBpOdav); 
plot(E3_KINEMATICS.Time, E3_KINEMATICS.bestImfBpOdav); 
plot(E4_KINEMATICS.Time, E4_KINEMATICS.sumImfBpOdav); 
plot(E4_KINEMATICS.Time, E4_KINEMATICS.bestImfBpOdav); 
plot(E5_KINEMATICS.Time, E5_KINEMATICS.sumImfBpOdav); 
plot(E5_KINEMATICS.Time, E5_KINEMATICS.bestImfBpOdav); 
hold off; 
grid;
s1d = subplot(414);
plot(E1_KINEMATICS.Time,E1_KINEMATICS.bestImfBpOdav);
hold on; 
plot(E2_KINEMATICS.Time, E2_KINEMATICS.bestImfBpOdav); 
plot(E3_KINEMATICS.Time, E3_KINEMATICS.bestImfBpOdav); 
plot(E4_KINEMATICS.Time, E4_KINEMATICS.bestImfBpOdav); 
plot(E5_KINEMATICS.Time, E5_KINEMATICS.bestImfBpOdav); 
hold off; 
grid;
linkaxes([s1a s1b s1c s1d],'x');