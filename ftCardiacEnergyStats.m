%% ftCardiacEnergyStats.m
%
%   Written by Dave Haas
%   18 September 2021

figure('Color', white, 'Position', [0 400 800 500] );

sCCa = subplot(211);
plot(E1_KINEMATICS.Time, E1_KINEMATICS.hrTimeConsensus,'Color',green);
hold on;
plot(E3_KINEMATICS.Time, E3_KINEMATICS.hrTimeConsensus, 'Color', blue);
plot(E5_KINEMATICS.Time, E5_KINEMATICS.hrTimeConsensus, 'Color', red);
plot(E1_CC_ENERGY.Time, E1_CC_ENERGY.iHR, 'o', 'Color', green);
plot(E3_CC_ENERGY.Time, E3_CC_ENERGY.iHR, 'o', 'Color', blue);
plot(E5_CC_ENERGY.Time, E5_CC_ENERGY.iHR, 'o', 'Color', red);
xline(E1_VENT_DATA.Time, '-', 'Color', green);
xline(E5_VENT_DATA.Time, '-', 'Color', red);
hold off;
xlabel('Time, local');
ylabel('Heart rate, beat · min^-^1');
titleTxt = sprintf('Heart rate profile for %s', TRIAL.SUBJECT.id);
title(ti
ylim([0 150]);
grid; 
legend('baseline wsstHR','apnea wsstHR','recovery wsstHR', ...
    'baseline vents','recovery vents');

sCCb = subplot(212);
plot(E1_CC_ENERGY.Time, E1_CC_ENERGY.iTotalKE, 'x', 'Color', green);
hold on;
plot(E3_CC_ENERGY.Time, E3_CC_ENERGY.iTotalKE, 'x', 'Color', blue);
plot(E5_CC_ENERGY.Time, E5_CC_ENERGY.iTotalKE, 'x', 'Color', red);
xline(E1_VENT_DATA.Time, '--', 'Color', green);
xline(E5_VENT_DATA.Time, '--', 'Color', red);
hold on;
xlabel('Time, local');
ylabel('J · s');
grid;
legend('baseline','apnea','recovery', ...
    'baseline vents','recovery vents');

linkaxes([sCCa sCCb],'x');

%%

figure('Color', white, 'Position', [0 600 1200 500] );
plot(E1_CC_ENERGY.Time, E1_CC_ENERGY.iTotalKE, 'x', 'Color', green);
hold on;
plot(E3_CC_ENERGY.Time, E3_CC_ENERGY.iTotalKE, 'x', 'Color', blue);
plot(E5_CC_ENERGY.Time, E5_CC_ENERGY.iTotalKE, 'x', 'Color', red);
xline(E1_VENT_DATA.Time, '-', 'Color', green);
xline(E5_VENT_DATA.Time, '-', 'Color', red);
hold off;
xlabel('Time, seconds');
ylabel('J · s');
titleTxt = 'Change in cardiac energy during pre-apnea (baseline), apnea, and post-apnea (recovery)';
title(titleTxt);
grid;
legend('baseline','apnea','recovery','baseline vents','recovery vents');