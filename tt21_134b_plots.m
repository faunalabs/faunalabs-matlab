%% optical density plots for Kolohe trial - tt21_134b

load('/Users/dave/Desktop/tt21_134b.mat');

blue = [0 0.4470 0.7410];           % '#0072BD';
red = [0.8500 0.3250 0.0980];       % '#D95319'	
goldenrod = [0.9290 0.6940 0.1250];	% '#EDB120'	
purple = [0.4940 0.1840 0.5560];    % '#7E2F8E'	
green = [0.4660 0.6740 0.1880];     % '#77AC30'	
cyan = [0.3010 0.7450 0.9330];      % '#4DBEEE'	
maroon = [0.6350 0.0780 0.1840];    % '#A2142F'
black = [0 0 0];                    % '#000000'
white = [1 1 1];                    % '#FFFFFF'
grey = [0.6 0.6 0.6];               % '#999999'


%% convert raw optical channel data to OD

od1050_e1 = convertRaw2OD(KOLOHE_E1_DECIMATED.led2);
od1200_e1 = convertRaw2OD(KOLOHE_E1_DECIMATED.led3);
odAmb_e1  = convertRaw2OD(KOLOHE_E1_DECIMATED.led1);

od1050_e3 = convertRaw2OD(KOLOHE_E3_DECIMATED.led2);
od1200_e3 = convertRaw2OD(KOLOHE_E3_DECIMATED.led3);
odAmb_e3  = convertRaw2OD(KOLOHE_E3_DECIMATED.led1);

od1050_e5 = convertRaw2OD(KOLOHE_E5_DECIMATED.led2);
od1200_e5 = convertRaw2OD(KOLOHE_E5_DECIMATED.led3);
odAmb_e5  = convertRaw2OD(KOLOHE_E5_DECIMATED.led1);


%% E1 (baseline) OD plot w/ vent times and wsstHR

ventTimes = KOLOHE_E1_VENT_DATA.Time;

figure('Color', white);
s1a = subplot(211);
plot(KOLOHE_E1_DECIMATED.Time, od1050_e1, 'Color', blue, 'LineStyle', '-');
hold on;
plot(KOLOHE_E1_DECIMATED.Time, od1200_e1, 'Color', red, 'LineStyle', '-');
plot(KOLOHE_E1_DECIMATED.Time, odAmb_e1, 'Color', black, 'LineStyle', '-');
xline(ventTimes, 'Color', green, 'LineStyle', '--');
hold off;
xlabel('Time, local');
ylabel('OD');
title('Optical density - baseline (E1) epoch');
grid;
legend('1050 nm','1200 nm','ambient');

s1b = subplot(212);
plot(KOLOHE_E1_DECIMATED.Time, KOLOHE_E1_DECIMATED.wsstHr .* 60, ...
    'Color', goldenrod, 'LineWidth', 2);
hold on;
xline(ventTimes, 'Color', green, 'LineStyle', '--', 'LineWidth', 2);
hold off;
xlabel('Time, local');
ylabel('wsstHR, bpm');
ylim([0 120]);
title('Heart rate - Baseline (E1) epoch');
grid;


%% E3 (apnea) OD plot

figure('Color', white);
s3a = subplot(211);
plot(KOLOHE_E3_DECIMATED.Time, od1050_e3, 'Color', blue, 'LineStyle', '-');
hold on;
plot(KOLOHE_E3_DECIMATED.Time, od1200_e3, 'Color', red, 'LineStyle', '-');
plot(KOLOHE_E3_DECIMATED.Time, odAmb_e3, 'Color', black, 'LineStyle', '-');
hold off;
xlabel('Time, local');
ylabel('OD');
title('Optical density - apnea (E3) epoch');
grid;
legend('1050 nm','1200 nm','ambient');

s3b = subplot(212);
plot(KOLOHE_E3_DECIMATED.Time, KOLOHE_E3_DECIMATED.wsstHr .* 60, ...
    'Color', goldenrod,'LineWidth', 2);
xlabel('Time, local');
ylabel('wsstfH, beats · min^-^1');
ylim([0 120]);
grid;
title('Heart rate - apnea (E3) epoch');


%% E5 (recovery) OD plot

ventTimes = KOLOHE_E5_VENT_DATA.Time;

figure('Color', white);
s3a = subplot(211);
plot(KOLOHE_E5_DECIMATED.Time, od1050_e5, 'Color', blue, 'LineStyle', '-');
hold on;
plot(KOLOHE_E5_DECIMATED.Time, od1200_e5, 'Color', red, 'LineStyle', '-');
plot(KOLOHE_E5_DECIMATED.Time, odAmb_e5, 'Color', black, 'LineStyle', '-');
xline(ventTimes, 'Color', green, 'LineStyle', '--', 'LineWidth',2);
hold off;
xlabel('Time, local');
ylabel('OD');
title('Optical density - Recovery (E5) epoch');
grid;
legend('1050 nm','1200 nm','ambient','vents', 'Location','southeast');

s5b = subplot(212);
plot(KOLOHE_E5_DECIMATED.Time, KOLOHE_E5_DECIMATED.wsstHr .* 60, ...
    'Color', goldenrod, 'LineWidth', 2);
hold on;
xline(ventTimes, 'Color', green, 'LineStyle', '--', 'LineWidth',2);
hold off;
xlabel('Time, local');
ylabel('wsstHR, bpm');
ylim([0 120]);
title('Heart rate - Recovery (E5) epoch');
grid;
legend('wsstHR','vents','Location','southeast');


