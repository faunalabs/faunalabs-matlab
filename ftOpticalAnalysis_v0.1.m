%% FaunaTag BioOptics Analysis Playground
%  by Dave Haas
%  Based on ppgPlayground by Dave Haas, originally written on 8 May 2017
%  This version first ported for use with FaunaTag data sets on 1 May 2021
%  Dependency: csvproc from dtagtools written by Mark Johnson

clear;

%% Import tagOS CSV data (assumes the data has been 
%  previously processed with the tagOS_dataPrep tool...

% specify the filename here:
fileName = '/Users/dave/daveData_20210425.txt';

S = csvproc(fileName, []);

%% parse S into usable data components

% optics (relative insensity)           @ 250 Hz
led2 = str2double(S(1,:));
led3 = str2double(S(2,:));
led1 = str2double(S(3,:));
led4 = str2double(S(4,:));

% accel (meters per second squared)     @ 100 Hz
ax = str2double(S(5,:));
ay = str2double(S(6,:));
az = str2double(S(7,:));

% gyro (degrees per second)                     @ 100 Hz
gx = str2double(S(8,:));
gy = str2double(S(9,:));
gz = str2double(S(10,:));

% prh (degree)                                  @ 100 Hz
pitch   = str2double(S(11,:));
roll    = str2double(S(12,:));
heading = str2double(S(13,:));

% depth (meters)                                @ 100 Hz
depth = str2double(S(14,:));

% temp (degrees Celsius)                        @ 100 Hz
temperature = str2double(S(15,:));

% mag (degree)                                  @ 20 Hz
mx = str2double(S(16,:));
my = str2double(S(17,:));
mz = str2double(S(18,:));

% timestamp                                     @ 2 Hz
sampletime = str2double(S(19,:));

% power measurements                            @ 2 Hz
voltage             = str2double(S(20,:));
current             = str2double(S(21,:));
power               = str2double(S(22,:));
chargeState         = str2double(S(23,:));
remainingCapacity   = str2double(S(24,:));

%% Diagnostic optical plots

figure;

led1Plot = subplot(411);
plot(led1, 'k-');
title('LED1');
ylabel('Intensity');

led2Plot = subplot(412);
plot(led2, 'b-');
title('LED2');
ylabel('Intensity');

led3Plot = subplot(413);
plot(led3, 'g-');
title('LED3');
ylabel('Intensity');

led4Plot = subplot(414);
plot(led4, 'r-');
title('LED4');
ylabel('Intensity');

linkaxes([led1Plot led2Plot led3Plot led4Plot], 'x');


%% Invert raw PPG waveforms

%  what was I trying to do here once? 


%% get 1st and 2nd derivatives of the raw G, R  and IR PPG waveforms

d1Led1             = diff(led1(1,:), 1);
d1Led1(end + 1, 1) = 0;

d1Led2               = diff(led2(1,:), 1);
d1Led2(end + 1, 1)   = 0;

d1Led3               = diff(led3(1,:), 1);
d1Led3(end + 1, 1)   = 0;

d1Led4               = diff(led4(1,:), 1);
d1Led4(end + 1, 1)   = 0;

d2Led1             = diff(led1(1,:), 2);
d2Led1(end + 1, 1) = 0;
d2Led1(end + 1, 1) = 0;

d2Led2             = diff(led2(1,:), 2);
d2Led2(end + 1, 1) = 0;
d2Led2(end + 1, 1) = 0;

d2Led3             = diff(led3(1,:), 2);
d2Led3(end + 1, 1) = 0;
d2Led3(end + 1, 1) = 0;

d2Led4             = diff(led4(1,:), 2);
d2Led4(end + 1, 1) = 0;
d2Led4(end + 1, 1) = 0;



%% plot 1st ppg erivatives

figure; 
fig2a = subplot(411);
plot(d1Led1(1,:), 'g-');
ylabel('LED1 Intensity');
title('1st Derivative PPG Waveforms');
fig2b = subplot(412);
plot(d1Led2(1,:), 'r-');
ylabel('LED2 Intensity');
fig2c = subplot(413);
plot(d1Led3(1,:), 'b-');
ylabel('LED3 Intensity');
fig2d = subplot(414);
plot(d1Led4(1,:),'k-');
xlabel('Time, ms');
ylabel('LED4 Intensity');

linkaxes( ([fig2a, fig2b, fig2c, fig2d]), 'x');

%% plot 2nd ppg derivatives

figure; 
fig2a = subplot(411);
plot(d2Led1(1,:), 'g-');
ylabel('LED1 Intensity');
title('2nd Derivative PPG Waveforms');
fig2b = subplot(412);
plot(d2Led2(1,:), 'r-');
ylabel('LED2 Intensity');
fig2c = subplot(413);
plot(d2Led3(1,:), 'b-');
ylabel('LED3 Intensity');
fig2d = subplot(414);
plot(d2Led4(1,:),'k-');
xlabel('Time, ms');
ylabel('LED4 Intensity');

linkaxes( ([fig2a, fig2b, fig2c, fig2d]), 'x');

%% plot LED2 and its 1st and 2nd derivatives

figure; 
fig4a = subplot(311);
plot(led2, 'b-');
fig4b = subplot(312);
plot(d1Led2(1,:), 'k-');
fig4c = subplot(313);
plot(d2Led2(1,:),'m-');
linkaxes( ([fig4a, fig4b, fig4c]), 'x');
