%% ftBandpassWsstridgeFinder.m

% filtLed2 and filtLed were decimated in the first example I did, but it
% may be worth trying this with full resolution (250 Hz) data, then compare
% it with the decimated by 5 (50 Hz) data.

%% I. start from scratch using a baseline 50 Hz data set

% for optics, that's decimate by a factor of 5

o50Hz = 5;
dLed1 = decimate(filtLed1, o50Hz);
dLed2 = decimate(filtLed2, o50Hz);
dLed3 = decimate(filtLed3, o50Hz);
dLed4 = decimate(filtLed4, o50Hz);

% for kinematic signals in E, decimate by a factor of 2

k50Hz = 2;
dBpOdba = decimate(E.bpOdba, k50Hz);
dBpOdav = decimate(E.bpOdav, k50Hz);
dBestImfBpOdba = decimate(E.bestImfBpOdba, k50Hz);
dBestImfBpOdav = decimate(E.bestImfBpOdav, k50Hz);
dConsensusHr = decimate(E.consensusHr, k50Hz);

dT = downsample(oT, o50Hz);
dT_s = downsample(oT_s, o50Hz);

%% II. Create first derivatives

% take first derivative to make it peakier

diffLed2 = diff(dLed2 * -1);
diffLed3 = diff(dLed3 * -1);

% pad a zero at last sample for diff loss by one

diffLed2(end+1) = 0;
diffLed3(end+1) = 0;

% have a look at the diffLeds

figure;
subplot(211);
plot(dT, diffLed2, 'Color', blue);
subplot(212);
plot(dT, diffLed3, 'Color', red);
hold off;
grid;


%% III. Compute some VMD IMFs for the two active LEDs

tic
imfDiffLed2 = vmd(diffLed2, 'NumIMFs', 10);
imfDiffLed3 = vmd(diffLed3, 'NumIMFs', 10);
toc

% make the waveform that exposes heart rate peaks (IMF 1 worked well, but
% it's probably worth trying sumImf instead of guessing that bestImf is 1

waveLed2 = smooth(abs(imfDiffLed2(:,1)));
waveLed3 = smooth(abs(imfDiffLed3(:,1)));

% look at the waveLed* plot

figWavePlots = figure('Position',[100 100 800 800], 'Color', white);
figWavePlots.Name = 'Wave plots for active optical channels';
s1a = subplot(211);
plot(dT, waveLed2, 'Color', blue);
xlabel('Time, local');
ylabel('A.U.');
title('Waveform or optical IMF LED2');
grid;
s1b = subplot(212);
plot(dT, waveLed3, 'Color', red);
xlabel('Time, local');
ylabel('A.U.');
title('Waveform or optical IMF LED2');
grid;
linkaxes([s1a s1b],'x');

% make the wssts for the two waveForms

[wsstWaveLed2, fWaveLed2] = wsst(waveLed2, 50);
[wsstWaveLed3, fWaveLed3] = wsst(waveLed3, 50);

% do the simple band-pass of wavelet results, screening in 24 - 210 BPM

fWaveLed2HrRange = find(fWaveLed2 < 3.5 & fWaveLed2 > 0.6);
fWaveLed3HrRange = find(fWaveLed3 < 3.5 & fWaveLed3 > 0.6);

% we want the first and last element to bound wsstridge searches

fLed2Start = fWaveLed2HrRange(1);
fLed2End   = fWaveLed2HrRange(end);

fLed3Start = fWaveLed3HrRange(1);
fLed3End   = fWaveLed3HrRange(end);

% now try finding some ridges based on these bounds

[ ridgeWaveLed2, iridgeWaveLed2 ] = wsstridge( ...
    wsstWaveLed2(fLed2Start:fLed2End,:), 20, ...
    fWaveLed2(1,fLed2Start:fLed2End), 'NumRidges', 2);

[ ridgeWaveLed3, iridgeWaveLed3 ] = wsstridge( ...
    wsstWaveLed3(fLed3Start:fLed3End,:), 20, ...
    fWaveLed3(1,fLed3Start:fLed3End), 'NumRidges', 2);

% preview these ridges

figure; 
subplot(211)
plot(ridgeWaveLed2);
grid;
subplot(212);
plot(ridgeWaveLed3);
grid;


%% do the summImf, bestImf, thresh & peaks sqwave shuffles w/ abs(diffLed*)

% make absDiffLed*

absDiffLed2 = abs(diffLed2);
absDiffLed3 = abs(diffLed3);

% get 10 IMFs of VMDs for absDiffLed*

numImfs = 10;
numRows = numImfs / 2;
numCols = numImfs / numRows;

fprintf('Computing VMDs for abs(diffLed*)...\n');

tic

[ imfAbsDiffLed2, residualsAbsDiffLed2, infoAbsDiffLed2 ] = vmd( ...
    absDiffLed2, 'NumIMFs', numImfs);

[ imfAbsDiffLed3, residualsAbsDiffLed3, infoAbsDiffLed3 ] = vmd( ...
    absDiffLed3, 'NumIMFs', numImfs);

toc

fprintf('Completed compute of VMDs for abs(diffLed*).\n');

% plot the IMFs for absDiffLed2 and absDiffLed3

figAbsDiffLed2Imfs = figure('Color', white);

t = tiledlayout(numRows,numCols,'TileSpacing','compact','Padding','compact');
for n = 1:numImfs
    ax(n) = nexttile(t);
    plot(dT,imfAbsDiffLed2(:,n)')
    xlim([dT(1) dT(end)])
    txt = ['IMF',num2str(n)];
    title(txt)
    xlabel('Time, local')
end
title(t,'Variational Mode Decomposition - 1050 nm optical channel (absDiffLed2)');


figAbsDiffLed3Imfs = figure('Color', white);

t = tiledlayout(numRows,numCols,'TileSpacing','compact','Padding','compact');
for n = 1:numImfs
    ax(n) = nexttile(t);
    plot(dT,imfAbsDiffLed3(:,n)');
    xlim([dT(1) dT(end)]);
    txt = ['IMF',num2str(n)];
    title(txt)
    xlabel('Time, local')
end
title(t,'Variational Mode Decomposition - 1200 nm optical channel (absDiffLed3)');

sumImfAbsDiffLed2  = sum(imfAbsDiffLed2(:,2:9),2);
bestImfAbsDiffLed2 = imfAbsDiffLed2(:,9);

sumImfAbsDiffLed3  = sum(imfAbsDiffLed3(:,2:9),2);
bestImfAbsDiffLed3 = imfAbsDiffLed3(:,9);

dFs = 50;

[wsstSumImfLed2, fSumImfLed2] = wsst(sumImfAbsDiffLed2, dFs);
[wsstBestImfLed2, fBestImfLed2] = wsst(bestImfAbsDiffLed2, dFs);

[wsstSumImfLed3, fSumImfLed3] = wsst(sumImfAbsDiffLed3, dFs);
[wsstBestImfLed3, fBestImfLed3] = wsst(bestImfAbsDiffLed3, dFs);

% have a quick look

figWsstLed2 = figure('Color', white);
figWsstLed.Name = '1050 nm optical channel';
subplot(211);
pcolor(dT, fSumImfLed2, abs(wsstSumImfLed2));
shading interp;
ylim([0 4]);
subplot(212);
pcolor(dT, fBestImfLed3, abs(wsstBestImfLed3));
shading interp;
ylim([0 4]);

figWsstLed3 = figure('Color', white);
figWsstLed.Name = '1200 nm optical channel';
subplot(211);
pcolor(dT, fSumImfLed3, abs(wsstSumImfLed3));
shading interp;
ylim([0 4]);
subplot(212);
pcolor(dT, fBestImfLed3, abs(wsstBestImfLed3));
shading interp;
ylim([0 4]);

% find ridges for sumImf and bestImf LED2

[ridgeSumImfLed2, iridgeSumImfLed2] = wsstridge( wsstSumImfLed2, 20, ...
    fSumImfLed2, 'NumRidges', 2);
[ridgeBestImfLed2, iridgeBestImfLed2] = wsstridge( wsstBestImfLed2, 20, ...
    fBestImfLed2, 'NumRidges', 2);

% find ridges for sumImf and bestImf LED3

[ridgeSumImfLed3, iridgeSumImfLed3] = wsstridge( wsstSumImfLed3, 20, ...
    fSumImfLed3, 'NumRidges', 2);
[ridgeBestImfLed3, iridgeBestImfLed3] = wsstridge( wsstBestImfLed3, 20, ...
    fBestImfLed3, 'NumRidges', 2);

%% Ridge plots

% 1050 nm

figure;
subplot(211);
plot(dT, ridgeSumImfLed2);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('1050nm SumImf');
grid;
legend('r1','r2');

subplot(212);
plot(dT, ridgeBestImfLed2);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('1050nm BestImf');
grid;
legend('r1','r2');

% 1200 nm

figure;
subplot(211);
plot(dT, ridgeSumImfLed3);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('1200nm SumImf');
grid;
legend('r1','r2');

subplot(212);
plot(dT, ridgeBestImfLed3);
xlabel('Time, local');
ylabel('Frequency, Hz');
title('1200nm BestImf');
grid;
legend('r1','r2');



figCorr = figure('Position', [100 100 600 600], 'Color', white);

[r, pVal, figCorr ] = corrplot( ...
    [ dConsensusHr , ridgeBestImfLed2(:,1) , ridgeBestImfLed3(:,1) ] , ...
    'type', 'Pearson', 'testR', 'on');

%%

minDiffLed2 = min(diffLed2);
maxDiffLed2 = max(diffLed2);
rangeDiffLed2 = round(maxDiffLed2 - minDiffLed2);

DiffLed2 = rescale(diffLed2, 1, rangeDiffLed2+1);




