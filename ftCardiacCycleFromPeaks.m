iterations = numel(E1_PEAKS_ODBA.width);

before = E1_PEAKS_ODBA.width ./ 1000;

%%

maxlags = numel(linearVxyz) * 0.5;
[xc,lag] = xcov(linearVxyz(:,1), maxlags);
figure; plot(xc);
figure; plot(xc(df));
figure; plot(lag, xc, 'k-', lag(df), xc(df),'kv','MarkerFaceColor','r');
[~,df] = findpeaks(xc, 'MinPeakDistance',75);
figure; plot(lag, xc, 'k-', lag(df), xc(df),'kv','MarkerFaceColor','r');
[~,mf] = findpeaks(xc,'MinPeakDistance',37);
figure; plot(lag, xc, 'k-', lag(df), xc(df),'kv','MarkerFaceColor','r');
figure; plot(diff(df); grid;
figure; plot(diff(df)); grid;
figure; plot(diff(mf)); grid;
[alignVx,alignVy] = alignsignals(linearVxyz(:,1),linearVxyz(:,2));
figure; plot(alignVx);
figure; plot(alignVy);
doc amd
doc median
doc peak2rms
peak2rms(linearVxyz)
peak2peak(linearVxyz)
rms(linearVxyz)
1.4806 / 0.0860
figure; dtw(linearVxyz(:,1), linearVxyz(:,2))
figure; boxplot(linearVxyz); grid;
figure; midcross(linearVxyz(:,1),100,'tolerance',25);
figure; midcross(bpAx,100,'tolerance',25);
figure; midcross(linearVxyz(:,1),100,'tolerance',50); grid;
figure; midcross(linearVxyz(:,1),100,'tolerance',49); grid;
figure; midcross(linearVxyz(:,1),100,'tolerance',10); grid;
figure; midcross(linearVxyz(:,1),100,'tolerance',5); grid;
figure; midcross(bpAx(:,1),100,'tolerance',5); grid;
midxBpAx = midcross(bpAx,100,'tolerance',1);
size(midxBpAx)
figure; plot(midxBpAx);
figure; plot(thisT_s, midxBpAx);
midxBpAx(1:5)
midxBpAx(1:15)
figure; plot(thisT_s, bpAx);
figure; plot(downsample(thisT_s,4), decimate(bpAx,4) );
figure; plot(thisT_s, bpAx);
doc midcross
statelevels(bpAx)
figure; midcross(bpAx(:,1),100,'tolerance',5); grid;
figure; plot(downsample(thisT_s,4), downsample(bpAx,4) );
statelevels(bpAx,4)
figure; midcross( downsample(bpAx,4),25,'tolerance',5); grid;
figure; midcross( downsample(bpAx,4),25,'tolerance',25); grid;
figure; midcross(bpAx,100,'tolerance',25); grid;
figure; midcross(linearVxyz,100,'tolerance',25); grid;
figure; midcross(linearVxyz(:,2), 100, 'tolerance', 25); grid;
figure; midcross(linearVxyz(:,3), 100, 'tolerance', 25); grid;
figure; midcross(bpGy, 100, 'tolerance', 25); grid;
figure; midcross(diff(bpGy), 100, 'tolerance', 25); grid;
figure; midcross(diff(bpGy), 100, 'tolerance', 5); grid;
figure; midcross( iLinKE(:,1), 100, 'tolerance', 5); grid;
figure; midcross( iLinKE(:,1), 100, 'tolerance', 25); grid;
figure; midcross( iLinKE(:,2), 100, 'tolerance', 25); grid;
figure; midcross( iLinKE(:,3), 100, 'tolerance', 25); grid;
figure; midcross( iRotKE(:,1), 100, 'tolerance', 25); grid;
figure; midcross( iRotKE(:,2), 100, 'tolerance', 25); grid;
figure; midcross( iRotKE(:,4), 100, 'tolerance', 25); grid;
figure; midcross( iRotKE(:,3), 100, 'tolerance', 25); grid;
figure; midcross( linKer, 100, 'tolerance', 25); grid;
figure; midcross( iRotKE(:,3), 100, 'tolerance', 25); grid;
figure; midcross( iRotKE(:,3), thisT_s, 'tolerance', 25); grid;
figure; midcross( iRotKE(:,3), thisT_s(1:end-1), 'tolerance', 25); grid;
midcross( iRotKE(:,3), thisT_s(1:end-1), 'tolerance', 25); grid;
[cIRotKeGz,MIDLEV] = midcross( iRotKE(:,3), thisT_s(1:end-1), 'tolerance', 25);
cIRotKeGz
[cIRotKeGz,MIDLEV] = midcross( iRotKE(:,3), 100, 'tolerance', 25);
cIRotKeGz
0.403 * 0.001
0.403 * 0.1
rms(bpAx)
peak2peak(bpAx)
peak2peak(bpAx) / rms(bpAx)
peak2peak(bpAy) / rms(bpAy)
peak2peak(bpAz) / rms(bpAz)
snr(bpAx)
snr(bpAy)
snr(bpAz)
snr(bpAx(1000:1500))
snr(bpAx(2000:2500))
snr(bpAx(2000:3500))
snr(bpAx(2000:4000))
snr(bpAx(2000:5000))
snr(bpAx(100:1000))
snr(bpAx,100,7);
[outIkERotZ,MIDLEV] = midcross( iRotKE(:,3), 100, 'tolerance', 25);
[outIkERotY,MIDLEV] = midcross( iRotKE(:,2), 100, 'tolerance', 25);
[outIkERotX,MIDLEV] = midcross( iRotKE(:,1), 100, 'tolerance', 25);
outIkERotX
outIkERotY
outIkERotZ
oFs / 50
1/5
kFs
1/ (kFs / 10)
1/ (kFs / 15)
kWS = kFs / 15;
gfBk = ones(1, kWS) / kWS;
kGF = kFs / 15;
kGF
kGF = kFs / 12.5
1/8
kGF = kFs / 12.5;
gfBk = ones(1, kGF) / kGF;
gfA = [1];
gfBpAx = filtfilt(gfBk, gfA, bpAx);
figure; plot(gfBpAx);
figure; plot(thisT_s, gfBpAx);
figure; plot(thisT_s, gfBpAx, thisT_s, bpAx);
doc designfilt
local_max
which local_max
doc local_max
local_max(linearVxys(1000:1100))
local_max(linearVxyz(1000:1100))
doc maxima
figure; wsst(iRotKE(:,2), 100); shading interp;
figure; wsst(iRotKE(:,2), 100); shading interp; ylim([0 4]);
figure; wsst(iLinKE(:,1), 100); shading interp; ylim([0 4]);
[wsstThis, fThis] = wsst(iLinKE(:,1), 100);
[ridge, iridge] = wsstridge(wsstThis, 20, kFs);
figure; plot(ridge); grid;
[ridge, iridge] = wsstridge(wsstThis, 20, fThis);
figure; plot(ridge); grid;
[wsstThis, fThis] = wsst(iRotKE(:,2), 100);
[ridge, iridge] = wsstridge(wsstThis, 20, fThis);
figure; plot(ridge); grid;
iterations = numel(E1_PEAKS_ODBA.width);
iterations
before = E3_PEAKS_ODBA.width ./ 1000;
before
Ebefore = E1_PEAKS_ODBA.width ./ 1000;
before = E1_PEAKS_ODBA.width ./ 1000;
before
E1_PEAKS_ODBA.width
figure; plot(diff(bpGy); grid;
figure; plot(diff(bpGy)); grid;
figure; plot(diff(bpGx)); grid;
figure; plot(diff(bpGz)); grid;
figure; plot(linearVxyz); grid;
figure; plot(iRotKE); grid;
figure; findpeaks(iRotKE(:,1));
figure; findpeaks(iRotKE(:,1),65);
figure; findpeaks(iRotKE(:,1),kFs);
figure; findpeaks(iRotKE(:,1),350);
figure; findpeaks(iRotKE(:,1),'MinPeakDistance',65);
figure; findpeaks(iRotKE(:,2),'MinPeakDistance',65);
figure; findpeaks(iRotKE(:,2),'MinPeakDistance',65,'MinPeakHeight',0.7);
figure; findpeaks(iRotKE(:,2),'MinPeakDistance',65,'MinPeakHeight',0.0007);
median(iRotKE(:,2))
meann(iRotKE(:,2))
mean(iRotKE(:,2))
std(iRotKE(:,2))


%% 

lenKinematics = length(bpAx);       % figure out how long to run the loop
slidingWin = 500;                    % 50 samples = 500 ms on a 100 Hz rate

thisFig = figure('Color', white);

for i = 1:slidingWin:lenKinematics-50

    iEnd = i + 499;
    
    rmsBpAx = rms(bpAx(i:iEnd));
    rmsBpAy = rms(bpAy(i:iEnd));
    rmsBpAz = rms(bpAz(i:iEnd));

    rmsBpGx = rms(bpGx(i:iEnd));
    rmsBpGy = rms(bpGy(i:iEnd));
    rmsBpGz = rms(bpGz(i:iEnd));

    fprintf('\tRMS - Ax: %.4f | Ay: %.4f | Az: %.4f\n', rmsBpAx, rmsBpAy, rmsBpAz);
    fprintf('\tRMS - Gx: %.4f | Gy: %.4f | Gz: %.4f\n', rmsBpGx, rmsBpGy, rmsBpGz);

    rmsBpAxThresh = 3 * rmsBpAx;
    rmsBpAyThresh = 3 * rmsBpAy;
    rmsBpAzThresh = 3 * rmsBpAz;

    rmsBpGxThresh = 3 * rmsBpGx;
    rmsBpGyThresh = 3 * rmsBpGy;
    rmsBpGzThresh = 3 * rmsBpGz;
    
    figure(thisFig);
    
    pAx = subplot(321);
    plot(KINEMATICS.Time(i:iEnd), bpAx(i:iEnd), 'Color', blue);
    hold on; 
    yline(rmsBpAxThresh, '--', 'Color', blue);
    hold off;
    xlabel('Time, local');
    ylabel('m/s^2');
    title('Band-passed Ax');
    grid;

    pAy = subplot(323);
    plot(KINEMATICS.Time(i:iEnd), bpAy(i:iEnd), 'Color', red);
    hold on; 
    yline(rmsBpAyThresh, '--', 'Color', red);
    hold off;
    xlabel('Time, local');
    ylabel('m/s^2');
    title('Band-passed Ax');
    grid;

    pAz = subplot(325);
    plot(KINEMATICS.Time(i:iEnd), bpAz(i:iEnd), 'Color', goldenrod);
    hold on; 
    yline(rmsBpAzThresh, '--', 'Color', goldenrod);
    hold off;
    xlabel('Time, local');
    ylabel('m/s^2');
    title('Band-passed Ax');
    grid;

    pGx = subplot(322);
    plot(KINEMATICS.Time(i:iEnd), bpGx(i:iEnd), 'Color', blue);
    hold on; 
    yline(rmsBpGxThresh, '--', 'Color', blue);
    hold off;
    xlabel('Time, local');
    ylabel('°/s');
    title('Band-passed Gx');
    grid;

    pGy = subplot(324);
    plot(KINEMATICS.Time(i:iEnd), bpGy(i:iEnd), 'Color', red);
    hold on; 
    yline(rmsBpGyThresh, '--', 'Color', red);
    hold off;
    xlabel('Time, local');
    ylabel('°/s');
    title('Band-passed Gx');
    grid;

    pGz = subplot(326);
    plot(KINEMATICS.Time(i:iEnd), bpGz(i:iEnd), 'Color', goldenrod);
    hold on; 
    yline(rmsBpGzThresh, '--', 'Color', goldenrod);
    hold off;
    xlabel('Time, local');
    ylabel('°/s');
    title('Band-passed Gx');
    grid;

    linkaxes([pAx pAy pAz pGx pGy pGz],'x');
    
    fprintf('i = %d - click left mouse button to continue...\n', i);
    
    [~,~,b] = ginput(1);
    
    switch(b)
        case 'n'
            fprintf('Moving to next sliding window increment...\n');
        case 'q'
            fprintf('Caught a quit. Ending.\n');
            break;
        case '1'
            fprintf('Moving to next sliding window of data...\n');
        otherwise
            fprintf('Invalid input.\n');
    end
    
end




