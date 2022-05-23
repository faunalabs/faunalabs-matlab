%% ftWsstRidgeStitching.m


%% stitch wsstridge

testO3 = rLed3;

figStitch = figure('Position',[100 500 1000 500]);

maxO = round( length(testO3), 3, 'significant') - oFs; 

lastSample = maxO - oFs + 1;

numRidges = 4;

ridgeTestO3(1:lastSample,numRidges) = zeros;

% advance in 5 second windows

windowSize = 10;
windowLength = windowSize * oFs;

for i = 1:windowLength:maxO-windowLength

    startSample = i;
    endSample = i + windowLength;
    
    fprintf('startSample: %d | endSample: %d\n', ...
        startSample, endSample );

    sample = testO3(startSample:endSample,1);
    
    [wsstThis, fThis] = wsst(sample, oFs);
    
    [thisRidge, thisIridge] = wsstridge(wsstThis, 20, ...
        fThis, 'NumRidges', numRidges);
    
    ridgeTestO3(startSample:endSample,1:numRidges) = ...
        thisRidge(:,1:numRidges);
    %ridgeTestO3(startSample:endSample,2) = thisRidge(2,:)';
    %ridgeTestO3(startSample:endSample,3) = thisRidge(3,:)';
    %ridgeTestO3(startSample:endSample,3) = thisRidge(3,:)';
    
    figure(figStitch);
    plot(oT_s(startSample:endSample), thisRidge);
    grid;
    
%     fprintf('[c] to continue to next sample, [q] to quit\n');
%     [~,~,button] = ginput(1);
%     if (button == 'q')
%         fprintf('Caught a [q] - quitting now!\n');
%         return;
%     end
        
end


figure;
plot(ridgeTestO3);
grid;




