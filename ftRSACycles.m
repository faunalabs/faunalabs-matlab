figure('Color', white');

kFs = 100;

time_s = 1:(1/kFs):60;

for vi = 1:numel(CC)-1
    
    startAt = CC.Time(vi);
    endAt   = CC.Time(vi+1) - seconds(2);
    
    thisDuration = seconds(endAt - startAt);
    
    v = timerange(CC.Time(vi), CC.Time(vi+1) - seconds(2));
    
    thisRsa = thisEpoch.consensusHr(v) * 60;
    
    lenThisRsa = length(thisRsa);
    meanThisRsa = mean(thisRsa);
    stdErrThisRsa = std(thisRsa) / sqrt(lenThisRsa);
    
    fprintf('RSA cycle %d - duration: %3.1f | mean(SE) iHR: %3.1f (%3.1f)\n\n', ...
        vi, thisDuration, meanThisRsa, stdErrThisRsa);
    
    
    hold on;
    
    plot(time_s(1:lenThisRsa), thisRsa);
    


end

hold on;
xlabel('Time, seconds');
ylabel('iHR');
grid;
hold off;