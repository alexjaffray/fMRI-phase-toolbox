function respPhase = calculateRespPhase(respTrace)

rmax = max(respTrace) - min(respTrace);

respTrace = (respTrace - min(respTrace))/rmax;

[counts, edges, bins] = histcounts(respTrace,100);

H0 = sum(counts);
Hb = cumsum(counts);

diffTrace = diff(respTrace);
diffTrace = [diffTrace; diffTrace(end)];

Hb = zeros(size(respTrace));

for i = 1:size(respTrace)
    
    stopbin = bins(i);
    
    Hb(i) = sum(counts(1:stopbin));
    
end

respPhase = pi*Hb.*sign(diffTrace)/H0;

end

