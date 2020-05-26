function rejectTrials = cog_scads_1_4(paraMat, lambda, polarAng)
medSensor = getCorrectedMedians(paraMat, polarAng, 1);
grandMed = nanmedian(medSensor, 2);
% This loops over number of sensors
for i = 1:size(paraMat,1)
    limit(i) = lambda*rms(paraMat(i,:) - grandMed(i));
    lb(i) = grandMed(i) - limit(i);
    ub(i) = grandMed(i) + limit(i);
end
% Trials to reject for each sensor
for i = 1:size(paraMat,1)
    rejectTrials{i, 1} = find(paraMat(i, :) < lb(i) | paraMat(i, :) > ub(i));
end
%     keyboard;
end