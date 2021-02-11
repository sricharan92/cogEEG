% Pass the median of the parameter across trials here
% 
% This is a small change to accomodate detrending as well as no detrending
function rejectSensors = cog_scads_1_3(paraMat, lambda, polarAng)
    medSensor = getCorrectedMedians(paraMat, polarAng, 1); 
    grandMed = nanmedian(medSensor);
    
    limit = lambda*nt_rms(medSensor - grandMed);
    lowerBound = grandMed - limit;
    upperBound = grandMed + limit;
    
    rejectSensors = find(medSensor < lowerBound | medSensor > upperBound);
end