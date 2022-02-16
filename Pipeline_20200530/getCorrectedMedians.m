%{ 
    It's really your choice to remove the trend or not
    
    INPUT: 
    paramVec - numTrials x numElectrodes matrix for one of the 3 measures
    polarAng - polar angles of the electrodes wrt reference
    toRemoveTrend - pass as 1 to remove quadratic trend

    OUTPUT: 
    medParamVec - Corrected medians for the parameter passed into this
                  function

    You can now use these corrected medians to compute the confidence
    interval and hence remove bad electrodes. 

    Author - SC, September 12th, 2017

%}
function medParamVecR = getCorrectedMedians(paramVec, polarAng, toRemoveTrend)

% A transformation I need to do for EGI
% Doesn't really matter as a trend removal is a trend removal
polarAng = abs(polarAng - nanmax(polarAng)); 

% Seems to be doing a median across trials
if size(paramVec,1) > 1
    medParamVec = nanmedian(paramVec, 1); 
else
    medParamVec = paramVec; % If there is only 1 epoch
end

% De-trending
if toRemoveTrend == 1
    coeffs = polyfit(polarAng(~isnan(medParamVec)), medParamVec(~isnan(medParamVec)), 2); 
    y_pred = coeffs(1).*power(polarAng, 2) + coeffs(2).*polarAng + coeffs(3); 
    medParamVecR = medParamVec - y_pred; 
end
% figure; 
% plot(polarAng, medParamVec, 'ko'); hold on; 
% text(polarAng, medParamVec, split(num2str(1:128))); 
% plot(polarAng, medParamVecR, 'ro'); 
% text(polarAng, medParamVecR, split(num2str(1:128))); 
%pause; 
% close all; 
end
