% Returns a nElec array of editing matrices for samples x elecs
% Returns a nElec X nTrials (or nTrials X nElec?) matrix of editing matrices. 
% The medians reuqired are calculated in later cog functions. 
function [pMaxAmp, pStdDev, pGradient] = cog_scads_1_2(data3d)
    pMaxAmp = []; pStdDev = []; pGradient = [];
    if iscell(data3d)
        for i = 1:length(data3d)
            % THESE TAKE THE AVERAGE ACROSS numSamples
            pMaxAmp(:,i) = nanmax(abs(data3d{i}), [], 1);
            pStdDev(:,i) = nanstd(data3d{i}, 0, 1);
            pGradient(:,i) = nanmax(abs(gradient(data3d{i})), [], 1);
        end
    else
        pMaxAmp = squeeze(nanmax(abs(data3d), [], 1))';
        pStdDev = squeeze(nanstd(data3d, 0, 1))';
        pGradient = squeeze(nanmax(abs(gradient(data3d)), [], 1))';
    end
end
