function [comp, dataSmallEpoch] = callICA(dataSmallEpoch)

% High pass filter the data from 1 Hz -- reference?
cfgF.hpfilter = 'yes'; cfgF.hpfreq = 1; cfgF.demean = 'yes'; cfgF.hpfilttype = 'but';
cfgF.hpfiltord = 4;
dataFilt = ft_preprocessing(cfgF, dataSmallEpoch);

% redefine the trials? 
for tr = 1:length(dataSmallEpoch.trial)
    indices = (dataSmallEpoch.trlSamples{tr}(1, 1):dataSmallEpoch.trlSamples{tr}(end, 2)) - dataSmallEpoch.samples{tr}(1, 1) + 1; 
    dataSmallEpoch.newTrlSamples{tr}(:, 1:2) = dataSmallEpoch.trlSamples{tr}(:, 1:2) - dataSmallEpoch.trlSamples{tr}(1, 1) + 1; 
    dataFilt.trial{tr} = dataFilt.trial{tr}(:, indices);
    dataFilt.time{tr} = dataFilt.time{tr}(:, indices); 
end

for tr = 1:length(dataSmallEpoch.trial)
    fprintf('Running ICA on epoch %d \n', tr); 
    cfgICA = [];
    cfgICA.trials = tr;
    comp(tr) = ft_componentanalysis(cfgICA, dataFilt);
end
% Save the ICA components
end
