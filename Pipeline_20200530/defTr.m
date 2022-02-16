function [cfgTr, cfgLoc] = defTr(info, fileIdx, data)

%% Fieldtrip epoching
eventTrial = 40; eventLocalizer = 230;

if ismember(info.subjectNumber, [15 16, 18:24 58]) % member of 4ADC-reduced Selection cohort
    eventTrial = 80; eventLocalizer = []; % No SSVEP localizer run
end

cfg = [];
cfg.dataset = sprintf('%s/%s', info.files(fileIdx).folder, info.sortedFileNames{fileIdx});
cfg.trialdef.eventtype  = 'STATUS';
cfg.trialdef.eventvalue = eventTrial;
cfg.trialdef.prestim    = 1;
cfg.trialdef.poststim   = 2;
cfgTr = ft_definetrial(cfg);
cfgTr.trl_rs = [round(cfgTr.trl(:,1:3) * (512/4096)), cfgTr.trl(:,4)];
cfgTr.trl_rs( (cfgTr.trl_rs(:,1) <= 0) ,1) = 1;
cfgTr.trl_rs( (cfgTr.trl_rs(:,3) > size(data.trial{1},2)) ,1) = size(data.trial{1},2);
cfgTr.trl_rs(:,5) = fileIdx; 

cfg.trialdef.eventvalue = eventLocalizer;
cfg.trialdef.prestim    = 0.5;
cfg.trialdef.poststim   = 5.5;

cfgLoc = []; % No localizer session run
if ((fileIdx == 1) || (fileIdx == 10)) && (~ismember(info.subjectNumber, [15 16, 18:24 58]))
    % Localizer was run
    cfgLoc = ft_definetrial(cfg);
    cfgLoc.trl_rs = [round(cfgLoc.trl(:,1:3) * (512/4096)), cfgLoc.trl(:,4)];
    cfgLoc.trl_rs( (cfgLoc.trl_rs(:,1) <= 0) ,1) = 1;
    cfgLoc.trl_rs( (cfgLoc.trl_rs(:,3) > size(data.trial{1},2)) ,1) = size(data.trial{1},2);
    cfgLoc.trl_rs(:,5) = fileIdx; 
end

end