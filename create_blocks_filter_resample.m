function [dat, trl, trlInBetween] = create_blocks_filter_resample(dataset, events, eventsInBetween, bp_freq, resample_fs, info)

% Information
s           = info.subject;
session     = info.session;
results_dir = info.results_dir;
block       = info.block; 

% ------------
flag128 = 1; 
dispEvts = 1; 
% ------------

% Find pause indices
disp('Loading events')
cfg1 = [];
cfg1.dataset = dataset;
cfg1.trialdef.eventtype = '_DINs'; 
cfg1.trialdef.eventvalue = {events{1}, events{2}};

% Check if events file already exists
file_name = ['events' num2str(s, '%.2d') '_session_' num2str(session, '%.2d')];
if ~exist([results_dir file_name '.mat'], 'file')
    disp(['Saving events as ' file_name '.'])
    
    % ------------------
    % If you get a 'could not define trials' error, and you're on EGI, 
    % change cfg1.trialdef.eventtype to '_DINs'
    cfgTr_ref = ft_definetrial(cfg1); % Loads events, takes time
    % ------------------
    
    % Shows you all the unique DIN events and their counts
    if dispEvts == 1
        listEventsEGI(cfgTr_ref.event, cfg1.trialdef.eventtype); 
    end
    
    % Converts all events to the 127 system 
    if flag128 == 1
        % Change all the events to 127 base
        cfgTr_ref.event = changeEvtBase(cfgTr_ref.event, cfg1.trialdef.eventtype);  
    end
    
    save([results_dir file_name], 'cfgTr_ref');
else
    disp([file_name ' already exists. Loading.'])
    load([results_dir file_name], 'cfgTr_ref')
end

% Get epoch indices
epochInd = []; epCount = 1; midEpCount = 1; 
for ev = 1:length(cfgTr_ref.event)
    if strcmp(cfgTr_ref.event(ev).type, 'epoch')
        epochInd(epCount) = ev; 
        epCount = epCount + 1; 
    end
    
    if strcmp(cfgTr_ref.event(ev).type, '_DINs')
        if cfgTr_ref.event(ev).sample - cfgTr_ref.event(ev-1).sample > 10000
            midEpInd(midEpCount) = ev;
            midEpCount = midEpCount + 1; 
        end
    end
    
end

if s == 1 && block == 1
    evtToTake = [epochInd(1):midEpInd(2)-1, epochInd(3):epochInd(4), midEpInd(4):length(cfgTr_ref.event)]; 
elseif s == 1 && block == 2
    evtToTake = [epochInd(1):midEpInd(2)-1, epochInd(3):length(cfgTr_ref.event)]; 
elseif s == 5 && block == 2
    evtToTake = [midEpInd(1):midEpInd(2)-1, epochInd(3):length(cfgTr_ref.event)]; 
else
    evtToTake = 1:length(cfgTr_ref.event);
end

cfgTr_ref.trl = getCorrectTrl_2(cfgTr_ref, cfg1.trialdef.eventvalue, evtToTake);
cfgTr_ref.trl = [cfgTr_ref.trl zeros(size(cfgTr_ref.trl, 1), 1)]; 
trlInBetween = getCorrectTrl_2(cfgTr_ref, eventsInBetween, evtToTake); 
trlInBetween(trlInBetween == 0) = nan; 
for i = 1:size(trlInBetween, 1)
    trlInBetween(i, :) = trlInBetween(i, :) - cfgTr_ref.trl(i, 1) + 1; 
end

% Epoch data into the above sections
disp('Loading data')
cfg1 = [];
cfg1.dataset = dataset;
cfg1.channel = 1:128;
cfg1.demean = 'yes';
dat = ft_preprocessing(cfg1);

% Downsampling and De-meaning the data
cfg1 = [];
cfg1.demean = 'yes'; % Demean the data
cfg1.resamplefs = resample_fs;
dat = ft_resampledata(cfg1, dat);

orig_fs = dat.cfg.origfs;

% Filtering and reref
[dat] = filtering(dat, bp_freq);

% Downsampling 'event'
ds_factor = orig_fs/resample_fs;
trl = round(cfgTr_ref.trl/ds_factor); 
trlInBetween = round(trlInBetween/ds_factor); 

% If there were 2 events. Otherwise with prestim/poststim, the original would be good enough 
% MAKE CHANGES HERE TOMORROW -- SC 20201104
%{
if length(events) == 2
    
    indToTake = []; 
    
    % Taking 500 samples before first event and 500 samples after the 2nd
    % event
    if ((s == 5 && info.block == 2) || s == 1) && strcmp(info.task, 'afc')
        if s == 5 && info.block == 2
            indToTake = [95:294, 374:size(trl, 1)]; 
        elseif s == 1 && info.block == 1
            indToTake = [1:300, 380:579, 647:size(trl, 1)]; 
        elseif s == 1 && info.block == 2
            indToTake = [1:300, 393:size(trl, 1)]; 
        end
    else
        indToTake = 1:size(trl, 1); 
    end
    trl = trl(indToTake, :); 
    
    % Instead of doing this, any inconsistencies can 
    trl = [trl(1:2:end, 1) trl(2:2:end, 1) trl(1:2:end, 3)]; % Assuming 0 offset
end
%}

% Rounding the downsampled sample numbers
% trl = round(trl); 

end