function [dat, trl] = create_blocks_filter_resample(dataset, events, bp_freq, resample_fs, info)

% Information
s           = info.subject;
session     = info.session;
results_dir = info.results_dir;

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
    cfgTr_ref = ft_definetrial(cfg1); 
    % ------------------
    
    % Shows you all the unique DIN events and their counts
    if dispEvts == 1
        listEventsEGI(cfgTr_ref.event, cfg1.trialdef.eventtype); 
    end
    
    % Converts all events to the 127 system 
    if flag128 == 1
        % Change all the events to 128 base
        cfgTr_ref.event = changeEvtBase(cfgTr_ref.event, cfg1.trialdef.eventtype); 
    end
    
    save([results_dir file_name], 'cfgTr_ref');
else
    disp([file_name ' already exists. Loading.'])
    load([results_dir file_name], 'cfgTr_ref')
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
[dat] = filtering(dat, bp_freq, alpha_flag);

% Downsampling 'event'
ds_factor = orig_fs/resample_fs;
trl = cfgTr_ref.trl/ds_factor; 

% If there were 2 events. Otherwise with prestim/poststim, the original would be good enough 
if length(events) == 2
    trl = [trl(1:2:end, 1) trl(2:2:end, 1) trl(1:2:end, 3)]; % Assuming 0 offset
end
    
% Rounding the downsampled sample numbers
trl = round(trl); 

end