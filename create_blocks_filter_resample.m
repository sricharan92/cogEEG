function [dat, event, trl] = create_blocks_filter_resample(dataset, events, trial_end, bp_freq, alpha_flag, resample_fs, info)

% Information
s           = info.subject;
session     = info.session;
results_dir = info.results_dir;

% Find pause indices
disp('Loading events')
cfg1 = [];
cfg1.dataset = dataset;
cfg1.trialdef.eventtype = '255_DINs';

cfg1.trialdef.eventvalue = {events{1}, events{2}};

% Check if events file already exists
file_name = ['events' num2str(s) '_session_' num2str(session)];
if ~exist([results_dir file_name '.mat'], 'file')
    disp(['Saving events as ' file_name '.'])
    cfgTr_ref = ft_definetrial(cfg1);
    save([results_dir file_name], 'cfgTr_ref');
else
    disp([file_name ' already exists. Loading.'])
    load([results_dir file_name], 'cfgTr_ref')
end

iEvents = cfgTr_ref.event;
% iEvents = ft_read_event(cfg1.dataset);

disp('Sorting events based on pauses')
pause_idx = [];
for i = 1:length(iEvents)
    if isempty(iEvents(i).value)
        pause_idx = [pause_idx i]; %#ok<AGROW>
    end
end

% Separate event flags of the above sections (each segment before a pause)
for i = 1:length(pause_idx)
    sectName = strcat('section', num2str(i));
    if i ~= length(pause_idx)
        event.(sectName) = iEvents(pause_idx(i):pause_idx(i+1)-1);
    else
        event.(sectName) = iEvents(pause_idx(i):end);
    end
end

% Correct for sample number
sample_offset = zeros(1, length(fieldnames(event)));
for i = 1:length(fieldnames(event))
    sectName = strcat('section', num2str(i));
    if i ~= 1
        sample_offset(i) = event.(sectName)(1).sample - 1;
        for k = 1:length(event.(sectName))
            event.(sectName)(k).sample = event.(sectName)(k).sample - sample_offset(i);
        end
    end
end

% Number of trials in each section
trl_idx = zeros(1, length(fieldnames(event)));
for j = 1:length(fieldnames(event))
    sectName = strcat('section', num2str(j));
    for i = 1:length(event.(sectName))
        if strcmp(event.(sectName)(i).value, trial_end)
            trl_idx(j) = trl_idx(j) + 1;
        end
    end
end

% Separate trial samples based on pauses
for i = 1:length(fieldnames(event))
    sectName = strcat('section', num2str(i));
    if i ~= 1
        idx1 = sum(trl_idx(1:i-1))*length(events) + 1;
    else
        idx1 = 1;
    end
    
    if i ~= length(fieldnames(event))
        idx2 = sum(trl_idx(1:i))*length(events);
    else
        idx2 = length(cfgTr_ref.trl);
    end
    
    trl.(sectName) = cfgTr_ref.trl(idx1:idx2, :);
    if i ~= 1
        trl.(sectName)(:, 1:2) = trl.(sectName)(:, 1:2) - sample_offset(i);
    end
end

% Epoch data into the above sections
disp('Loading data')
cfg1 = [];
cfg1.dataset = dataset;
cfg1.channel = 1:128;
cfg1.demean = 'yes';
dat = ft_preprocessing(cfg1);

disp('Epoching data at pauses')
for i = 1:length(pause_idx)
    time_offset = dat.time{1}(iEvents(pause_idx(i)).sample);
    if i ~= length(pause_idx)
        trial{i} = dat.trial{1}(:, iEvents(pause_idx(i)).sample:iEvents(pause_idx(i+1)).sample-1);
        time{i}  = dat.time{1}(:, iEvents(pause_idx(i)).sample:iEvents(pause_idx(i+1)).sample-1) - time_offset;
    else
        trial{i} = dat.trial{1}(:, iEvents(pause_idx(i)).sample:end);
        time{i}  = dat.time{1}(:, iEvents(pause_idx(i)).sample:end) - time_offset;
    end
end

dat.trial = trial;
dat.time  = time;

% Downsampling and De-meaning the data
cfg1 = [];
cfg1.demean = 'yes'; % Demean the data
cfg1.resamplefs = resample_fs;
dat = ft_resampledata(cfg1, dat);

orig_fs = dat.cfg.origfs;

% Filtering
[dat] = filtering(dat, bp_freq, alpha_flag);

% Downsampling 'event'
ds_factor = orig_fs/resample_fs;
sections = fieldnames(event);
for i = 1:length(sections)
    for k = 1:length(event.(sections{i}))
        event.(sections{i})(k).sample = floor((event.(sections{i})(k).sample + (ds_factor - 1))/4);
    end
    trl.(sections{i}) = floor((trl.(sections{i}) + (ds_factor - 1))/4);
end

% Removing few samples from the start and end of each block to avoid artifacts (maybe atrifacts appear because epochs end right before the response window)
del = 5 + 1; % window (number of samples) to be truncated at the beginning and end of the block with respect to first and last relevant events
for i = 1:length(dat.trial)
    dat.trial{i}(:, trl.(sections{i})(end, 1) + (del):end) = []; % truncate end of the block first
    dat.trial{i}(:, 1:trl.(sections{i})(1, 1) - (del))     = []; % truncate beginning of the block
    dat.time{i}(:, trl.(sections{i})(end, 1) + (del):end)  = []; % truncate end of the block first
    dat.time{i}(:, 1:trl.(sections{i})(1, 1) - (del))      = []; % truncate beginning of the block
    for k = 1:length(event.(sections{i}))
        event.(sections{i})(k).sample = event.(sections{i})(k).sample - trl.(sections{i})(1, 1) + del;
    end
end

end