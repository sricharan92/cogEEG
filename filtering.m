function [data] = filtering(data, bp_freq, alpha_flag, dataset)

disp('Filtering data')
cfg = [];

if nargin < 3
    alpha_flag = 0;
end

if nargin > 3
    cfg.dataset = dataset;
end

cfg.channel = 1:128;
cfg.demean = 'yes';

%% Filter
cfg.bpfilter   = 'yes';
cfg.bpfreq     = bp_freq; % [0.1, 35];
cfg.bpfiltord  = 4;
cfg.bpfilttype = 'but';

%% Line noise removal
cfg.bsfilter    = 'yes';
cfg.bsfreq      = [49 51];
cfg.bsfiltord   = 4;
cfg.bsfilttype  = 'but';

%% Alpha filtering
if alpha_flag == 1
    cfg.bsfreq = [cfg.bsfreq; 8 12];
end

if nargin > 3
    data = ft_preprocessing(cfg);
else
    data = ft_preprocessing(cfg, data);
end

end