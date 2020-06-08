function [data] = filtering(data, bp_freq)

disp('Filtering data')
cfg = [];

cfg.channel = 1:128;
cfg.demean = 'yes';

cfg.reref = 'yes'; 
cfg.refchannel = {'E57', 'E100'}; % Check if this works 

%% Does padding work for data already epoched? In the help page it says that feature used only if data loaded from file

%% Filter
cfg.bpfilter   = 'yes';
cfg.bpfreq     = bp_freq; % [0.1, 35];
cfg.bpfiltord  = 3;
cfg.bpfilttype = 'but';

%% Line noise removal
cfg.bsfilter    = 'yes';
cfg.bsfreq      = [49 51];
cfg.bsfiltord   = 3;
cfg.bsfilttype  = 'but';

% Removing line noise using dftfilter
%cfg.dftfilter = 'yes'; cfg.dftfreq = [50, 100, 150]; cfg.dftreplace = 'neighbour'; 

data = ft_preprocessing(cfg, data);


end