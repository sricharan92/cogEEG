close all;
clear all; %#ok<CLALL>

dbstop if error;

%% Give the subject and block details here
subjects = [11]; % Enter subjects to analyze
mff_keyword = 'afc';
dashes = '----------------';

%% Give preprocessing details here
acquisition_system = 'EGI'; % 'Biosemi'or 'EGI'
orig_fs            = 1000; % Original sampling frequency
resample_fs        = 250; % Downsample frequency
ICA_flag           = 0; % 0 - run ICA; 1 - skip ICA
events_req         = {'DI60', 'D120'}; % Events to epoch
layout             = 'GSN-HydroCel-128.sfp'; % layout filename
samp_omit_scads    = 16*resample_fs; % samples to omit before electrode rejection using SCADS to avoid filter artifacts' effects on rejection - (start and end of data)
bp_freq            = [0.1 80]; % band pass filtering frequency range.
alpha_flag         = 0; % alpha suppression flag: 0 - do not suppress; 1 - suppress.

%% Add all relevant toolboxes to path
ftPath = '/Volumes/Untitled/attend-expect/server_copy/fieldtrip';
ntPath = '/Volumes/Untitled/attend-expect/server_copy/toolboxes/NoiseTools20191024';
homeDir = '/Volumes/Untitled/attend-expect/EEG_Data'; % The directory in which raw EEGData should be put in 'Raw' and where the processed data will be stored as well

% Adding FieldTrip, NoiseTools, SCADS
pd = pwd; cd(ftPath); ft_defaults; cd(pd); % Fieldtrip toolbox
addpath(genpath(ntPath)); % NoiseTools toolbox
addpath('./SCADS'); % Adding SCADS code to path

data_dir = [homeDir '/Raw/'];
save_dir = [homeDir '/Processed_new1/'];

%% EEG electrodes' layout
dist = electrode_dist_and_polar_ang(acquisition_system);

%% Loop for all subjects
for s = subjects
    subject_dir = [data_dir 'subject' num2str(s, '%.2d') '/' mff_keyword '/'];
    allCont = dir(subject_dir);
    folders = allCont([allCont(:).isdir]);
    folders = folders(~ismember({folders(:).name}, {'.', '..'})); % Ignore '.' and '..'
    numBlocks = length(folders);
    
    for block = 1:numBlocks
        
        disp([dashes ' Subject ' num2str(s, '%.2d') ' ' dashes])
        block_dir = [subject_dir 'block' num2str(block, '%.2d'), '/'];
        results_dir = [save_dir 'subject' num2str(s, '%.2d') '/' mff_keyword '/' 'block' num2str(block, '%.2d') '/'];
        
        files = dir([block_dir '*' mff_keyword '*.mff']);
        if ~exist(results_dir, 'dir')
            mkdir(results_dir);
        end
        
        session_list = 1:length(files);
        blocks_list  = [];
        
        %% Creating blocks
        for session = session_list
            dataset = [block_dir files(session).name]; % EEG file name
            
            % Information
            info.subject = s;
            info.session = session;
            info.results_dir = results_dir;
            
            % Create data epoched at pauses as blocks, filter and downsample data
            if ~exist([results_dir '/0_S' num2str(s, '%.2d') '_session_' num2str(session, '%.2d') '.mat'], 'file')
                disp('Saving sessions as blocks by separating at pauses')
                [dat, trl] = create_blocks_filter_resample(dataset, events_req, bp_freq, alpha_flag, resample_fs, info);
                save([results_dir '/0_S' num2str(s, '%.2d') '_session_' num2str(session, '%.2d')], 'dat', 'trl');
            end
            
            clearvars dat trl
            
        end
        
        %% Preprocessing
        datAll = {}; evtAll = {};
        
        for session = session_list
            
            datEpoched = {}; count = 1;
            dataset = [block_dir files(session).name]; % EEG file name
            
            % Load data epoched at pauses as blocks
            disp('Loading saved events and downsampled data')
            load([results_dir '/0_S' num2str(s, '%.2d') '_session_' num2str(session, '%.2d')])
            cfg = []; cfg.trl = trl;
            
            datAll{session} = ft_redefinetrial(cfg, dat);
            
        end
        
        % Append the data if the folder had many .mff files
        if length(datAll) > 1
            cfg = [];
            cfg.keepsampleinfo = 'no';
            datAll = ft_appenddata(cfg, datAll{:});
        else
            datAll = datAll{1}; 
        end
        
        if ~exist([results_dir '/artfct_elecs_ICA.mat'])
            % ------ nt_find_bad_channels -----
            concatDat = cat(2, datAll.trial{:});
            concatDat = concatDat'; % nt needs time x channels
            thresh1 = 3; thresh2 = 100; thresh3 = 3; proportion = 0.33;
            [badElecs, ~] = nt_find_bad_channels(concatDat, proportion, thresh1, thresh2, thresh3);
            clear concatDat;
            
            % ------ ft_artifact_zvalue -------
            cfg = [];
            cfg.continuous = 'no';
            cfg.artfctdef.jump.channel = setdiff(1:128, badElecs);
            cfg.artfctdef.jump.cutoff =20;
            cfg.artfctdef.jump.interactive = 'no';
            cfg.artfctdef.zvalue.cumulative = 'yes';
            [cfg_artifact_jump, artifact_jump] = ft_artifact_jump(cfg, datAll);
            
            % Reject artifacts -- At present removes the trials with artifacts completely
            cfg_artifact_jump.artfctdef.reject = 'partial'; 
            datAllRej = ft_rejectartifact(cfg_artifact_jump, datAll);
            
            % Run ICA
            datForICA = cat(2, datAllRej.trial{:});
            cfg = [];
            cfg.channel      = setdiff(1:128, badElecs);
            n_ic = rank(datForICA(cfg.channel, :)); % Rank of the matrix without bad elecs
            cfg.method       = 'runica';
            cfg.numcomponent = n_ic;
            comp = ft_componentanalysis(cfg, datAllRej);
            
            save([results_dir '/artfct_elecs_ICA'], 'comp', 'cfg_artifact_jump', 'badElecs');
        else
            load([results_dir '/artfct_elecs_ICA']);
        end
        
        if ~exist([results_dir '/1_S' num2str(s, '%.2d') '.mat'])
            % Visualize and Remove ICA components here
            visualizeICAcomponents;
            
            % Recreate the data after rejecting components
            % ISSUE: Artifact rejection -- The 'partial' rejection increases
            % number of trials and the 'complete' rejection decreases it. Thus,
            % we're unable to recreate the data using JUST the components. i.e.
            % I have passed 'datAll' also here. But then the backprojection has
            % 128 electrodes when we actually gave ICA lesser electrodes.
            % It's probably okay to interpolate those electrodes but have to
            % think about it.
            datICARej = ft_rejectcomponent(cfgICA, comp, datAll);
            
            % ===== Interpolation =====
            % For this, we need a 'neighbours' structure.
            cfg1 = []; cfg1.method = 'triangulation'; cfg1.elec = ft_read_sens('./GSN-HydroCel-128.sfp');
            cfg1.channel = datAll.label;
            neighbours = ft_prepare_neighbours(cfg1);
            cfg = []; cfg.badchannel = datAll.label(badElecs); cfg.neighbours = neighbours;
            cfg.elec = ft_read_sens('./GSN-HydroCel-128.sfp'); cfg.method = 'spline';
            datInterp = ft_channelrepair(cfg, datICARej);
            
            % We will be saving some other variable once we do the visual
            % rejection (or maybe overwrite it)
            save([results_dir '/1_S' num2str(s, '%.2d')], 'datInterp');
        else
            load([results_dir '/1_S' num2str(s, '%.2d')]);
        end
        
        % ======= Flag trials here =========
        flagTr = artifact_rejection_jfos(datInterp.trial, 1000/datInterp.fsample);
        
        % ======= Visual rejection of trials here ========
        cfg = []; cfg.method = 'trial'; cfg.neighbours = datInterp.cfg.neighbours; 
        cfg.trials = setdiff(1:length(datInterp.trial), flagTr); cfg.keeptrial = 'yes';
        ft_rejectvisual(cfg, datInterp); 
        
        %disp('Clearing variables');
        %clearvars -except subjects mff_keyword dashes results_folder acquisition_system orig_fs ...
        %    resample_fs ICA_flag events_req layout samp_omit_scads home_dir data_dir save_dir dist ...
        %    ntPath trial_end bp_freq alpha_flag;
        
    end
end

% Cleaning up
rmpath(genpath(ntPath)); % NoiseTools toolbox
rmpath('./SCADS');