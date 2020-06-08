close all;
clear all; %#ok<CLALL>

dbstop if error;

%% Give the subject and block details here
subjects = [6, 7, 8, 10, 11, 12, 13, 14, 16]; % Enter subjects to analyze
mff_keyword = 'afc';
dashes = '----------------';

%% Give preprocessing details here
acquisition_system = 'EGI'; % 'Biosemi'or 'EGI'
orig_fs            = 1000; % Original sampling frequency
resample_fs        = 250; % Downsample frequency
events_req         = {'DI60', 'D120'}; % Events to epoch
layout             = 'GSN-HydroCel-128.sfp'; % layout filename
bp_freq            = [0.1 45]; % band pass filtering frequency range.

%% Add all relevant toolboxes to path
ftPath = '/Volumes/Untitled/attend-expect/server_copy/fieldtrip';
ntPath = '/Volumes/Untitled/attend-expect/server_copy/toolboxes/NoiseTools20191024';
homeDir = '/Volumes/Untitled/attend-expect/EEG_Data'; % The directory in which raw EEGData should be put in 'Raw' and where the processed data will be stored as well

% Adding FieldTrip, NoiseTools, SCADS
pd = pwd; cd(ftPath); ft_defaults; cd(pd); % Fieldtrip toolbox
addpath(genpath(ntPath)); % NoiseTools toolbox
addpath('./SCADS'); % Adding SCADS code to path

data_dir = [homeDir '/Raw/'];
save_dir = [homeDir '/Processed/'];

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
        if ~exist(results_dir, 'dir')
            mkdir(results_dir);
        end
        files = dir([block_dir '*' mff_keyword '*.mff']);
        
        if ~exist([results_dir '/artfct_elecs_ICA.mat'])
            
            session_list = 1:length(files);
            blocks_list  = [];
            
            % ========  Creating blocks ==============
            for session = session_list
                
                dataset = [block_dir files(session).name]; % EEG file name
                
                % Information
                info.subject = s;
                info.session = session;
                info.results_dir = results_dir;
                
                % Create data epoched at pauses as blocks, filter and downsample data
                if ~exist([results_dir '/0_S' num2str(s, '%.2d') '_session_' num2str(session, '%.2d') '.mat'], 'file')
                    disp('Saving sessions as blocks by separating at pauses')
                    [dat, trl] = create_blocks_filter_resample(dataset, events_req, bp_freq, resample_fs, info);
                    save([results_dir '/0_S' num2str(s, '%.2d') '_session_' num2str(session, '%.2d')], 'dat', 'trl');
                end
                
                clearvars dat trl
                
            end
            
            % ========== Preprocessing ===============
            datAll = {}; evtAll = {};
            
            for session = session_list
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
            
            clearvars comp datAll datForICA datAllRej;
        end
    end
end

%% STEP 2: ICA rejection, Visual trial rejection
disp('==== Running STEP 2 -- ICA component rejection, interpolation, visual rejection ===='); 
for s = subjects
    for block = 1:numBlocks
        
        disp([dashes ' Subject ' num2str(s, '%.2d') ' ' dashes]);
        results_dir = [save_dir 'subject' num2str(s, '%.2d') '/' mff_keyword '/' 'block' num2str(block, '%.2d') '/'];
        block_dir = [subject_dir 'block' num2str(block, '%.2d'), '/'];
        files = dir([block_dir '*' mff_keyword '*.mff']);
        session_list = 1:length(files); 
        
        % Can make the following cleaner
        % For this we might have to save the appended data file 'datAll' in
        % the previous section. Let's do that -- TODO 
        % ========== Loading and appending data ===============
        datAll = {}; evtAll = {};
        
        for session = session_list
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
        
        % Load ICA components
        load([results_dir 'artfct_elecs_ICA']); 
        % ==========================================
        
        % ========= ICA component rejection, interpolation =============
        if ~exist([results_dir '/1_S' num2str(s, '%.2d') '.mat'])
            % Visualize and Remove ICA components here
            visualizeICAcomponents;
            
            % Recreate the data after rejecting components
            cfg = []; cfg.channel = datAll.label(setdiff(1:128, badElecs));
            datForICA = ft_preprocessing(cfg, datAll); % choosing good elecs
            datICARej = ft_rejectcomponent(cfgICA, comp, datForICA);
            
            % ===== Interpolation =====
            % For this, we need a 'neighbours' structure.
            cfg1 = []; cfg1.method = 'triangulation'; cfg1.elec = ft_read_sens('./GSN-HydroCel-128.sfp');
            cfg1.channel = datAll.label;
            neighbours = ft_prepare_neighbours(cfg1);
            cfg = []; cfg.missingchannel = datAll.label(badElecs); cfg.neighbours = neighbours;
            cfg.elec = ft_read_sens('./GSN-HydroCel-128.sfp'); cfg.method = 'spline';
            datInterp = ft_channelrepair(cfg, datICARej);
            
            % Manual reordering of channels
            datInterp = reorderChannelEGI(datInterp);
            
            % We will be saving some other variable once we do the visual
            % rejection (or maybe overwrite it)
            save([results_dir '/1_S' num2str(s, '%.2d')], 'datInterp', 'badElecs');
        else
            load([results_dir '/1_S' num2str(s, '%.2d')]);
        end
        % =========================================================
        
        % ================ VISUAL TRIAL REJECTION =======================
        if ~exist([results_dir '/2_S' num2str(s, '%.2d') '.mat'])
            % ========== detrend here [Should we do it before?] ========
            cfg = []; cfg.detrend = 'yes';
            datInterp = ft_preprocessing(cfg, datInterp);
            
            % ======= Flag trials here =========
            [flagTr, vsMat, hfMat] = artifact_rejection_jfos(datInterp.trial, 1000/datInterp.fsample);
            
            % Select bad trials using ft_preprocessing
            cfg =[]; cfg.trials = flagTr;
            badTrDat = ft_selectdata(cfg, datInterp);
            
            %
            cfg = []; cfg.method = 'trial';
            cfg.keeptrial = 'yes';
            veryBadTrDat = ft_rejectvisual(cfg, badTrDat);
            
            % Find very bad trials
            veryBadTrInd = find(ismember(veryBadTrDat.sampleinfo(:, 1), veryBadTrDat.cfg.artfctdef.trial.artifact(:, 1)));
            veryBadTr = flagTr(veryBadTrInd);
            
            % Make very bad trials nan
            for tr = 1:length(veryBadTr)
                datInterp.trial{veryBadTr(tr)} = nan(size(datInterp.trial{veryBadTr(tr)}));
            end
            
            % Make select electrodes from other trials nan
            badTr = setdiff(flagTr, veryBadTr);
            
            % Use vsMat and hfMat here to get the electrodes that need
            % rejection.
            artfctMat = vsMat + hfMat;
            artfctBin = artfctMat > 0;
            % Get electrodes which are bad for the badTr
            for tr = 1:length(badTr)
                datInterp.trial{badTr(tr)}(artfctBin(badTr(tr), :) > 0, :) = ...
                    nan(sum(artfctBin(badTr(tr), :) > 0), size(datInterp.trial{badTr(tr)}, 2));
            end
            
            save([results_dir '/2_S' num2str(s, '%.2d')], 'datInterp', 'flagTr',...
                'vsMat', 'hfMat', 'badTr', 'veryBadTr', 'badElecs');
        else
            disp([results_dir '/2_S' num2str(s, '%.2d') '   ALREADY EXISTS']);
        end
        % ===========================================================
        
        clear datInterp comp datAll datForICA; 
    end
end

% Cleaning up
rmpath(genpath(ntPath)); % NoiseTools toolbox
rmpath('./SCADS');