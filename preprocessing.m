close all;
clear all; %#ok<CLALL>
% clc;

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
%trial_end          = 'DI72'; % Corresponds to last flag in a trial.
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
        datAll = {}; evtAll = {}; comp = {};
        
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
        
        % ------ nt_find_bad_channels -----
        
        % ------ ft_artifact_zvalue -------
        
        % Run ICA
        datForICA = cat(2, datAll.trial{:});
        cfg = [];
        n_ic = rank(datForICA);
        %cfg.channel      = info.PreprocessingCommets.Sensors.AllGoodSensors;
        cfg.method       = 'runica';
        cfg.numcomponent = n_ic;
        comp{session} = ft_componentanalysis(cfg, datAll);
        
        % Need to save ICA here 
        
        %{
        for block = blockNoRef + 1:blockNoRef + length(dat.trial)
            
            fprintf(['\n' dashes ' Subject - %d, Block - %d ' dashes '\n\n'], s, block);
            sectName    = strcat('section', num2str(block - blockNoRef));
            cfgTr.event = event.(sectName);
            cfgTr.trl   = trl.(sectName);
            
            data       = dat;
            data.trial = data.trial(block - blockNoRef);
            data.time  = data.time(block - blockNoRef);
            
            % electrode rejection (whole block and trial-by-trial) [SCADS - visual - epoch - SCADS]
            %             [data_out, rej_elecs{block}, tr_rej_elecs{block}] = electrode_rejection(data.trial{1}', cfgTr.event, events_req, polar_ang, block, samp_omit_scads);
            if ~exist([results_dir '/1_S' num2str(s, '%.2d') '_block_' num2str(block, '%.2d') '.mat'], 'file')
                disp(['Saving preprocessed block ' num2str(block, '%.2d')])
                [data_out, rej_elecs, tr_rej_elecs] = electrode_rejection(data.trial{1}', cfgTr.event, events_req, polar_ang, block, samp_omit_scads);
                save([results_dir '/1_S' num2str(s, '%.2d') '_block_' num2str(block, '%.2d')], 'data_out', 'rej_elecs', 'tr_rej_elecs');
            end
            
        end
        %}
        
        disp('Clearing variables');
        clearvars -except subjects mff_keyword dashes results_folder acquisition_system orig_fs ...
            resample_fs ICA_flag events_req layout samp_omit_scads home_dir data_dir save_dir dist ...
            ntPath trial_end bp_freq alpha_flag;
    end
end

% Cleaning up
rmpath(genpath(ntPath)); % NoiseTools toolbox
rmpath('./SCADS');