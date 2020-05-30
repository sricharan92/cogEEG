function eeg_rejectICAComp(info_ica)


%% Add functions to path
addpath('/home/cog-lab/Documents/MATLAB/Toolboxes/eeglab2019_1'); % EEG lab
eeglab; 

load(sprintf('%s/Sub%d/1_Filt.mat', info_ica.savePath, info_ica.subjectNumber), 'eegFilt', 'info');

for idx = 1:length(eegFilt)
    
    % Convert FT data structure to eeglab data structure
    data_ft = eegFilt{idx}; 
    
    % Apply ICA
    cfg1 = [];
cfg1.method = 'runica'; % this is the default and uses the implementation from EEGLAB
comp = ft_componentanalysis(cfg1, data_ft); 

% Run SASICA toolbox
eeg.srate = comp.fsample; 
eeg.icawinv = comp.topo;
eeg.chanlocs = table2struct(elecLocs);
eeg.data = comp.trial{1,1}; 

cfg1 = [];
cfg1.autocorr.enable = 1;
cfg1.focalcomp.enable = 1;
cfg1.chancorr.enable = 1;
cfg1.ADJUST.enable = 1;

out = eeg_SASICA(eeg, cfg1);
    
    
    hdr.Fs = data_ft.fsample;
    hdr.nChans = size(data_ft.label,1); 
    hdr.label = data_ft.label;
    hdr.nTrials = size(data_ft.trial,1);
    hdr.nSamples = size(data_ft.trial{1,1},2);
    hdr.nSamplesPre = 0;
    hdr.chantype = cell(128,1);
    hdr.chantype(:) = {'eeg'};
    hdr.chanunit = cell(128,1);
    hdr.chanunit(:) = {'uV'};
    
    data_eeg = fieldtrip2eeglab(hdr, data_ft.trial{1}); 
    ica = pop_runica(data_eeg, 'icatype', 'runica', 'reorder', 'on','interrupt','off');
    ica.chanlocs = table2struct(elecLocs); % table2struct(BiosemiEEGLocs);
%     [ica.chanlocs(:).radius] = deal(sqrt(BiosemiEEGLocs.X(1).^2 + BiosemiEEGLocs.Y(1).^2 + BiosemiEEGLocs.Z(1).^2));
    
    %% Visualize components and remove the noisy ones
    out = SASICA(ica);
    
end

end