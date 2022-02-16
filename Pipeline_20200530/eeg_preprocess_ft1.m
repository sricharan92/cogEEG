function eeg_preprocess_ft1(info) 

trialNo = 1;
for block = info.blocksList
    if info.ICA_flag
        load(sprintf('%s/2_S%d_B%d', info.resultsFolder, info.subjectNumber, block));
    else
        load(sprintf('%s/1_S%d_B%d', info.resultsFolder, info.subjectNumber, block),'dataSmallEpoch');
    end
    %% Combine the epoched data (combine the data obtained after epoching of 200s)
    x = []; w = []; trl = [];
    if size(dataSmallEpoch.trial,2) ==1
        x = dataSmallEpoch.trial{1,1}';
        w = dataSmallEpoch.w_inpaint{1,1};
        trl = dataSmallEpoch.trlSamples{1,1}; 
    elseif size(dataSmallEpoch.trial,2) >= 2
        for epochNo = 1: size(dataSmallEpoch.trial,2)-1
            x = [x; dataSmallEpoch.trial{1,epochNo}'];
            w = [w; dataSmallEpoch.w_inpaint{epochNo,1}];
            trl = [trl; dataSmallEpoch.trlSamples{1,epochNo}];
        end
        ll = max(dataSmallEpoch.samples{end,1}) - size(x,1);
        temp = dataSmallEpoch.trial{1,epochNo + 1}';
        x = [x; temp(end - ll + 1 : end, :)];
        w = [w; dataSmallEpoch.w_inpaint{epochNo + 1,1}(end - ll + 1 : end, :)];
        ll_trl = find(trl(end,1)==dataSmallEpoch.trlSamples{1,epochNo+1}) + 1;
        trl = [trl; dataSmallEpoch.trlSamples{1,epochNo+1}(ll_trl:end,:)];
    end
    
    %% Robust re-referencing
    x=nt_rereference(x,w);
    
    if info.fig_flag
        figNo = 1; figure(figNo); clf
        plot((1:size(x,1))/info.resamplefs, x); title('re-referenced');
        xlabel('Time (s)'); ylabel('\muV')
    end
    
    x = x - median(x); %Demean x
        %% Epoch the data according to trial start and trial end
    cfg1 = [];
    cfg1.dataset = sprintf('%s/S%d_B%d.bdf',info.dataFolder,info.subjectNumber,block); % EEG file name
    cfg1.origfs = info.origfs;
    cfg1.resamplefs = info.resamplefs;
    cfg1.acquisitionSystem = info.acquisitionSystem;
    eventStart = 10; eventEnd = 100;
    cfgTr = defTr(eventStart, eventEnd, cfg1);

    %% Epoch the data according to trial start and trial end
    nchans = size(x,2);
    ntrials = size(trl,1);
    nsamples = max(trl(:,2) - trl(:,1)) + 1;

    %Do some cleaning on the epoched data
    if ~exist('xx_detrended')
        xx_detrended.data = nan(3*nsamples,nchans, ntrials*max(info.blocksList)*2); % Assign a nan matrix of a large size. Trim it later
        if info.events_flag
            xx_detrended.events = nan(3*nsamples, ntrials*max(info.blocksList)*2);
        end
    end
    for tr = 1:ntrials
        % temp = nt_detrend(x(trl(tr,1) : trl(tr,2),:), 1); % Fist order detrend
        temp = x(trl(tr,1) : trl(tr,2),:); % Fist order detrend
        xx_detrended.data(1:trl(tr,2) - trl(tr,1) + 1,:,trialNo) = nt_demean(temp);
        xx_detrended.samples(trialNo,:) = [trl(tr,1), trl(tr,2)];
        xx_detrended.blockNumber(trialNo,1) = block; 
        
        if info.events_flag
            fromIdx = find(cfgTr.event_ds.sample == cfgTr.trl(tr,1));
            toIdx = find(cfgTr.event_ds.sample == cfgTr.trl(tr,2));
            ev = []; 
            for id = fromIdx : (toIdx -1)
                ev = [ev; cfgTr.event_ds.value{id}*ones(cfgTr.event_ds.sample(id + 1) - cfgTr.event_ds.sample(id),1)];
            end
            ev = [ev; cfgTr.event_ds.value{id}];
            if size(ev,2) ~=0 
            xx_detrended.events(1:cfgTr.trl(tr,2) - cfgTr.trl(tr,1) + 1,trialNo) = ev; 
            else
                xx_detrended.events(1:cfgTr.trl(tr,2) - cfgTr.trl(tr,1) + 1,trialNo)  = nan; 
            end
        end 
        trialNo = trialNo + 1; 
    end
end

% Trim the xx_detrended structure to save memory
aa = isnan(xx_detrended.data); 
idxTrials = squeeze(mean(mean(aa,1),2)) == 1; % one implies all nan and hence non-existent trials
idxTimePoints = squeeze(mean(mean(aa,2),3)) == 1; % one implies all nan and hence non-existent time points
xx_detrended.data(:,:,idxTrials) = []; 
xx_detrended.data(idxTimePoints,:,:) = []; 
if info.events_flag
    xx_detrended.events(:,idxTrials) = [];
    xx_detrended.events(idxTimePoints,:) = [];
end

%% Remove noisy trials - SCADS

xx_scads = xx_detrended; 
xx_toDSS = xx_detrended.data;  

addpath('SC_Scads')
lambda = 2; % As per the SCADS paper
% load the polar angle variable depending on the acquisition system
if strcmp(info.acquisitionSystem, 'Biosemi')
    load('elecLocsBiosemi.mat')
    polarAng = elecLocs.sph_phi; 
elseif strcmp(info.acquisitionSystem, 'EGI')
    load('elecLocsEGI.mat')
    polarAng = elecLocs.sph_phi; 
else
    disp('elecLocs file not available for acquisition systems other than Biosemi and EGI.')
end

[maxAmp, stDev, maxGrad] = cog_scads_1_2(xx_detrended.data); 
rejMaxAmpEp = cog_scads_1_4(maxAmp, lambda, polarAng); 
rejStDevEp = cog_scads_1_4(stDev, lambda, polarAng); 
rejMaxGradEp = cog_scads_1_4(maxGrad, lambda, polarAng); 
trialRej = [rejMaxAmpEp; rejStDevEp; rejMaxGradEp];

rejectedTrials = []; 
for i = 1:size(trialRej, 1)
    for elec = 1:size(trialRej, 2)
        for tr = trialRej{i, elec}
            xx_scads.data(:,elec,tr) = NaN;
            rejectedTrials = [rejectedTrials; tr];
        end
    end
end
rejectedTrials = unique(rejectedTrials);
xx_detrended.data(:,:,rejectedTrials) = []; % Remove rejected Trials
rmpath('SC_Scads')

%% DSS to emphasize repeatablity (Optional Step)
if info.DSS_flag
    temp = mean(xx_toDSS(:,1,:),3);
    xx_toDSS = xx_toDSS(~isnan(temp),:,:); % Use the minimum epoch size to apply DSS
    
    [todss,pwr0,pwr1]=nt_dss1(xx_toDSS);
    fromdss=pinv(todss);
    
    if info.fig_flag
        z=nt_mmat(xx_toDSS,todss);
        
        figure(8); clf;
        subplot(1,2,1);
        plot(pwr1./pwr0,'.-'); xlabel('component'); ylabel('score'); title ('repeatability DSS');
        subplot(1,2,2);
        nt_bsplot(z(:,1,:)); title('best DSS component');
        
        % denoise by selecting best components and projecting back to sensor space
        NKEEP=7;
        figNo = figNo+1; figure(figNo); clf;
        subplot 121; plot(mean(xx_toDSS,3)); title('before DSS denoising')
        xx=nt_mmat(xx_toDSS,todss(:,1:NKEEP)*fromdss(1:NKEEP,:));
        subplot 122; plot(mean(xx,3)); title('after');
    end
end

%% Save the data
% Save the data after loading, filtering and downsampling
if info.DSS_flag
    save(sprintf('%s/3_S%d_filt', info.resultsFolder, info.subjectNumber), 'xx_scads', 'todss', 'fromdss'); 
else
    save(sprintf('%s/3_S%d', info.resultsFolder, info.subjectNumber), 'xx_scads');
end

if info.fig_flag
    % Save all images
    FolderName = info.resultsFolder;   % Your destination folder
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        FigName   = num2str(get(FigHandle, 'Number'));
        set(0, 'CurrentFigure', FigHandle);
        saveas(FigHandle, fullfile(FolderName, ['PreprocessFig', FigName '.png']));
    end
end
close all; % Close all figures 

end