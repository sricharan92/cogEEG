
figure;
compNo = 1;
dat = cat(2,comp.trial{:});
timeRel = cat(2,comp.time{:});
%trEndIdx = find(diff(timeRel)<=0) + 1;
%trEndIdx = [trEndIdx, size(dat,2)];
%trEndIdxTime = trEndIdx/comp.fsample;
%trCueIdx = find(timeRel ==0);
%trCueIdxTime = trCueIdx/comp.fsample;
timeAbs = [1:size(dat,2)]/comp.fsample;

while compNo <= size(comp.unmixing, 1)
    
    clf; disp(compNo);
    componentTc = dat(compNo,:);
    
    d = movstd(componentTc, 5*comp.fsample); % moving standard deviation
    
    % Plotting STD, IC and spectrum
    subplot(2, 3, 4);
    histogram(d); xlabel('Standard Deviation'); ylabel('#');
    
    subplot(2, 3, 5);
    cfg1 = [];
    cfg1.component = compNo;
    cfg1.colorbar = 'yes';
    cfg1.layout = './GSN-HydroCel-128.sfp';
    ft_topoplotIC(cfg1, comp)
    
    subplot(2, 3, 6);
    f = comp.fsample*(0:length(componentTc)-1)/length(componentTc);
    plot(f,abs(fft(componentTc)));
    xlim([0.5 45]);
    xlabel('Frequency (Hz)');
    title('Frequency Spectrum');
    
    subplot(2, 3, [1, 2 3]);
    
    % Plotting the time course
    ymax = max(componentTc); ymin = min(componentTc);
    
    %xVals = [trCueIdxTime; trEndIdxTime];
    %yVals = repmat([ymin ymax]',1, size(xVals, 2));
    %line(xVals,yVals,'Color','blue' );
    
    %xVals = [repmat(trCueIdxTime, 2, 1); repmat(trEndIdxTime, 2, 1)];
    %yVals = repmat([ymin ymax ymax ymin]', 1, size(xVals, 2));
    %patch(xVals, yVals, [1 1 1]*0.8);
    
    hold on; 
    plot(timeAbs, componentTc);
    xlabel('Time (s)'); ylabel('micro volts');
    title(sprintf('Component plotted - %d, Subject - %d, blockNo - %d', compNo, s, block));
    
    pause;
    
    compNo = compNo + 1;
    
end

% Give components to reject
disp('---------------- PREPROCESSING - REJECTING ICA COMPONENTS ----------------');
c = input('Enter Components To Remove: ', 's');
cfgICA.component = str2num(c);
