% A script to test EGI events if they are all in 1 system or not 
% Right now only handles '_DINs' or '255_DINs'
%
% SC

function [] = listEventsEGIEpoched(events, evtType)

evtArr = {}; epoch = 0; evtInd = {}; 
evtStr = {'DI60', 'DI80', 'D120', 'DI12', 'DI32'}; 

% Get unique events
for evtNum = 1:length(events)
    if strcmp(events(evtNum).type, evtType)
        if contains(events(evtNum).value, 'DI')
            evtArr{epoch} = [evtArr{epoch}; str2num(events(evtNum).value(3:end))]; 
        else
            evtArr{epoch} = [evtArr{epoch}; str2num(events(evtNum).value(2:end))]; 
        end
        evtInd{epoch} = [evtInd{epoch}; evtNum];
    elseif strcmp(events(evtNum).type, 'epoch')
        disp('Incrementing Epoch'); 
        epoch = epoch + 1; 
        epochSample(epoch) = events(evtNum).sample; 
        evtArr{epoch} = []; 
        evtInd{epoch} = [];
    end
end

for ep = 1:length(evtArr)
    fprintf('EPOCH - %d \n', ep); 
    uniqEvt = unique(evtArr{ep});
    evtLt = [];
    for evt = 1:length(uniqEvt)
        evtLt(evt) = length(find(evtArr{ep} == uniqEvt(evt)));
        fprintf('Event Number %d: %d \n',  uniqEvt(evt), evtLt(evt));
    end
end

end