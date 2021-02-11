% A script to test EGI events if they are all in 1 system or not 
% Right now only handles '_DINs' or '255_DINs'
%
% SC

function [] = listEventsEGI(events, evtType)

evtArr = [];
% Get unique events
for evtNum = 1:length(events)
    if strcmp(events(evtNum).type, evtType)
        if contains(events(evtNum).value, 'DI')
            evtArr = [evtArr; str2num(events(evtNum).value(3:end))]; 
        else
            evtArr = [evtArr; str2num(events(evtNum).value(2:end))]; 
        end
    end
end

uniqEvt = unique(evtArr); 
evtLt = []; 
for evt = 1:length(uniqEvt)
    evtLt(evt) = length(find(evtArr == uniqEvt(evt))); 
    fprintf('Event Number %d: %d \n',  uniqEvt(evt), evtLt(evt)); 
end 

end