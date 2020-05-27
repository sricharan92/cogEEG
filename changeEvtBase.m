% Changes the events in '255' system to '127' system
%
% SC

function events = changeEvtBase(events, evtType)

for evtNum = 1:length(events)
    if strcmp(events(evtNum).type, evtType)
        if contains(events(evtNum).value, 'DI')
            evtVal = str2num(events(evtNum).value(3:end));
        else
            evtVal = str2num(events(evtNum).value(2:end));
        end
        if evtVal >= 128
            newEvtVal = evtVal - 128;
            if newEvtVal < 100
                startStr = 'DI'; 
            else
                startStr = 'D'; 
            end
            events(evtNum).value = [startStr num2str(newEvtVal)]; 
        end
    end
end

end