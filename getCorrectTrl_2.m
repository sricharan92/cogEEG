%{

    ft_definetrial gives the trl according to events in
    'trialdef.eventtype'. After this, we convert the events to base 128.
    We do this to homogenize since in our data collection, the triggers
    switched from D255 to D127 randomly for some subjects. 

    The first cfgTr_ref.trl that we get from ft_definetrials will miss some
    events since the event numbers specified are in either D255 or D127

    This code gets the trl _after_ the homogenization. 

    SC - November 4th, 2020

%}

function trl = getCorrectTrl_2(cfgTr_ref, relEvents, evtToTake)

count = 1; trl = []; prevEvtInd = 0; 

for i = evtToTake
    if ~isempty(cfgTr_ref.event(i).value) && ~strcmp(cfgTr_ref.event(i).type, 'BAD')
        if sum(ismember(relEvents, cfgTr_ref.event(i).value))
            evtInd = find(ismember(relEvents, cfgTr_ref.event(i).value));
            
            trl(count, evtInd) = cfgTr_ref.event(i).sample;
            
            %{
            if evtInd == 1 && prevEvtInd > 0
                count = count+1; % Next trial... Assumes fixation (or first event) not lost for all trials
            end
            %}
            
            if evtInd == length(relEvents)
                count = count+1; % Next trial... Assumes fixation (or first event) not lost for all trials
            end
            
            prevEvtInd = evtInd;
        end
    end
    %{
    if strcmp(cfgTr_ref.event(i).value, relEvents{1})
        trl(count, 1) = cfgTr_ref.event(i).sample;
        count = count + 1;
    elseif strcmp(cfgTr_ref.event(i).value, relEvents{2})
        trl(count, 1) = cfgTr_ref.event(i).sample;
        count = count + 1;
    end
    %}
end

%trl = repmat(trl, [1, 2]); 
%trl(:, 3) = zeros(size(trl, 1), 1); 

end