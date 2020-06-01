function epoch_sample_idx = epoch_data_1(event, event_flags, dat)

evtNew = {};
evtNew{1} = event.section1;
evtNew{2} = event.section2;
evtNew{3} = event.section3;
evtNew{4} = event.section4;
evtNew{5} = event.section5;

epoch_sample_idx = {};
for trial = 1:length(dat.trial)
    events = evtNew{trial};
    count = 1; 
    for i = 1:length(events)
        if strcmp(event_flags{1}, events(i).value)
            epoch_sample_idx{trial}(count, 1) = events(i).sample; %#ok<AGROW>
        elseif strcmp(event_flags{2}, events(i).value)
            epoch_sample_idx{trial}(count, 2) = events(i).sample;
            epoch_sample_idx{trial}(count, 3) = 0; 
            count = count + 1; 
        end
    end
end

end