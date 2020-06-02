function epoch_sample_idx = epoch_data_1(evtNew, event_flags, dat)

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
