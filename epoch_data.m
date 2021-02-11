function [epochs] = epoch_data(data, events, event_flags, block, cell_flag)
%%%% cell_flag - 0 - output 3D matrix; 1 - output cell.

epoch_sample_idx = [];
for i = 1:length(events)
    if strcmp(event_flags{1}, events(i).value)
        epoch_sample_idx(end+1, 1) = events(i).sample; %#ok<AGROW>
    elseif strcmp(event_flags{2}, events(i).value)
        epoch_sample_idx(end, 2) = events(i).sample;
    end
end

no_of_trials = length(epoch_sample_idx);
disp(['Number of trials in Block ' num2str(block) ' is ' num2str(no_of_trials)])

epoch_lengths = epoch_sample_idx(:, 2) - epoch_sample_idx(:, 1) + 1;

if cell_flag == 0
    max_length = max(epoch_lengths);

    epochs = nan(max_length, size(data, 2), no_of_trials); % samples x electrodes x trials
    for i = 1:no_of_trials
        epochs(1:epoch_lengths(i), :, i) = data(epoch_sample_idx(i, 1):epoch_sample_idx(i, 2), :);
    end
elseif cell_flag == 1
    epochs = cell(no_of_trials, 1);
    for i = 1:no_of_trials
        epochs{i} = data(epoch_sample_idx(i, 1):epoch_sample_idx(i, 2), :); % samples x electrodes
    end
end

end