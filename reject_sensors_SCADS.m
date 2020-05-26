function [data_out, rejected_sensors] = reject_sensors_SCADS(data, polar_ang, samp_omit_scads)
%%%% data - 2 dimentional (samples x electrodes) then electrode rejection
%%%% data - cell (trails x 1), always for trial-specific electrode rejection. Each cell (samples x electrodes).

if ~iscell(data)
    data_scads = data(samp_omit_scads:end-samp_omit_scads, :);
    [maxAmp, stdDev, gradient] = cog_scads_1_2(data_scads); % Get the editing matrices
else
    [maxAmp, stdDev, gradient] = cog_scads_1_2(data); % Get the editing matrices
end

lambda = 2; % as per the paper
    
if ~iscell(data)
    % Reject bad electrodes from complete time-series
    disp('---------------- REJECT SENSORS ----------------');
    
    rejMaxAmp = cog_scads_1_3(maxAmp, lambda, polar_ang);
    rejStdDev = cog_scads_1_3(stdDev, lambda, polar_ang);
    rejGrad = cog_scads_1_3(gradient, lambda, polar_ang);
    
    rejected_sensors = unique([rejMaxAmp; rejStdDev; rejGrad]);
    rejected_sensors = setdiff(rejected_sensors, [126; 127]);
    
%     for i = 1:size(data, 2)
%         for elec = 1:length(rejected_sensors)
%             data(:, rejected_sensors(elec)) = nan;
%         end
%     end
    data_out = data;
    data_out(:, rejected_sensors) = nan;
    
elseif iscell(data)
    % Reject trial specific bad electrodes
    disp('---------------- REJECT BAD TRIALS ----------------');
    
    rejMaxAmpEp = cog_scads_1_4(maxAmp, lambda, polar_ang);
    rejStdDevEp = cog_scads_1_4(stdDev, lambda, polar_ang);
    rejGradEp = cog_scads_1_4(gradient, lambda, polar_ang);
    
    rejected_sensors = [rejMaxAmpEp rejStdDevEp rejGradEp];
    
    for ix = 1:length(rejMaxAmpEp)
        rejected_sensors{ix, 1} = unique([rejMaxAmpEp{ix} rejStdDevEp{ix} rejGradEp{ix}]);
    end
    
    for elec = 1:length(rejected_sensors)
        for tr = 1:length(rejected_sensors{elec})
            data{rejected_sensors{elec}(tr)}(:, elec) = nan;
        end
    end
    
    size_tr = cellfun(@size, data, 'UniformOutput', false); size_tr = cat(1, size_tr{:});
    max_idx = max(size_tr, [], 1); max_length = max_idx(1); elecs = max_idx(2);
    data_out = nan(max_length, elecs, length(data));
    for ix = 1:length(data)
        data_out(1:size(data{ix},1), :, ix) = data{ix};
    end
end