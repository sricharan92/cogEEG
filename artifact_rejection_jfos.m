function [trial_idx] = artifact_rejection_jfos(data_in, sf)
%% Input
% data_in: EEG data trials' cell array. Each trial - electrodes x samples
% sf     : Downsampling factor

%% Output
% trial_idx: Rejected trial indices

%% Code
% Identifies EEG artifacts as described in Joshua J. Foster et al., 2020
tic;
% Voltage drift rejections
vd_rejects = zeros(length(data_in), size(data_in{1}, 1));
for tr = 1:length(data_in)
    for el = 1:size(data_in{1}, 1)
        tr_len = size(data_in{tr}, 2); % Length of trial
        quart1 = mean(data_in{tr}(el, 1:round(tr_len/4))); % Mean value of first quarter of trial
        quart4 = mean(data_in{tr}(el, round(3*tr_len/4):end)); % Mean value of last quarter of trial
        if abs(quart4 - quart1) > 40 % Check if absolute difference exceeds 40 uV
            vd_rejects(tr, el) = 1;
        end
    end
end
vd_rejects_idx = find(sum(vd_rejects, 2) ~= 0); % Indices of trials with at least one rejected electrode
toc;

% Sudden voltage step rejections
tic;
win_len = round(250/sf); % Length of window in samples
win_stp = round(20/sf); % Step size in samples
vs_rejects_marg = []; vs_rejects = {}; 
for tr = 1:length(data_in)
    disp(tr);
    no_of_steps = floor((size(data_in{tr}, 2) - win_len + 1)/win_stp); % Number of time steps for the sliding window
    vs_rejects{tr} = zeros(size(data_in{1}, 1), no_of_steps);
    for el = 1:size(data_in{1}, 1)
        t1 = 1;
        for ti = 1:no_of_steps
            dat = data_in{tr}(el, t1:t1+win_len-1); % Data in current window
            half1 = mean(dat(1:floor(win_len/2))); % Mean value of first half of window
            half2 = mean(dat(floor(win_len/2)+1:end)); % Mean value of second half of window
            if abs(half2 - half1) > 60 % Check if absolute difference exceeds 60 uV
                vs_rejects{tr}(el, ti) = 1;
            end
            t1 = t1 + win_stp;
        end 
    end
    % For this trial, does any time point for any elec exceed the threshold?
    vs_rejects_marg(tr, 1) = nansum(nansum(vs_rejects{tr}, 2), 1); 
end
vs_rejects_idx = find(vs_rejects_marg ~= 0); % Indices of trials with at least one rejected electrode
toc; 

% High frequency noise rejections
tic;
win_len = round(15/sf); % Length of window in samples
win_stp = round(15/sf); % Step size in samples
hf_rejects_marg = []; hf_rejects = {}; 
for tr = 1:length(data_in)
    no_of_steps = floor((size(data_in{tr}, 2) - win_len + 1)/win_stp); 
    hf_rejects{tr} = zeros(size(data_in{1}, 1), no_of_steps);
    for el = 1:size(data_in{1}, 1)
        t1 = 1;
        for ti = 1:no_of_steps
            dat = data_in{tr}(el, t1:t1+win_len-1); % Data in current window
            ext1 = min(dat); % Min value of window
            ext2 = max(dat); % Max value of window
            if abs(ext2 - ext1) > 120 % Check if absolute difference exceeds 120 uV
                hf_rejects{tr}(el, ti) = 1;
            end
            t1 = t1 + win_stp;
        end
    end
    hf_rejects_marg(tr, 1) = nansum(nansum(hf_rejects{tr}, 2), 1); 
end
hf_rejects_idx = find(hf_rejects_marg ~= 0); % Indices of trials with at least one rejected electrode
toc; 

trial_idx = unique([vd_rejects_idx; vs_rejects_idx; hf_rejects_idx]);

end