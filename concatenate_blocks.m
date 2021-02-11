close all;
clear all; %#ok<CLALL>
% clc;

dbstop if error;

%% Give the subject and block details here
subjects = 20;
dashes = '----------------';

%% Add all relevant toolboxes to path
hd = mfilename('fullpath'); ix = find(hd == '\'); home_dir = hd(1:ix(end-2)); % Home Directory
addpath(genpath([home_dir 'Scripts\'])); % Scripts

data_dir = [home_dir 'Data\Processed\'];

%% Loop for all subjects
for s = subjects
    disp([dashes ' Subject ' num2str(s, '%.2d') ' ' dashes])
    subject_dir = [data_dir 'subject' num2str(s, '%.2d') '\Preprocessed\'];
    block_names = dir([subject_dir '1_S*']);
    for block = 1:length(block_names)
        disp([dashes ' Block ' num2str(block, '%.2d') ' ' dashes])
        load([subject_dir block_names(block).name], 'data_out')
        
        % Concatenate blocks
        if ~exist('data', 'var')
            data = data_out;
        else
            if size(data, 1) > size(data_out, 1)
                data_out(end:size(data, 1), :, :) = nan;
            elseif size(data, 1) < size(data_out, 1)
                data(end:size(data_out, 1), :, :) = nan;
            end
            data = cat(3, data, data_out);
        end
    end
    save([subject_dir '\2_S' num2str(s)], 'data');
    clear data
end