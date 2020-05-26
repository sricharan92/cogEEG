function [data, elec_str] = visual_electrode_rejection(data, rejected_sensors, polar_ang)
%%%% data - electrodes x samples

disp('---------------- MANUAL ELECTRODE REMOVAL ----------------'); 
disp('You will see 3 figures. Based on those, choose the electrodes to be excluded'); 
disp('Once done with all the electrodes, you will be asked about the electrodes you want removed. Enter the number with SPACES in between'); 
disp('For example:'); 
disp('Enter electrodes to remove: 1 10 25 50'); 

% Get the editing matrices
elecs = 1:128; % might want to make this global
editing_matrix = [];
[editing_matrix(:, :, 1), editing_matrix(:, :, 2), editing_matrix(:, :, 3)] = cog_scads_1_2(data);

med_sensor = []; para_str = {'Maximum Amplitude', 'Standard Deviation', 'Gradient'}; 
for i = 1:3 
    med_sensor(:, i) = getCorrectedMedians(editing_matrix(:, :, i), polar_ang(elecs), 1); %#ok<AGROW>
    med_sensor(rejected_sensors, i) = NaN; % Since for the standard deviation measure (nanstd from fieldtrip), it takes 0 for a vector of NaNs
    
    figure;
    plot(polar_ang(elecs), med_sensor(:, i), '.'); hold on; 
    title(sprintf('%s editing matrix after average rereferencing', para_str{i})); 
    xlabel('Polar Angle'); ylabel('Editing matrix values after de-trending'); 
    text(polar_ang(elecs), med_sensor(:, i), strsplit(num2str(elecs)), 'FontSize', 14); % if it doesn't work take transpose of elecs
end

% Display and remove electrodes
figure; 
title('Trial mean time course data for each electrode'); 
xlabel('Samples'); ylabel('microVolts'); 
for i = 1:length(elecs)
    plot(data(:, i)); hold on;
    fprintf('Electrode plotted - %d \n\n', i);
    pause;
end

elec_str = input('Enter electrodes to remove: ', 's');
elec_str = str2num(elec_str)';  %#ok<ST2NM>
elec_str = setdiff(elec_str, [126; 127]);

rejected_sensors = [rejected_sensors; elec_str];  %#ok<NASGU>

% for elec = 1:length(elec_str)
%     data(:, elec_str(elec)) = nan;
% end
data(:, elec_str) = nan;

end