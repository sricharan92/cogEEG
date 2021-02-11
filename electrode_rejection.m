function [data_out, rejected_sensors, trial_rejected_sensors] = electrode_rejection(data, events, events_req, polar_ang, block, samp_omit_scads)

% Reject bad electrodes and replace with NaNs
[data, rejected_sensors] = reject_sensors_SCADS(data, polar_ang, samp_omit_scads);

% Visual rejection
%[data, elec_str] = visual_electrode_rejection(data, rejected_sensors, polar_ang);
elec_str = []; 
rejected_sensors = [rejected_sensors elec_str];

%%%% ICA

% Epoch (create temp_data. Replace with NaNs after filtering)
[epochs] = epoch_data(data, events, events_req, block, 1);

% Identify trial-specific bad electrodes
[data_out, trial_rejected_sensors] = reject_sensors_SCADS(epochs, polar_ang);

%%%% Interpolation

end