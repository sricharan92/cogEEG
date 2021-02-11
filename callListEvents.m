subjects = [1, 4:8, 10:21]; 
dataPath = '/4tb/sricharan/EEG_Data/Processed_postChg_try2';
task = 'afc'; 
for s = subjects
    
    fprintf('============= SUBJECT %.2d =============== \n\n', s); 
    for block = 1:2
        
        fprintf('--------- BLOCK %.2d ---------- \n\n', block); 
        
        % Load events
        load(sprintf('%s/subject%.2d/%s/block%.2d/events%.2d_session_%.2d.mat', dataPath, s, task, block, s, 1)); 
        
        % List events
        listEventsEGI(cfgTr_ref.event, '_DINs');
    end
end
