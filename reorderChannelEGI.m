function dat = reorderChannelEGI(dat)

% Assumption1 -- electrode labels are of the form 'EXX' etc. 
elecArr = cellfun(@(x) str2num(x(2:end)), dat.label); 
[~, ind] = sort(elecArr, 'asc'); 
dat.trial = cellfun(@(x) x(ind, :), dat.trial, 'UniformOutput', false);
dat.label = dat.label(ind); 

end