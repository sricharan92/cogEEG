function visualizeTr(badTrDat, vsMat, hfMat, flagTr)

vsMat(vsMat > 0) = 1; 
hfMat(hfMat > 0) = 2; 
bothMat = vsMat + hfMat; 
bothMat(bothMat < 3) = 0; % Both conditions satisfied

figure; tr = 1; 
while tr <= length(badTrDat.trial)
    
    bothElec = find(bothMat(flagTr(tr), :) == 1); 
    stepElec = setdiff(find(vsMat(flagTr(tr), :) == 1), bothElec); 
    hfElec = setdiff(find(hfMat(flagTr(tr), :) == 2), bothElec);
    remElec = setdiff(1:128, unique([bothElec, stepElec, hfElec])); 
    d = badTrDat.trial{tr}'; 
    
    plot(d(:, remElec), 'k--'); hold on; 
    plot(d(:, stepElec), 'r'); hold on; 
    plot(d(:, hfElec), 'b'); hold on; 
    plot(d(:, bothElec), 'm'); hold on; 
    
    title(sprintf('Trial - %d/%d, Step - %d, HF - %d, BOTH - %d, Total - %d\n RED - STEP, BLUE - HF, MAG - BOTH',...
        tr, length(badTrDat.trial), length(stepElec), length(hfElec), length(bothElec), length(unique([bothElec, stepElec, hfElec])))); 
    xlabel('Sample'); ylabel('\muV');
    ylim([-150, 150]); 
    
    pause; 
    clf; 
    
    tr = tr + 1; 
    
end

end