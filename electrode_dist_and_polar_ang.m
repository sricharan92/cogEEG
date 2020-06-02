function [dist, polar_ang] = electrode_dist_and_polar_ang(system)
%{
layout_path - Full path, including filename and format.
system - Data acquisition system. Either 'Biosemi'or 'EGI'.
%}

dist = nan(128, 128); polar_ang = nan(128, 1);

if strcmpi(system, 'EGI')
    
    %elecLocs = readlocs('GSN-HydroCel-128.sfp');
    elecLocs = ft_read_sens('./GSN-HydroCel-128.sfp'); 
    for i = 1:128
        for j = 1:128
            %dist(i, j) = sqrt(power(elecLocs(i).X - elecLocs(j).X, 2) + power(elecLocs(i).Y - elecLocs(j).Y, 2) + power(elecLocs(i).Z - elecLocs(j).Z, 2));
            dist(i, j) = sqrt(sum(power(elecLocs.elecpos(i, :) - elecLocs.elecpos(j, :), 2))); 
        end
    end
    
    %{
    for i = 1:128
        polar_ang(i, 1) = elecLocs(i).sph_phi;
    end
    %}
elseif strcmpi(system, 'BioSemi')
    load('BioSemi_128_polarAng'); load('BioSemi_128_cartesian'); %#ok<LOAD>
    
    for i = 1:128
        for j = 1:128
            dist(i, j) = sqrt(sum(power(cartesian(i, :) - cartesian(j, :), 2))); %#ok<NODEF>
        end
    end
    
    polar_ang = sph_phi;
end

end