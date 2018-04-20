function [results, return_code] = score_well_red11(dir_name)

if isunix || ismac
    slash = '/';
else
    slash = '\';
end

if dir_name(end)~=slash
    dir_name = [dir_name slash];
end

aligned_nd2_ = dir([dir_name '*Aligned.nd2']);
if isempty(aligned_nd2_)
    return_code = -1;
    return
end
    
IMG = {};
spore_count={};
Coords = {}; GerminationTime={}; GerminationFrame = {};
dist_to_nearest_neighbor = {};
Area = {}; times = {};
%AverageSporeCount={};
IMG_area={};
SporeDensity={}
Eccentricity = {};
Ellipticity = {};
results = struct;

for pos = 1 : length(aligned_nd2_) %%%%%%%%%%%%%%%%
    xy = aligned_nd2_(pos).name;
    tmp = findstr(xy, ' - Aligned');
    xy1 = xy(1:tmp(1)-1);
    
    results(pos).field_name = xy1;
    
    try
        
        [IMG{pos},results(pos).spore_count,results(pos).IMG_area,results(pos).SporeDensity,results(pos).Coords, results(pos).GerminationTime, results(pos).GerminationFrame,...
            results(pos).dist_to_nearest_neighbor, results(pos).Area, results(pos).times,results(pos).Eccentricity, results(pos).Ellipticity] = ...
            score_field_of_view_red11([dir_name xy]);
        
      
        save([dir_name 'results11'], 'results');
        save([dir_name 'IMGs11'], 'IMG');
    end
end

save([dir_name 'results11'], 'results');
save([dir_name 'IMGs11'], 'IMG'); 