function score_all11(DD, forced_redo)
%score all sub-directories containing alighned nd2 files
%DD is the parent directory

if nargin < 2
    forced_redo = 0; %do not redo by default
end

disp(DD);

if isunix || ismac
    slash = '/';
else
    slash = '\';
end

if DD(end) ~= slash
    DD = [DD slash];
end

%check if the directory contains nd2 files
aligned_nd2_ = dir([DD '*Aligned.nd2']);
if isempty(aligned_nd2_)
   A = dir(DD);
   A = A([A.isdir]);
   A(1:2) = [];
   for di = 1 : length(A)
       score_all11([DD A(di).name], forced_redo);
   end
else
    if (forced_redo==1) || ~exist([DD 'result.mat'])
        %score_well_red(DD);
        score_well_red11(DD);
        close all
    end
end




