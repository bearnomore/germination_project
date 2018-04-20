% contains
% 1) manual scoring of doublets
% 2) manual scorind of a random subset of singles
%% INITIALIZE
%f = 'X:\Ye Xu\Microscopy\interaction_HP\VsM145\M145VsM145\SupConditioned\20140730\TrisControlM145\rep1\';
%f='X:\Ye Xu\Microscopy\interaction_HP\VsM145\Aste7VsM145\SupConditioned\20140826\TrisControlM145\rep1\'; %super dramatic
%f= '/Volumes/ws/sysbio/vetsigiangroup/shared/BACKUPS of Yes on VETSIGIAN03/130725_Ye/Ye Xu/Microscopy/interaction_HP/SefInhibition/ISP5230/20150227_dilut_condition/germsup/rep1/'
%f = '/Volumes/ws/sysbio/vetsigiangroup/shared/BACKUPS of Yes on VETSIGIAN03/130725_Ye/Ye Xu/Microscopy/interaction_HP/SefInhibition/ISP5230/20150306_sup_on0.02NBAgar/0.02NB/rep1/';
f = 'G:\Dropbox (Vetsigian lab)\Vetsigian lab Team Folder\Ye\interaction_HP\other codes\sample data\';
if exist([f 'results_with_manual_doublets.mat'])
    error('Careful: manual data already exists');
end
load ([f 'results11']);



Dbls = cell(length(results),1);
for fi = 1 : length(results)
    if isempty(results(fi).Coords)
        continue
    end
    C = results(fi).Coords{1};
    Y = pdist(C);
    %Z = linkage(C, 'single');
    Z = linkage(Y, 'average'); % try different methods 'single', 'ward', 'average' and verify dissimilarity by cophenet(Z,Y) 
    T = cluster(Z,'cutoff',15, 'criterion', 'distance');
    cSize = accumarray(T,1);
    cD = find(cSize==2);
    Dbls{fi} = zeros(length(cD),2);
    for di = 1 : length(cD)
        Dbls{fi}(di,:) = find(T==cD(di))';
    end
end



for fi = 1 : length(results)
    results(fi).Doublets = Dbls{fi};
    results(fi).DoubletsGerminationTime=nan(size(Dbls{fi}));
    results(fi).DoubletsGerminationFrame=nan(size(Dbls{fi}));
end

save([f 'results_with_manual_doublets'], 'results');
% save([f 'results_with_speed_BIG'],'results');
%% To start or continue manual scoring
%f='X:\Ye Xu\Microscopy\interaction_HP\VsM145\Aste7VsM145\SupConditioned\20140826\TrisControlM145\rep1\'; %super dramatic
%f='/Volumes/ws/sysbietsigiangroup/shared/BACKUPS of Yes on VETSIGIAN03/130725_Ye/Ye Xu/Microscopy/interaction_HP/SefInhibition/ISP5230/20150422ISPdTvsB1511eG/ISP5230dTcontrol/rep1/'
% load ([f 'results_with_manual_doublets']);
% load ([f 'IMGs11']);



load ([f 'results_with_manual_doublets']);
load ([f 'IMGs11']);
%% MANUAL SCORING

    %  find first non nan and not an empty image
    
    for fi = 1 : length(results)
        if isempty(results(fi).Coords)
            continue
        end
        indx = find(isnan(results(fi).DoubletsGerminationFrame(:,1)), 1, 'first');
%         indx = find(iszero(results(fi).DoubletsGerminationFrame(:,1)), 1, 'first');
        if isempty(indx)
            continue
        else
            break
        end
        
    end
    while 1 < 2
        [incr, GF] = movie_strip_score_doublet_red_MANUALLY(IMG{fi}, results(fi).Doublets(indx,:), results(fi).Coords, results(fi).DoubletsGerminationFrame(indx, :), [num2str(fi) ': ' num2str(indx)]);
        if incr ~= -1
            GT = sign(GF) .* results(fi).times(abs(GF));
            GT(GF==1) = Inf;
            GF(GF==1) = Inf;
            results(fi).DoubletsGerminationFrame(indx,:) = GF;
            results(fi).DoubletsGerminationTime(indx,:) = GT;
            
            save([f 'results_with_manual_doublets'], 'results');
%             save([f 'results_with_speed_BIG'],'results');
            
        end
        indx = indx + incr;
        
    end
    
 %% ==== manual score singles =====
 %=======================================
 % randomize the singles
 R = results; SS = [];
 for fi = 1 : length(results)
     d = find(R(fi).dist_to_nearest_neighbor > 60);
     SS =[SS [fi*ones(size(d)); d]]; 
 end 
 RandomizedSingles = randperm(size(SS,2));
 SS = SS(:, RandomizedSingles);
 %% score a random subset of singles
 %score 100 random spores
 % can safely interrupt at any time. To score more, run with "for i = 101 : 200" below,
 % or simply rerun the randomization module above.
 for i = 1:size(SS,2)
     fi = SS(1,i);
     sp = SS(2,i);
     [incr, GF] = movie_strip_score_germination_red_MANUALLY(IMG{fi}, sp, 0);
     %[incr, GF] = movie_strip_score_germination_green_MANUALLY(IMG{fi}, sp, 0);
        if incr ~= -1
            GT = sign(GF) .* results(fi).times(abs(GF));
            GT(GF==1) = Inf;
            GF(GF==1) = Inf;
            results(fi).IsolatedGerminationFrame(sp) = GF;
            results(fi).IsolatedGerminationTime(sp) = GT;
            
            save([f 'results_with_manual_doublets'], 'results');
        end
        i = i + incr;
 end
