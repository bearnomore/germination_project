function [IMG,spore_count,IMG_area,SporeDensity,Coords, GerminationTime, GerminationFrame, dist_to_nearest_neighbor, Area, times, Eccentricity, Ellipticity] = ...
    score_field_of_view_red11 (nd2_fname)

Zmin_red = 3; %std above the mean
Area_min_red = 10;
hg = fspecial('gaussian', [5 5], 0.8); %works great for spores
%hg = fspecial('gaussian', [3 3], 0.8);

max_offset = 20;

Coords = {}; 

max_fig = get(0, 'ScreenSize');
%figure;
IMG_ = {};
IMG_R_ = {};

rad = 40;

try
    meta = imreadBFmeta(nd2_fname); %failed in one case for unknown reason
catch
    disp('Error: meta = imreadBFmeta(nd2_fname);')
    SporeInfo = struct('AbsoluteCoords', {}, 'IsDoublet', {}, 'Area', {}, ...
        'Eccentricity', {}, 'Solidity', {}, ...
        'NumSporesInNeighborhood', {}, 'LocalCoords_AllNeighbors',  {});
    times = [];
    
    return
end
IMG_area=meta.width*meta.height;
num_frames = meta.nframes;
tmp =strmatch('Series 0 timestamp', meta.parameterNames);
T = meta.parameterValues(tmp);
T_=[];
for ti = 1 : length(T)
    T_(ti) = T(ti);
end
T_ = sort(T_);
times = T_;
times = times-times(1);

if length(T_) ~= num_frames
    error('time stamps')
end

R_ = imreadBF(nd2_fname, 1, 1:num_frames, 1);
A_ = imreadBF(nd2_fname, 1, 1:num_frames, 2);
LL = nan(num_frames,1); AA = cell(num_frames,1);
CC = cell(num_frames,1);
spore_id = cell(num_frames,1);
for ti = 1 : num_frames
    A = A_(:,:,ti);
    R = R_(:,:,ti);
    rg = imfilter(R,hg);
    
    %thresholding
    Thresh = median(R(:))+Zmin_red*std(R(:)); % this works great for identifying spores
    Thresh2 = median(R(:))+1.5*std(R(:)); %use weaker threshold for identifying elongated mycellium - large parameter relative the area
    
    r1 = rg;
    r1(r1<Thresh)=0;
    L = watershed(-r1);
    r1(L==0) = 0;
    
    r2 = rg;
    r2(r2<Thresh2)=0;
    
    BW = r1>Thresh;
    STATS = regionprops(BW, 'Centroid', 'Area', 'Solidity', 'Eccentricity', 'MajorAxisLength', 'Orientation', 'Perimeter');
    C = cat(1, STATS.Centroid);
    Ar=[STATS.Area];
    C(Ar<Area_min_red,:)=[];
    STATS(Ar<Area_min_red)=[];
    Ar(Ar<Area_min_red) = [];
    
    
    CC{ti} = C;
    
    if ti > 1
        %compute distance matrix to the previous time point
        
        
        C1 = CC{ti-1};
        
        
        D =pdist2(C1, C);
        
        [S,indx] = sort(D);
        
        spore_id{ti} = spore_id{ti-1}(indx(1,:));
        
        
        %this is with respect to t-1 - not universal ids
        indx_usage = accumarray(indx(1,:)',1);
        unused_indx = find (indx_usage==0);
        unused_ids = spore_id{ti-1}(unused_indx);
        conflict_id = spore_id{ti-1}(indx_usage>1);
        
        %There are cases when one spore splits into two.
        
        for s1 = 1 : length(conflict_id)
            pos = find(spore_id{ti}==conflict_id(s1));
            Sc = S(1,pos);
            Sc(Sc > max_offset) = nan;
            [~, tmp] = nanmin(Sc);
            spore_id{ti}(setdiff(pos, pos(tmp))) = nan;
        end
        
        spore_id{ti}(S(1,:) > max_offset) = nan;
        
        orphans = find(isnan(spore_id{ti}));
        if ~isempty(unused_indx)
            %try to match unused ids with nearby orphans
            % i.e., find orphans that are close to unused ids.
            D1 = D(unused_indx, orphans);
            [mD1, tmp] = min(D1, [], 1);
            tmp1 = tmp(mD1 < 20);
            O1 = orphans(mD1<20);
            if ~isempty(tmp1)
                Q = unused_ids(tmp1);
                mD1 = mD1(mD1 < 20);
                
                %problem: multiple orphans can find the same unused_id
                conflicts = find(accumarray(tmp1',1)>1);
                for s1 = 1 : length(conflicts)
                    pos = find(tmp1==conflicts(s1));
                    Sc = mD1(1,pos);
                    [~, qq] = nanmin(Sc);
                    Q(setdiff(pos, pos(qq))) = nan;
                end
                spore_id{ti}(O1) = Q;
            end
        end
        
        % how about remove all orphans for now
        bad = isnan(spore_id{ti});
        spore_id{ti}(bad) = [];
        CC{ti}(bad,:) = [];
        C(bad,:) = [];
        %       %remove orphans that are very close to existing spores
        %       %they are either spurious or indicate germinations
        %
        %       %remove orphans that are part of extended low-thres clusters
        %
        
        
        % how aboput orphans that match previous orphans (remove those too)
        
        
    else
        Area_ = numel(A);
        spore_id{1} = 1 : size(C,1);
    
        %print a validation image from ti=1;
        r2(r2>70) = 70;
        figure('Position', max_fig); imagesc(r2); axis equal; hold on; colormap(gray); axis ij;
        %set(gcf, 'KeyPressFcn', 'figure(gcf+1)')
        plot(C(:,1), C(:,2), 'r.');
        for qz = 1 : length(spore_id{ti})
            text(C(qz,1)+5, C(qz,2)+5, num2str(spore_id{ti}(qz)), 'Color', 'white', 'FontSize', 8)
        end
        title(ti)
    end
    
    %cutout the images
    Coords{ti} = zeros(length(spore_id{1}),2);
    for w = 1 : length(spore_id{1})
        
        ind = find(spore_id{ti}==w);
        if isempty(ind)
             Coords{ti}(w,:) = Coords{ti-1}(w,:);
        else            
            Coords{ti}(w,:) = C(ind,:);
        end

        
        x = round(Coords{ti}(w,1));
        y = round(Coords{ti}(w,2));
        
        try
            IMG_{ti,w} =  uint16(A(y-rad:y+rad, x-rad:x+rad));
            IMG_R_{ti,w} =uint16(R(y-rad:y+rad, x-rad:x+rad));
        catch
            IMG_{ti,w} =  [];
            IMG_R_{ti,w} = [];
        end
        
        %  save([ELEMENTS_dir num2str(w) '\T' num2str(ti)], ['T'  num2str(ti)])
    end
    
end
IMG = cat(3,IMG_, IMG_R_);
%===================================
%score germination events
%eventually, we can put this in the above loop

%hg = fspecial('gaussian', [5 5], 0.8); %this works fine
%hg = fspecial('gaussian', [5 5], 2);
hg = fspecial('gaussian', [5 5], 1);


GerminationFrame = zeros(1, length(spore_id{1}));
GerminationTime = Inf(size(GerminationFrame));
Ellipticity = cell(num_frames,1);
Eccentricity = cell(num_frames,1);
Area = cell(num_frames, 1);
for spore_num = 1 : length(spore_id{1})
    
    P = []; Ar = []; mal=[];
    for ti = 1 : size(IMG,1)
        
        sn = find(spore_id{ti}==spore_num);
        if isempty(sn)
            GerminationFrame(spore_num) = -ti;
            break
        end
        
        A =  double(IMG{ti,spore_num, 1});
        R =  double(IMG{ti,spore_num, 2});
        
        if isempty(A)
            GerminationFrame(spore_num) = -ti;
            break
        end
        mid = median(1:size(A,1));
        full = size(A,1);
        
        
        C = CC{ti};
        C0 = C(sn,:);
        

        D =pdist2(C, C0);
        
        neighbors = find(D<20 & D>0.001);
        
        C = bsxfun(@minus, C, C0) +mid;
        
        R = imfilter(R,hg);
        maxR = median(R(:))+6*std(R(:));
        R(R>maxR) = maxR;
        
        
        % %    Thresh = median(R(:))+3*std(R(:));
        % %    r1 = R;
        % %     r1(r1<Thresh)=0;
        % %     L = watershed(-r1);
        % %     r1(L==0) = 0;
        % %     r1 = r1>0;
        
        %Erase neighboring spores
        for z = 1 : length(neighbors)
            C1 = C(neighbors(z),:);
            
            %             Delta = C1-C0;
            %
            %             % if norm(Delta) > 7
            %             R(round(max(1,mid+Delta(2)-6): min(full, mid+Delta(2)+6)),...
            %                 round(max(1,mid+Delta(1)-6): min(full, mid+Delta(1)+6)))=0;
            %             %   end
            
            M = (C1 + mid)/2; % midpoint
            ev = C1 - mid; ev = ev/norm(ev); nv = zeros(1,2);
            nv(1) = -ev(2); nv(2) = ev(1);
            
            lambda = -8:8;
            for l = lambda
                Dline = round(M + nv.*l);
                R(Dline(2), Dline(1)) = 0;
            end
            
            R(mid,mid) = maxR;
            
            
        end
        
        % %     L = watershed(-R);
        % %     R(L==0) = 0;
        

        %  R = double(R > median(R(:))+2*std(R(:)));
        R = double(R > median(R(:))+2*std(R(:))); %both 3*std and 2*std work well
        
        %R=imclearborder(R);
        
        R = bwselect(R, mid, mid, 4);
        
        area_= bwarea(R);
%        Ar(ti) = sum(area_(:));
%        tmp = bwperim(R);
%        P(ti) = sum(tmp(:));
        
        STATS = regionprops(R,'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', 'Area');
        Eccentricity{ti}(spore_num) = STATS.Eccentricity;
        Ellipticity{ti}(spore_num) = STATS.MinorAxisLength/STATS.MajorAxisLength;
        Area{ti}(spore_num) = STATS.Area;
        if ~isempty(STATS)
            mal(ti) = STATS.MajorAxisLength;
            
            %  if ti > 3 & P(ti)/Ar(ti) > median(P./Ar)+4
            if ti > 3 && mal(ti) > median(mal(mal>0))+1.4
                GerminationFrame(spore_num) = ti;
                 %1) Too sudden an increase is scored as unreliable 
                if mal(ti) > 40
                     GerminationFrame(spore_num) = -ti;
                end
                break
            end
        else
             mal(ti) = 0; %ignore later
        end
              
       
    end
    
    
    
end


tmp = GerminationFrame~=0;
GerminationTime(tmp) = sign(GerminationFrame(tmp)).*times(abs(GerminationFrame(tmp)));


%%=============
%compute the distance to a nearest neigbour for all time points and get the
%smallest
%use only the first time point for now
 D = squareform(pdist(CC{1}));
 D(1:size(D,1)+1:end) = 10000;
dist_to_nearest_neighbor = min(D,[],1);
spore_count=length(spore_id{1});
SporeDensity=spore_count/IMG_area*10^6;

 
