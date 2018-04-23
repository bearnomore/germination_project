%% Loading data
R = {};
D3 = 'S:\Dropbox (Vetsigian lab)\Vetsigian lab Team Folder\Ye\interaction_HP\SefInhibition\ISP5230\20141119\D2_237sp\';
load([D3 '\results_with_manual_doublets'], 'results')
R{3} = results;
% D1 = 'S:\Dropbox (Vetsigian lab)\Vetsigian lab Team Folder\Ye\interaction_HP\SefInhibition\ISP5230\20141119\D2_233sp\';
% load([D1 '\results_with_manual_doublets'], 'results')
% R{1} = results;

D2 = 'S:\Dropbox (Vetsigian lab)\Vetsigian lab Team Folder\Ye\interaction_HP\SefInhibition\ISP5230\20141119\D2_155sp\';
load([D2 '\results_with_manual_doublets'], 'results')
R{2} = results;

D1 = 'S:\Dropbox (Vetsigian lab)\Vetsigian lab Team Folder\Ye\interaction_HP\SefInhibition\ISP5230\20141119\D2_71sp\';
load([D1 '\results_with_manual_doublets'], 'results')
R{1} = results;

%% old ISP5230 spores on 0.1NB
R = {};
D1 = 'S:\Dropbox (Vetsigian lab)\Vetsigian lab Team Folder\Ye\interaction_HP\SefInhibition\ISP5230\20150618ISP5230dT\rep1\';
load([D1 '\results_with_manual_doublets'], 'results')
R{1} = results;

D1 = 'S:\Dropbox (Vetsigian lab)\Vetsigian lab Team Folder\Ye\interaction_HP\SefInhibition\ISP5230\20150618ISP5230dT\rep2\';
load([D1 '\results_with_manual_doublets'], 'results')
R{1} = [R{1} results];

D1 = 'S:\Dropbox (Vetsigian lab)\Vetsigian lab Team Folder\Ye\interaction_HP\SefInhibition\ISP5230\20150618ISP5230dT\rep3\';
load([D1 '\results_with_manual_doublets'], 'results')
R{1} = [R{1} results];

D1 = 'S:\Dropbox (Vetsigian lab)\Vetsigian lab Team Folder\Ye\interaction_HP\SefInhibition\ISP5230\20150618ISP5230dT\rep4\';
load([D1 '\results_with_manual_doublets'], 'results')
R{1} = [R{1} results];

D1 = 'S:\Dropbox (Vetsigian lab)\Vetsigian lab Team Folder\Ye\interaction_HP\SefInhibition\ISP5230\20150618ISP5230dT\rep5\';
load([D1 '\results_with_manual_doublets'], 'results')
R{1} = [R{1} results];

D1 = 'S:\Dropbox (Vetsigian lab)\Vetsigian lab Team Folder\Ye\interaction_HP\SefInhibition\ISP5230\20150618ISP5230dT\rep6\';
load([D1 '\results_with_manual_doublets'], 'results')
R{1} = [R{1} results];
%% For isolated spores
N = [];
T ={};
Germ_prob = {};
Spore_density =[];

for k= 1:length(R)
    Spore_density(k) = mean([R{k}.SporeDensity]);
    num_frames = length(R{k}(1).Coords);
    T{k} = 1:num_frames;
    GF = [R{k}.IsolatedGerminationFrame];
    GF(GF<=0) = []; %unscored ones should be zero. defective ones should be -1
    N(k) = length(GF); % total num of scored spores
    GF = GF(isfinite(GF)); % germ frames of germinated spores
    dN = accumarray(GF',1); % num of newly germinated spores at each time point
    plat_time = length(dN); % timepoint when germination reach stationary phase
    if plat_time ~= num_frames
        dN = [dN; zeros(num_frames-plat_time,1)];
    end
    Ngerm = cumsum(dN); % accumulative num of spores along timeline
    Nungerm = N(k) - Ngerm; %num of ungerminated spores at each timepoint
    Germ_prob{k} = Ngerm/N(k);   
    imgA = [R{k}.IMG_area];
    sp_count = [R{k}.spore_count];
    Spore_density(k) = mean(sp_count./(imgA*0.0225));%1px^2 = 0.0225um^2
       
end

%% For touching spores
L =[];
P00 = {};
P10 = {};
P11 = {};
GFd00 = {};
GFd11 = {};
GFd10 = {};

for k = 1:length(R)
    GFd = cat(1, R{k}.DoubletsGerminationFrame);
    bad = any (GFd<0 | isnan(GFd),2);
    GFd(bad,:) = [];
    L(k) = size(GFd,1); %num of pairs of touching spores
    germID = sum(isfinite(GFd),2)~=0; %find spore pair of which at least one germinated
    ungermID = sum(isfinite(GFd),2)==0; %find spore pair of which neither germinated
    GFd00{k} = GFd(ungermID,:);%neither germinated spore pairs
    GFd_ = GFd(germID,:); %at least one germinated
    GFd11{k} = GFd_(isfinite(GFd_(:,2)),:); %both germinated spore pairs
    GFd10{k} = GFd_(~isfinite(GFd_(:,2)),:); %one germinated spore pairs
    kkk = 1;
    for ti = T{k}
        z = sum(GFd<=ti,2);
        P00{k}(kkk) = sum(z==0)/L(k); %zero germinated
        P10{k}(kkk) = sum(z==1)/L(k); %one germinated
        P11{k}(kkk) = sum(z==2)/L(k); %both germinated
        kkk = kkk+1;
    end

end

%% Calsulate prob errors for isolated spore pairs
Err11_iso = {};
Err11_iso_ = {};
Err10_iso = {};
Err10_iso_ = {};
Err00_iso = {};
Err00_iso_ = {};

for k = 1:length(R)
    for kk = T{k}
        E11 = [];
        for kkk = 1:1000
            tmp11 = binornd(N(k), Germ_prob{k}(kk));
            tmp11_ = binornd(N(k)-tmp11, Germ_prob{k}(kk));
            E11(kkk) = tmp11*tmp11_/(N(k)*(N(k)-tmp11));
        end
        E11 = sort(E11);
        Err11_iso{k}(kk,:) = [E11(50) E11(950)];
        Err11_iso_{k}(kk,:) = std(E11);%/sqrt(N(k)); %s.e.m
        
        E10 = [];
        for kkk = 1:1000
            tmp10 = binornd(N(k), Germ_prob{k}(kk));
            tmp10_ = binornd(N(k)-tmp10,1-Germ_prob{k}(kk));
            E10(kkk) = 2*(tmp10*tmp10_)/(N(k)*(N(k)-tmp10));
        end
        E10 = sort(E10);
        Err10_iso{k}(kk,:) = [E10(50) E10(950)];
        Err10_iso_{k}(kk,:) = std(E10);%/sqrt(N(k)); %s.e.m
        
        E00 = [];
        for kkk = 1:1000
            tmp00 = binornd(N(k), 1-Germ_prob{k}(kk));
            tmp00_ = binornd(N(k)-tmp00,1-Germ_prob{k}(kk));
            E00(kkk) = (tmp00*tmp00_)/(N(k)*(N(k)-tmp00));
        end
        E00 = sort(E00);
        Err00_iso{k}(kk,:) = [E00(50) E00(950)];
        Err00_iso_{k}(kk,:) = std(E00);%/sqrt(N(k)); %s.e.m
    
    end
end

%% Calculate prob errors for touching spore pairs
Err11_touch = {};
Err11_touch_ = {};
Err10_touch = {};
Err10_touch_ = {};
Err00_touch = {};
Err00_touch_ = {};

for k = 1:length(R)
    for kk = T{k}
        tmp11 = sort(binornd(L(k), P11{k}(kk), 1,1000));
        Err11_touch{k}(kk,:) = [tmp11(50) tmp11(950)]/L(k);
        Err11_touch_{k}(kk,:) = std(tmp11/L(k));%/sqrt(L(k)); %s.e.m.
        
        tmp10 = sort(binornd(L(k), P10{k}(kk), 1, 1000));
        Err10_touch{k}(kk,:) = [tmp10(50) tmp10(950)]/L(k);
        Err10_touch_{k}(kk,:) = std(tmp10/L(k));%/sqrt(L(k)); %s.e.m
        
        tmp00 = sort(binornd(L(k), P00{k}(kk), 1,1000));
        Err00_touch{k}(kk,:) = [tmp00(50) tmp00(950)]/L(k);
        Err00_touch_{k}(kk,:) = std(tmp00/L(k));%/sqrt(L(k)); %s.e.m.
    end
end
%% plot only 11
for k = 1:length(R)
    subplot(1,3,k)
    % plot(T{k}, Germ_prob{k}.^2, 'bo-', 'LineWidth', 2); hold on
    % plot(T{k}, Err11_iso{k}, 'b--', 'LineWidth', 1.5);
    % plot(T{k}, P11{k}, 'ro-', 'LineWidth',2);
    % plot(T{k}, Err11_touch{k}, 'r--','LineWidth', 1.5); hold off
    errorbar(T{k}, Germ_prob{k}.^2,Err11_iso_{k}, 'ko-', 'LineWidth', 2, 'markers', 8); hold on
    errorbar(T{k}, P11{k}, Err11_touch_{k}, 'ro-', 'LineWidth', 2, 'markers', 8);
    set(gca, 'FontSize', 15, 'XTick', 0:2:16);
    ylim([0 0.5])
    xlabel('Time(hours)')
    ylabel('Germination fraction')
    %legend(['Isolate (' sprintf('N = %d)', N(k))],['Touch (' sprintf('N = %d)', L(k))])
    legend('Isolated','Touching','Location', 'northwest')
    title({[num2str(round(100*Spore_density(k),2)), 'spores/', '100\mum','{^2}']},'Interpreter','Tex')
end

%% Plot

for k = 1:length(R)
    figure(k)
    subplot(1,3,1)
    % plot(T{k}, Germ_prob{k}.^2, 'bo-', 'LineWidth', 2); hold on
    % plot(T{k}, Err11_iso{k}, 'b--', 'LineWidth', 1.5);
    % plot(T{k}, P11{k}, 'ro-', 'LineWidth',2);
    % plot(T{k}, Err11_touch{k}, 'r--','LineWidth', 1.5); hold off
    errorbar(T{k}, Germ_prob{k}.^2,Err11_iso_{k}, 'bo-', 'LineWidth', 2, 'markers', 8); hold on
    errorbar(T{k}, P11{k}, Err11_touch_{k}, 'ro-', 'LineWidth', 2, 'markers', 8);
    set(gca, 'FontSize', 25, 'XTick', 0:2:16);
    ylim([0 0.6])
    xlabel('Time point (hours)')
    ylabel('Germination probability')
    %legend(['Isolate (' sprintf('N = %d)', N(k))],['Touch (' sprintf('N = %d)', L(k))])
    legend('Isolated','Touching')
    title({'11',[num2str(round(Spore_density(k),4)), 'spores/', '\mum','{^2}']},'Interpreter','Tex')
    
    subplot(1,3,2)
    % plot(T{k}, 2*Germ_prob{k}.*(1-Germ_prob{k}), 'bo-', 'LineWidth', 2); hold on
    % plot(T{k}, Err10_iso{k}, 'b--', 'LineWidth', 1.5);
    % plot(T{k}, P10{k}, 'ro-', 'LineWidth',2);
    % plot(T{k}, Err10_touch{k}, 'r--','LineWidth', 1.5); hold off
    errorbar(T{k}, 2*Germ_prob{k}.*(1-Germ_prob{k}),Err10_iso_{k}, 'bo-', 'LineWidth', 2, 'markers', 8); hold on
    errorbar(T{k}, P10{k}, Err10_touch_{k}, 'ro-', 'LineWidth', 2, 'markers', 8);
    set(gca, 'FontSize', 25, 'XTick', 0:2:16);
    ylim([0 0.6])
    xlabel('Time point (hours)')
    ylabel('Germination probability')
    %legend(['Isolate (' sprintf('N = %d)', N(k))], ['Touch (' sprintf('N = %d)', L(k))])
    legend('Isolated','Touching')
    title({'10', [num2str(round(Spore_density(k), 4)),'spores/', '\mum','{^2}']},'Interpreter','Tex')
    
    subplot(1,3,3)
    % plot(T{k}, (1-Germ_prob{k}).^2, 'bo-', 'LineWidth', 2); hold on
    % plot(T{k}, Err00_iso{k}, 'b--', 'LineWidth', 1.5);
    % plot(T{k}, P00{k}, 'ro-', 'LineWidth',2);
    % plot(T{k}, Err00_touch, 'r--','LineWidth', 1.5); hold off
    errorbar(T{k}, (1-Germ_prob{k}).^2, Err10_iso_{k}, 'bo-', 'LineWidth', 2, 'markers', 8); hold on
    errorbar(T{k}, P00{k}, Err00_touch_{k}, 'ro-', 'LineWidth', 2, 'markers', 8);
    set(gca, 'FontSize', 25,'XTick', 0:2:16);
    ylim([0 1])
    xlabel('Time point (hours)')
    ylabel('Germination probability')
    %legend(['Isolate (' sprintf('N = %d)', N(k))], ['Touch (' sprintf('N = %d)', L(k))])
    legend('Isolated','Touching')
    title({'00', [num2str(round(Spore_density(k),4)),'spores/', '\mum','{^2}']},'Interpreter','Tex')

end
