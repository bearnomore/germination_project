%% Germination curves of ISP5230dT, M145dT, B1511dT anf GrRd7dT
%  All spores were harvested after 7-10 days of incubation at 28 degree
%  All germinations were tested within 3 days after harvest on 0.02 NB HP
%  agar and at same apore densities

%% Load data 
R = {};


D_GrRd7 = 'Z:\Dropbox (Vetsigian lab)\Vetsigian lab Team Folder\Ye\interaction_HP\SefInhibition\GrRd7\20141201\Day3_40xD\';
load([D_GrRd7 '\results_with_manual_doublets'], 'results')
R{1} = results;
D_B1511 = 'Z:\Dropbox (Vetsigian lab)\Vetsigian lab Team Folder\Ye\interaction_HP\SefInhibition\B1511\20150123\Storage_sup_conditioned\rep1\';
load([D_B1511 '\results_with_manual_doublets'], 'results')
R{2} = results;
% ISP5230
D_ISP5230 = 'Z:\Dropbox (Vetsigian lab)\Vetsigian lab Team Folder\Ye\interaction_HP\SefInhibition\ISP5230\20141118\D1_119sp\';
load([D_ISP5230 '\results_with_manual_doublets'], 'results')
R{3} = results;

D_M145 = 'Z:\Dropbox (Vetsigian lab)\Vetsigian lab Team Folder\Ye\interaction_HP\SefInhibition\M145\20141110\D3_10xD\';
load([D_M145 '\results_with_manual_doublets'], 'results')
R{4} = results;



%% plot germination curve of each species with its germination rate and fitting
Prob_error = {};
Rate_error = {};
Germ_rate = {};
Germ_prob = {};
SporeDensity = [];
T = {};
plat_timepoint = [];
Nungerm = {}; 
NN = [];%total scored spore number
for kk =1:length(R)
    num_frames = length(R{kk}(end).Coords);
    GF = [R{kk}.IsolatedGerminationFrame];
    GF(GF<=0) = [];
    NN(kk) = length(GF); % total number of scored spores
    GF = GF(isfinite(GF)); % germination time of germinated spores
    dN = accumarray(GF', 1); %count number of germinated spores at each time point
    Ngerm = cumsum(dN); % accumulative count of numers of germinated spores along the timeline
    K = length(GF); % set up saturation number of germinated spores as the total number of germinated spores 
    % accumulative count of numbers of ungerminated spores along the timeline
    plat_timepoint = length(Ngerm);
    if plat_timepoint < num_frames
        Ngerm = [Ngerm; repmat(Ngerm(end), num_frames-plat_timepoint, 1)];
        dN = [dN; zeros(num_frames-plat_timepoint, 1)];
    end
    Germ_prob{kk} = Ngerm/NN(kk);
    Germ_rate{kk} = dN'/NN(kk);
    T{kk} = 1:length(Ngerm);
    for t = T{kk}
    % For the ungerminated spores at each time point (t), germination of ungerminated spores follows a binomial
    % distribution with the probability equal to the germination prob at that
    % point, Germ_prob
    Nungerm{kk}(t)= NN(kk) - Ngerm(t); %ungerminated spores at timepoint t
    tmp1=sort(binornd(NN(kk),Germ_prob{kk}(t), 1,1000)); 
    tmp2=sort(binornd(Nungerm{kk}(t),Germ_rate{kk}(t), 1,1000)); 
    % Set up a confidance interval of 80%. It could be st up to any other
    % range
%     Prob_error{kk}(t,:) = [tmp1(25) tmp1(975)]/NN(kk);
%     Rate_error{kk}(t,:) = [tmp2(25) tmp2(975)]/Nungerm{kk}(t);
%     
    Prob_error{kk}(t,:) = std(tmp1/NN(kk));%/sqrt(NN(kk));
    Rate_error{kk}(t,:) = std(tmp2/Nungerm{kk}(t));%/sqrt(NN(kk));
    
    end
    SporeDensity(kk) = round(mean([R{kk}.SporeDensity]));
    
end
%% Plot the germination curves with C.I = 95%
figure(1);
h1 = plot(T{1}, Germ_prob{1}, 'k-v','LineWidth', 2, 'markers', 12); hold on
plot(T{1}, Prob_error{1}, 'k--','Linewidth', 1.5);
h2 = plot(T{2}, Germ_prob{2}, 'k-^', 'LineWidth', 2, 'markers', 12);
plot(T{2}, Prob_error{2}, 'k--','Linewidth', 1.5);
h3 = plot(T{3}, Germ_prob{3}, 'k-o', 'LineWidth', 2, 'markers', 12)
plot(T{3}, Prob_error{3}, 'k--','Linewidth', 1.5);
h4 = plot(T{4}, Germ_prob{4}, 'k-sq', 'LineWidth', 2, 'markers', 12);
plot(T{4}, Prob_error{4}, 'k--','Linewidth', 1.5);
xlim([0 18])
xlabel('Germination time')
ylabel('Germination probability')
title('Germination of Streptomyces')
%% Plot with error bar
figure(2);
h1 = errorbar(T{1},Germ_prob{1}, Prob_error{1},  'k-v','LineWidth', 2, 'markers', 12); hold on
h2 = errorbar(T{2},Germ_prob{2}, Prob_error{2},  'k-^','LineWidth', 2, 'markers', 12) 
h3 = errorbar(T{3},Germ_prob{3}, Prob_error{3},  'k-o','LineWidth', 2, 'markers', 12) 
h4 = errorbar(T{4},Germ_prob{4}, Prob_error{4},  'k-sq','LineWidth', 2, 'markers', 12) 
xlim([0 18])
set(gca,'FontSize', 20, 'XTick', 0:2:18, 'YTick', 0:0.2:1)

%axis square
% xlabel('Time (hour)')
% ylabel('Germination fraction')
% title('Germination of Streptomyces')
% legend(sprintf('ISP5230 (%d sp/fv), N = %d', SporeDensity(1), NN(1)),sprintf('M145(%d sp/fv), N = %d', SporeDensity(2), NN(2)), sprintf('B1511(%d sp/fv), N = %d', SporeDensity(3), NN(3)), sprintf('GrRd7(%d sp/fv), N = %d', SporeDensity(4), NN(4)))
%legend(sprintf('ISP5230 (%d sp/fv)', SporeDensity(1)),sprintf('M145(%d sp/fv)', SporeDensity(2)), sprintf('B1511(%d sp/fv)', SporeDensity(3)), sprintf('GrRd7(%d sp/fv)', SporeDensity(4)))
%h = legend('1. {\itS.viridochromogenes} B1511', '2. {\itS.sp.} GrRd7', '3. {\itS.venezualae} ISP5230 ', '4. {\itS.coelicolor} M145', 'location', 'northoutside');
%LH = legend([h3 h4 h1 h2], sprintf('1. %s\n%s','{\itS.viridochromogenes}', 'B1511'), '2. {\itS.sp.} GrRd7', sprintf('3. %s\n%s', '{\itS.venezualae}', 'ISP5230'), '4. {\itS.coelicolor} M145')

% h = legend('{\itS.sp.} GrRd7','{\itS.viridochromogenes} B1511', '{\itS.venezualae} ISP5230 ', '{\itS.coelicolor} M145', 'location', 'northoutside')
% % set(h, 'FontSize', 20)
% legend boxoff
%%
figure(3); % B1511
[ax1, h1, h2]= plotyy(T{1},Germ_prob{1}, T{1}, Germ_rate{1}, @plot); 

% set(get(ax1(1), 'Ylabel'), 'String', 'Germination probability');
% set(get(ax1(2), 'Ylabel'), 'String', 'Germination rate');
set(h1, 'Color', 'k', 'LineWidth', 2.5, 'LineStyle', '-', 'Marker', '^', 'MarkerSize', 8)
set(h2, 'Color', 'r', 'LineWidth', 2.5, 'LineStyle', '-', 'Marker', '^', 'MarkerSize', 8)

set(ax1(1), 'YColor', 'k', 'FontSize', 24)
set(ax1(2), 'YColor', 'r', 'FontSize', 24)
set(ax1(1), 'YLim', [0 1])
set(ax1(2), 'YLim', [0 .6])
set(ax1(1), 'YTick', 0:0.2:1)
set(ax1(2), 'YTick', 0:0.2:.6)
set(ax1(1), 'XLim',[0 14])
set(ax1(2), 'XLim', [0 14])
set(ax1(1), 'XTick',0:2:14)
set(ax1(2), 'XTick', 0:2:14)
axes(ax1(1))
hold on; 
errorbar(T{1},Germ_prob{1}, Prob_error{1},'k-', 'LineWidth', 1.2)
% xlabel('Time (hour)')
axes(ax1(2))
hold on; 
errorbar(T{1},Germ_rate{1}, Rate_error{1}, 'r-', 'LineWidth', 1.2)
%title('Germination probability and germination rate of ISP5230')
xlabel('Time (hour)')
%%
figure(4); % GrRd7
[ax1, h1, h2]= plotyy(T{2},Germ_prob{2}, T{2}, Germ_rate{2}, @plot); 
xlabel('Time (hour)')
% set(get(ax1(1), 'Ylabel'), 'String', 'Germination probability');
% set(get(ax1(2), 'Ylabel'), 'String', 'Germination rate');
set(h1, 'Color', 'k', 'LineWidth', 2.5, 'LineStyle', '-', 'Marker', 'v', 'MarkerSize', 8)
set(h2, 'Color', 'r', 'LineWidth', 2.5, 'LineStyle', '-', 'Marker', 'v', 'MarkerSize', 8)

set(ax1(1), 'YColor', 'k', 'FontSize', 24)
set(ax1(2), 'YColor', 'r', 'FontSize', 24)
set(ax1(1), 'YLim', [0 1])
set(ax1(2), 'YLim', [0 .8])
set(ax1(1), 'YTick', 0:0.2:1)
set(ax1(2), 'YTick', 0:0.2:.8)
set(ax1(1), 'XTick',0:2:18)
set(ax1(2), 'XTick', 0:2:18)
set(ax1(1), 'XLim',[0 18])
set(ax1(2), 'XLim', [0 18])
axes(ax1(1))
hold on; 
errorbar(T{2},Germ_prob{2}, Prob_error{2}, 'k--', 'LineWidth', 1.2)
% xlabel('Time (hour)')
axes(ax1(2))
hold on; 
errorbar(T{2}, Germ_rate{2},Rate_error{2}, 'r--', 'LineWidth', 1.2)
%title('Germination probability and germination rate of M145')
% xlabel('Time (hour)')
%% 
figure(5); % ISP5230
[ax1, h1, h2]= plotyy(T{3},Germ_prob{3}, T{3}, Germ_rate{3}, @plot); 
xlabel('Time (hour)')
% set(get(ax1(1), 'Ylabel'), 'String', 'Germination probability');
% set(get(ax1(2), 'Ylabel'), 'String', 'Germination rate');
set(h1, 'Color', 'k', 'LineWidth', 2.5, 'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 8)
set(h2, 'Color', 'r', 'LineWidth', 2.5, 'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 8)

set(ax1(1), 'YColor', 'k', 'FontSize', 24)
set(ax1(2), 'YColor', 'r', 'FontSize', 24)
set(ax1(1), 'YLim', [0 1])
set(ax1(2), 'YLim', [0 .4])
set(ax1(1), 'YTick', 0:0.2:1)
set(ax1(2), 'YTick', 0:0.1:0.4)
set(ax1(1), 'XLim',[0 16])
set(ax1(2), 'XLim', [0 16])
set(ax1(1), 'XTick',0:2:16)
set(ax1(2), 'XTick', 0:2:16)

axes(ax1(1))
hold on; 
errorbar(T{3},Germ_prob{3},Prob_error{3}, 'k--', 'LineWidth', 1.2)
% xlabel('Time (hour)')
axes(ax1(2))
hold on; 
errorbar(T{3},Germ_rate{3}, Rate_error{3}, 'r--', 'LineWidth', 1.2)
%title('Germination probability and germination rate of B1511')
% xlabel('Time (hour)')
%% 
figure(6); % M145
[ax1, h1, h2]= plotyy(T{4},Germ_prob{4}, T{4}, Germ_rate{4}, @plot); 
xlabel('Time (hour)')
% set(get(ax1(1), 'Ylabel'), 'String', 'Germination probability');
% set(get(ax1(2), 'Ylabel'), 'String', 'Germination rate');
set(h1, 'Color', 'k', 'LineWidth', 2.5, 'LineStyle', '-', 'Marker', 'sq', 'MarkerSize', 8)
set(h2, 'Color', 'r', 'LineWidth', 2.5, 'LineStyle', '-', 'Marker', 'sq', 'MarkerSize', 8)

set(ax1(1), 'YColor', 'k', 'FontSize', 24)
set(ax1(2), 'YColor', 'r', 'FontSize', 24)
set(ax1(1), 'YLim', [0 1])
set(ax1(2), 'YLim', [0 0.2])
set(ax1(1), 'YTick', 0:0.2:1)
set(ax1(2), 'YTick', 0:0.1:0.2)
set(ax1(1), 'XLim',[0 18])
set(ax1(2), 'XLim', [0 18])
set(ax1(1), 'XTick',0:2:18)
set(ax1(2), 'XTick', 0:2:18)

axes(ax1(1))
hold on; 
errorbar(T{4},Germ_prob{4},Prob_error{4}, 'k--', 'LineWidth', 1.2)
% xlabel('Time (hour)')
axes(ax1(2))
hold on; 
errorbar(T{4},Germ_rate{4}, Rate_error{4}, 'r--', 'LineWidth', 1.2)
% xlabel('Time (hour)')
%title('Germination probability and germination rate of GrRd7')

%% Generate germination curve using isolatedspore data
% timesS = {};
% curveS = {};
% numFrame = [];
% SporeDensity = [];
% 
% for i = 1:length(R)
%     % Singles: generate germination probabilities along timeline
%     numFrame(i) = max(length(R{i}(end).Coords), length(R{i}(8).Coords));
%     GF = [R{i}.IsolatedGerminationFrame];
%     GF(GF<=0) = []; %unscored ones should be zero. defective ones should be -1
%     
%     T2 = sort(unique(GF));
%     GerminatedSpores2 = [];
%     kkk = 1;
%     for ti = T2(1:end-1)
%         GerminatedSpores2(kkk) = sum(GF<=ti);
%         kkk = kkk+1;
%     end
% %indx2 = find(diff(GerminatedSpores2)./diff(T2(1:end-1))<50);
% %indx2 = 1 : length(GerminatedSpores2);
%     T2 = T2(1 : end-1);
%     timesS{i} = 1 : numFrame(i);
%     GermSp = zeros(1, length(timesS{i}));
%     GermSp(T2) = GerminatedSpores2; 
% 
%     for k = 2 : length(GermSp)
%         if GermSp(k) < GermSp(k-1)
%             GermSp(k) = GermSp(k-1);
%         end
%     end
%     curveS{i} = GermSp/length(GF);
%     SporeDensity(i) = round(mean([R{i}.SporeDensity]));
%     minS = min(timesS{i});
%     if minS~=0
%         timesS{i} = [0:(minS-1), timesS{i}];
%         curveS{i} = [zeros(1,minS), curveS{i}];
%     end
% end
% % 