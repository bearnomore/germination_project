function [incr, germination_time] = movie_strip_score_germination_red_MANUALLY(IMG, spore_num, germination_time_in, ...
    title_txt, ref_time)

incr = 1;

if isnan(germination_time_in)
    germination_time_in = 0;
end

max_fig = [1          41        1400         988];
%max_fig = [ 1          41        1680         934];
fig_out = figure('Position', max_fig);
set(gcf, 'MenuBar', 'none');



H = []; kkk = 0;
for ti = 1 : size(IMG,1)
   A =  double(IMG{ti,spore_num, 1});
   R =  double(IMG{ti,spore_num, 2});
   mid = median(1:size(A,1)); 
   if ti == 1
       minA = min(A(:));
       maxA = max(A(:));
       maxR = max(R(:));
       medR = median(R(:));
       stdR = std(R(:));
   end
    
   if isempty(A)
       if ti == 1
           germination_time = -1;
           incr = 1;
           close;
           return
       end
       break
   end
   kkk = kkk + 1;
  % A = A - min(A(:));
   A = A/maxA/2;
   
   R = R - medR - 0.5*stdR;
  R(R<0) = 0;
   R = R/(4*stdR);
   
   A = cat(3, A,A,A);
  A(:,:,1) = A(:,:,1) + 0.6*R;
%  A(:,:,1) = A(:,:,1) + 0.6*R;
   kalin_subplot(4,5,ti);
   H(ti) = imshow(A); axis square; axis off; 
   hold on; plot(mid, mid, 'g.')
   set(H(ti), 'UserData', 0); 
end

if nargin >= 4 & ~isempty(title_txt);
    title(title_txt)
end

%load germination_time_manual germination_time
t1 = germination_time_in;
germination_time = germination_time_in;

if t1 ~= 0
    axes(get(H(abs(t1)), 'Parent'));
    if  t1 > 0
        recH = rectangle('Position', [1, 1, 60, 60], 'EdgeColor', 'green', 'LineWidth',2);
    elseif t1<0
        recH = rectangle('Position', [1, 1, 60, 60], 'EdgeColor', 'red', 'LineWidth',2);
    end
    set(H(abs(t1)), 'UserData', recH);
end

for ti = 1 : kkk
    set(H(ti), 'ButtonDownFcn', {@get_time_point, ti, H});
end

if nargin >= 5
    axes(get(H(abs(ref_time)), 'Parent'));
    recH = rectangle('Position', [3, 3, 57, 57], 'EdgeColor', 'cyan', 'LineWidth',2);
end

set(gcf, 'KeyPressFcn',{@key_pressed});

waitfor(gcf,'UserData');

if strcmp(get(gcf, 'UserData'), '-1')
    incr = -1;
else
    germination_time = get(gcf, 'UserData');
end

close(gcf)



    function key_pressed(src,evnt)
        % load germination_time_manual germination_time
        switch evnt.Character
            %           case ' ' %complete inhibition (an alternative to left-clicking on the first square)
            %               germination_time(isolate)
            %
            %           case char(27) %esc %defective
            
            case    char(8) %backspace
                set(src, 'UserData', '-1');
        end
        
       

function get_time_point(src,  eventdat,  t, H)

%load germination_time_manual germination_time


%left click on first square is complete inhibition
%righ click on any square is defective well


if strcmp(get(gcf, 'SelectionType'), 'normal')
    
    set(gcf,'UserData', t);
   
    recH = rectangle('Position', [1, 1, 60, 60], 'EdgeColor', 'green', 'LineWidth',2);
    set(src, 'UserData', recH);
    
else
    set(gcf,'UserData', -t); % deffective well
    
    recH = rectangle('Position', [1, 1, 60, 60], 'EdgeColor', 'red', 'LineWidth',2);
    set(src, 'UserData', recH);
end

%save germination_time_manual germination_time

%set(gcf,'UserData', 'done');



