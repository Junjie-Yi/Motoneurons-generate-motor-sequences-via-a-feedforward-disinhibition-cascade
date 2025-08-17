clear
MyColor = [62/255 43/255 109/255;
    0 0.45 0.74;
    0.49 0.18 0.56;
    48/255 151/255 164/255;
    240/255 100/255 73/255;
    255/255 170/255 50/255;
    5/255 80/255 91/255];

MyColorHot = [0.85 0.33 0.10;
    240/255 100/255 73/255;
    0.93 0.69 0.13;
    255/255 187/255 0];
MyColorCold = [0.49 0.18 0.56;
    5/255 80/255 91/255;
    48/255 151/255 164/255;
    0 0.45 0.74];
MyColorGray = [0 0 0;
    70/255 70/255 70/255;
    140/255 140/255 140/255;
    210/255 210/255 210/255];


clear res
load('fig1c_data.mat');
% mode 0 = plot single fly; 1 = plot population average
mode = 1;
i_bin = 0.995:-0.01:0.005;
i_bin = fliplr(i_bin);
t_bin_all = zeros(length(res),length(i_bin)*2);
t_bin_diff_all = zeros(length(res),length(i_bin)*2-1);
period = zeros(1,length(res));
phaseDuration_all = zeros(length(res),3);
i_min = zeros(1,length(res));
% figure
if mode == 1
    fig1 = figure('Units', 'centimeters', 'Position', [20, 10, 3, 2.5]);hold on
    set(gca, 'Position', [0.3, 0.2, 0.5, 0.5])
end

for jj = 1:length(res)
    p = res(jj).localMax;
    clear i;
    i = res(jj).I;

    period_tmp = diff(p);
    res(jj).period = period_tmp;
    period_mean_tmp = mean(period_tmp);
    period(jj) = period_mean_tmp;
    
    t_bin = zeros(length(p)-1,2*length(i_bin));
    t_bin_diff = zeros(length(p)-1,2*length(i_bin)-1);
    phaseDuration = zeros(length(p)-1,3);
    
    i_range = 0;
    range = zeros(1,length(p));
    for ii = 1:length(p)-1
        b = p(ii);e = p(ii+1);
        b_idx = find(res(jj).t == b);e_idx = find(res(jj).t == e);
        t = res(jj).t(b_idx:e_idx) - res(jj).t(b_idx); i_single = i(b_idx:e_idx);
        range = (max(i_single) - min(i_single));
        i_range = i_range + (max(i_single) - min(i_single))/(length(p)-1);
    end
    
    
    i_min_tmp  = zeros(1,length(p)-1);
    for ii = 1:length(p)-1
        b = p(ii);e = p(ii+1);
        b_idx = find(res(jj).t == b);e_idx = find(res(jj).t == e);
        t = res(jj).t(b_idx:e_idx) - res(jj).t(b_idx); i_single = i(b_idx:e_idx);

        i_norm = 1-(i_single - min(i_single)) / i_range;
        t_norm = t/max(t);
        
        t_inter1 = 0:0.001:max(t);
        i_inter1 = interp1(t,i_norm,t_inter1,'pchip');
        i_min_tmp(ii) = i_inter1(1);
        xidx = zeros(1,2*length(i_bin));
        for  kk = 1:length(i_bin)
            tmp = find(i_inter1>i_bin(kk));
            xidx(kk) = tmp(1);
            xidx(end-kk+1) = tmp(end);
        end
        t_bin(ii,:) = t_inter1(xidx);
        t_bin_diff(ii,:) = diff(t_inter1(xidx));  
    end
    
    i_min(jj) = mean(i_min_tmp);
    t_bin_mean = mean(t_bin);
    t_bin_diff_mean = mean(t_bin_diff);

    t_sum = zeros(1,length(t_bin_diff_mean));
    for ii = 1:length(t_sum)
        t_sum(ii) = sum(t_bin_diff_mean(1:ii));
    end
    
    if mode == 1
        c = MyColorGray(1,:);
        y = [i_bin(1:end),flip(i_bin(1:end-1))];
        if jj == 2
            h(jj) = plot(t_sum,[i_bin(1:end),flip(i_bin(1:end-1))],'color',[0 0.45 0.87],'linewidth',0.8);
        else
            h(jj) = plot(t_sum,[i_bin(1:end),flip(i_bin(1:end-1))],'color',c,'linewidth',0.5);
        end
        ax = gca;
        set(ax,'color','none')
        set(ax,'box','off')
    end
    
end