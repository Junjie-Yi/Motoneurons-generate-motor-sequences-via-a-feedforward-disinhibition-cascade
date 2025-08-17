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
load('fig1a,b_data.mat');
i_fh = array{1,9};
i_fh = i_fh(3*120:10*120);
i_lh = array{1,10};
i_lh = i_lh(3*120:10*120);
localMin = islocalmin(i_fh,'MinSeparation',10,'FlatSelection', 'first','MinProminence',3);
localMax = islocalmax(i_fh,'MinSeparation',10,'FlatSelection', 'first','MinProminence',3);
tf_lmin = find(localMin>0);
tf_lmax = find(localMax>0);
i_lmin = i_fh(localMin);
i_lmax = i_fh(localMax);
res.i = i_fh;
res.i_lh = i_lh;
t =  1/120:1/120:length(i_fh)/120;
res.t = 1/120:1/120:length(i_fh)/120;
res.localMax = t(tf_lmax);
res.tf_lmin = tf_lmin;
res.i_lmin = i_lmin;
res.i_lmax = i_lmax;


%%
% mode 0 = plot single fly; 1 = plot population average
mode = 0;
i_bin = 0.995:-0.01:0.005;
i_bin = fliplr(i_bin);
t_bin_all = zeros(length(res),length(i_bin)*2);
t_bin_diff_all = zeros(length(res),length(i_bin)*2-1);
period = zeros(1,length(res));
phaseDuration_all = zeros(length(res),3);
i_min = zeros(1,length(res));
% figure
for jj = 1:length(res)
    p = res(jj).localMax;
    clear i;
    i_fh = res(jj).i;
    i_lh = res(jj).i_lh;
    if mode ==0
        fig1 = figure('Units', 'centimeters', 'Position', [20, 10, 5, 3]);hold on;ylim([0 1])
    end
    
    period_tmp = diff(p)/120*1000;
    res(jj).period = period_tmp;
    period_mean_tmp = mean(period_tmp);
    period(jj) = period_mean_tmp;
    
    t_bin = zeros(length(p)-1,2*length(i_bin));
    t_bin_diff = zeros(length(p)-1,2*length(i_bin)-1);
    phaseDuration = zeros(length(p)-1,3);
    
    i_range = 0;
    i_range_l = 0;
    for ii = 1:length(p)-1
        b = p(ii);e = p(ii+1);
        b_idx = find(res(jj).t == b);e_idx = find(res(jj).t == e);
        t = res(jj).t(b_idx:e_idx) - res(jj).t(b_idx);
        i_single = i_fh(b_idx:e_idx);
        i_single_l = i_lh(b_idx:e_idx);
        range = (max(i_single) - min(i_single));
        i_range = i_range + (max(i_single) - min(i_single))/(length(p)-1);
        range_l = (max(i_single_l) - min(i_single_l));
        i_range_l = i_range_l + (max(i_single_l) - min(i_single_l))/(length(p)-1);
    end
    
    
    i_min_tmp  = zeros(1,length(p)-1);
    for ii = 1:length(p)-1
        b = p(ii);e = p(ii+1);
        b_idx = find(res(jj).t == b);e_idx = find(res(jj).t == e);
        t = res(jj).t(b_idx:e_idx) - res(jj).t(b_idx);
        i_single = i_fh(b_idx:e_idx);
        i_single_l = i_lh(b_idx:e_idx);
        % normalize trace to 0-1
        i_norm = 1-(i_single - min(i_single)) / i_range;
        i_norm_l = 1-(i_single_l - min(i_single_l)) / i_range_l;
        t_norm = t/max(t);
        
        t_inter1 = 0:0.001:max(t);
        i_inter1 = interp1(t,i_norm,t_inter1,'pchip');
        i_inter1_l = interp1(t,i_norm_l,t_inter1,'pchip');
        i_min_tmp(ii) = i_inter1(1);
        xidx = zeros(1,2*length(i_bin));
        for  kk = 1:length(i_bin)
            tmp = find(i_inter1>i_bin(kk));
            xidx(kk) = tmp(1);
            xidx(end-kk+1) = tmp(end);
        end
        t_bin(ii,:) = t_inter1(xidx);
        t_bin_diff(ii,:) = diff(t_inter1(xidx));
        
        if mode == 0
            plot(t_inter1,i_inter1,'color',[0,0.45,0.7]);
        end
        
        thre1 = 0.01;
        thre2 = 0.1;
        thre3 = 0.3;
        i_thre1 = find(i_inter1>thre1);
        i_thre2 = find(i_inter1>thre2);
        i_thre3 = find(i_inter1>thre3);
        idx_thre1_b = i_thre1(1);
        idx_thre1_e = i_thre1(end);
        idx_thre2_b = i_thre2(1);
        idx_thre2_e = i_thre2(end);
        idx_thre3_b = i_thre3(1);
        idx_thre3_e = i_thre3(end);
        
        phaseDuration(ii,:) = [t_inter1(idx_thre2_b) - t_inter1(idx_thre1_b)
            t_inter1(idx_thre3_e) - t_inter1(idx_thre3_b)
            t_inter1(idx_thre1_e) - t_inter1(idx_thre2_e)
            ]*1000;
        
        
        
    end
    for ii = 1:length(p)-1
        b = p(ii);e = p(ii+1);
        b_idx = find(res(jj).t == b);e_idx = find(res(jj).t == e);
        t = res(jj).t(b_idx:e_idx) - res(jj).t(b_idx);
        i_single = i_fh(b_idx:e_idx);
        i_single_l = i_lh(b_idx:e_idx);
        % normalize trace to 0-1
        i_norm = 1-(i_single - min(i_single)) / i_range;
        i_norm_l = 1-(i_single_l - min(i_single_l)) / i_range_l;
        t_norm = t/max(t);
        
        t_inter1 = 0:0.001:max(t);
        i_inter1 = interp1(t,i_norm,t_inter1,'pchip');
        i_inter1_l = interp1(t,i_norm_l,t_inter1,'pchip');
        i_min_tmp(ii) = i_inter1(1);
        xidx = zeros(1,2*length(i_bin));
        for  kk = 1:length(i_bin)
            tmp = find(i_inter1>i_bin(kk));
            xidx(kk) = tmp(1);
            xidx(end-kk+1) = tmp(end);
        end
        t_bin(ii,:) = t_inter1(xidx);
        t_bin_diff(ii,:) = diff(t_inter1(xidx));
        
        if mode == 0
            plot(t_inter1,i_inter1_l,'color',[0.85,0.33,0.1]);
        end
        
        thre1 = 0.01;
        thre2 = 0.1;
        thre3 = 0.3;
        i_thre1 = find(i_inter1>thre1);
        i_thre2 = find(i_inter1>thre2);
        i_thre3 = find(i_inter1>thre3);
        idx_thre1_b = i_thre1(1);
        idx_thre1_e = i_thre1(end);
        idx_thre2_b = i_thre2(1);
        idx_thre2_e = i_thre2(end);
        idx_thre3_b = i_thre3(1);
        idx_thre3_e = i_thre3(end);
        
        phaseDuration(ii,:) = [t_inter1(idx_thre2_b) - t_inter1(idx_thre1_b)
            t_inter1(idx_thre3_e) - t_inter1(idx_thre3_b)
            t_inter1(idx_thre1_e) - t_inter1(idx_thre2_e)
            ]*1000;
        
    end
    i_min(jj) = mean(i_min_tmp);
    t_bin_mean = mean(t_bin);
    t_bin_diff_mean = mean(t_bin_diff);
    t_bin_all(jj,:) = t_bin_mean;
    t_bin_diff_all(jj,:) = t_bin_diff_mean;
    
    phaseDuration_mean = mean(phaseDuration);
    phaseDuration_all(jj,:) = phaseDuration_mean;
    
    t_sum = zeros(1,length(t_bin_diff_mean));
    for ii = 1:length(t_sum)
        t_sum(ii) = sum(t_bin_diff_mean(1:ii));
    end
end

period_all_mean = mean(period);

