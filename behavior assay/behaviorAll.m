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
    180/255 180/255 180/255];
gene = ["12v";"12v_gal4";"9882_gal4";"11d";"11d_gal4";"10";"10_gal4";"8";"8_gal4";"5";"5_gal4"];

for gg = 1:length(gene)
    load(strcat('.\single fly tnt\data_',gene(gg),'.mat'));
    load(strcat('.\single fly tnt\b_',gene(gg),'.mat'));
    clear res
    for ii = 1:height(array)
        i = array{ii,7};
        i = i(bb(ii)*120:(bb(ii)+3)*120);
        localMin = islocalmin(i,'MinSeparation',10,'FlatSelection', 'first','MinProminence',3);
        localMax = islocalmax(i,'MinSeparation',10,'FlatSelection', 'first','MinProminence',3);
        tf_lmin = find(localMin>0);
        tf_lmax = find(localMax>0);
        i_lmin = i(localMin);
        i_lmax = i(localMax);
        res(ii).i = i;
        res(ii).tf_lmax = tf_lmax;
        res(ii).tf_lmin = tf_lmin;
        res(ii).i_lmin = i_lmin;
        res(ii).i_lmax = i_lmax;
        
        i_fh = array{ii,9};i_lh = array{ii,10};i_fh = i_fh(bb(ii)*120:(bb(ii)+3)*120);i_lh = i_lh(bb(ii)*120:(bb(ii)+3)*120);
        res(ii).i_fh = i_fh;res(ii).i_lh = i_lh;
        localMin_fh = islocalmin(i_fh,'MinSeparation',10,'FlatSelection', 'first','MinProminence',3);
        localMax_fh = islocalmax(i_fh,'MinSeparation',10,'FlatSelection', 'first','MinProminence',3);
        localMin_lh = islocalmin(i_lh,'MinSeparation',10,'FlatSelection', 'first','MinProminence',3);
        localMax_lh = islocalmax(i_lh,'MinSeparation',10,'FlatSelection', 'first','MinProminence',3);
        res(ii).tf_lmax_fh = find(localMax_fh>0);
        res(ii).tf_lmin_fh = find(localMin_fh>0);
        res(ii).tf_lmax_lh = find(localMax_lh>0);
        res(ii).tf_lmin_lh = find(localMin_lh>0);
        
    end

    i_bin = 0.995:-0.01:0.005;
    i_bin = fliplr(i_bin);
    t_bin_all = zeros(length(res),length(i_bin)*2);
    t_bin_diff_all = zeros(length(res),length(i_bin)*2-1);
    period = zeros(1,length(res));
    delay = zeros(1,length(res));
    phaseDuration_all = zeros(length(res),3);
    i_min = zeros(1,length(res));
            peak_all =zeros(1,length(res));
        peak_all_lh =zeros(1,length(res));
    % figure
    for jj = 1:length(res)
        p = res(jj).tf_lmax;
        clear i;
        i = res(jj).i;

        period_tmp = diff(p)/120*1000;
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
            t = 0:1/120:(e-b)/120; i_single = i(b:e);
            range = (max(i_single) - min(i_single));
            i_range = i_range + (max(i_single) - min(i_single))/(length(p)-1);
        end
        
        
        i_min_tmp  = zeros(1,length(p)-1);
        for ii = 1:length(p)-1
            b = p(ii);e = p(ii+1);
            t = 0:1/120:(e-b)/120; i_single = i(b:e);
            
            % normalize trace to 0-1
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
    
    % delay
    for jj = 1:length(res)
        p = res(jj).tf_lmax_fh;
        clear i;
        i = res(jj).i_fh;
        i_lh = res(jj).i_lh;
        i_norm = (i-movmean(i,120))./movmean(i,120);
        i_norm_lh = (i_lh-movmean(i_lh,120))./movmean(i_lh,120);      
        period_tmp = diff(p)/120*1000;
        res(jj).period = period_tmp;
        period_mean_tmp = mean(period_tmp);
        period(jj) = period_mean_tmp;
        
        t_bin = zeros(length(p)-1,2*length(i_bin));
        t_bin_diff = zeros(length(p)-1,2*length(i_bin)-1);
        
        delay_tmp = zeros(1,length(p)-1);
        i_range = 0;
        i_range_lh = 0;
        for ii = 1:length(p)-1
            b = p(ii);e = p(ii+1);
            t = 0:1/120:(e-b)/120; i_single = i(b:e);i_single_lh = i_lh(b:e);
            range = (max(i_single) - min(i_single));
            i_range = i_range + (max(i_single) - min(i_single))/(length(p)-1);
            range_lh = (max(i_single_lh) - min(i_single_lh));
            i_range_lh = i_range_lh + (max(i_single_lh) - min(i_single_lh))/(length(p)-1);
            peak_lh(ii) = max(i_norm_lh(b:e));
            peak(ii) = max(i_norm(b:e));
        end
        peak_all(jj) = mean(peak);
        peak_all_lh(jj) = mean(peak_lh);
        
        
        for ii = 1:length(p)-1
            b = p(ii);e = p(ii+1);
            t = 0:1/120:(e-b)/120; i_single = i(b:e);i_single_lh = i_lh(b:e);
            delay_tmp(ii) = (find(i_single_lh(10:end)==min(i_single_lh(10:end)),1,'first')+9 - find(i_single==min(i_single),1,'first'))/120;
            % normalize trace to 0-1
            i_norm = 1-(i_single - min(i_single)) / i_range;
            i_norm_lh = 1-(i_single_lh - min(i_single_lh)) / i_range_lh;
            t_norm = t/max(t);
            
            t_inter1 = 0:0.001:max(t);
            i_inter1 = interp1(t,i_norm,t_inter1,'pchip');
            i_inter1_lh = interp1(t,i_norm_lh,t_inter1,'pchip');
            
            xidx = zeros(1,2*length(i_bin));
            for  kk = 1:length(i_bin)
                tmp = find(i_inter1>i_bin(kk));
                xidx(kk) = tmp(1);
                xidx(end-kk+1) = tmp(end);
            end
            t_bin(ii,:) = t_inter1(xidx);
            t_bin_diff(ii,:) = diff(t_inter1(xidx));
        end
        delay(jj) = mean(delay_tmp);
    end
    
    period_all_mean = mean(period);
    t_bin_diff_all_mean = mean(t_bin_diff_all);
    t_bin_diff_all_std = std(t_bin_diff_all);
    t_sum_all = zeros(size(t_bin_diff_all));
    for ii = 1:length(t_bin_diff_all)
        t_sum_all(:,ii) = sum(t_bin_diff_all(:,1:ii),2);
    end
    
    % only for i_bin = 0.975:-0.05:0.025;
    p = [sum(t_bin_diff_all_mean(5:15)),sum(t_bin_diff_all_mean(16:22)),sum(t_bin_diff_all_mean(23:33))];
    
    p_mean = mean(phaseDuration_all);
    p_std = std(phaseDuration_all);
    
    c = MyColorHot(1,:);
    x = 1:length(t_bin_diff_all_mean);
    t_sum = zeros(1,length(t_bin_diff_all_mean));
    t_sum_std = zeros(1,length(t_bin_diff_all_mean));
    for ii = 1:length(t_sum)
        t_sum(ii) = sum(t_bin_diff_all_mean(1:ii));
    end
    t_sum_std = std(t_sum_all);
    
    behavior(gg).gene = gene(gg);
    behavior(gg).period = period;
    behavior(gg).interval = phaseDuration_all(:,1)+phaseDuration_all(:,3);
    behavior(gg).fwhm = phaseDuration_all(:,2);
    behavior(gg).start = i_min;
    behavior(gg).delay = delay;
    behavior(gg).peak = peak_all;
    behavior(gg).peak_lh = peak_all_lh;
    
end