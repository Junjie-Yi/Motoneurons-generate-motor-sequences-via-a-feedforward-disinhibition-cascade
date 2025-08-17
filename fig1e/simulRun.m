b = [4.8,6.9,10.7,7.8,5.8,6.2,10,12,10.5,11.5,8.5,6.6,10,7.4,14,14.2];
e = b+3;

thre = [-25,-30,-50,-50,-50,-50,-40,-40,-60,-60,-60,-60,-50,-40,-30,-40];
inputFile = 'simultaneousRecording.xlsx';
dataPath = readtable(inputFile);
for ii = 1:height(dataPath)
    resFile = strcat(dataPath.filePath{ii},'\',dataPath.fileName{ii},'.xlsx');
    ephy = xlsread(resFile,'Sheet1');
    video = xlsread(resFile,'Sheet2');
    t_e = ephy(:,1);
    ep = ephy(:,2);
    t_v = video(:,1);
    i = video(:,2);
    
    b_e = find(t_e>b(ii),1,'first');
    e_e = find(t_e>e(ii),1,'first');
    b_v = find(t_v>b(ii),1,'first');
    e_v = find(t_v>e(ii),1,'first');
    
    ep = ep(b_e:e_e);
    i = -i(b_v:e_v);
    t_e = t_e(b_e:e_e);
    t_v = t_v(b_v:e_v);
    
    lmax_v = find(islocalmax(i,'Minseparation',10,'MinProminence',0.02)>0);
    lmin_e = find(islocalmin(ep,'Minseparation',20,'MinProminence',10)>0);
    lmin_e = lmin_e(ep(lmin_e)<thre(ii));
    tf_spike = t_e(lmin_e);
    spike = zeros(1,1001);
    fr_mean = zeros(1,1001);
    v_mean = zeros(1,1001);
    
    number(ii) = length(lmax_v)-2;
    v_norm_all = zeros(length(lmax_v)-2,1001);
    fr = zeros(length(lmax_v)-2,1001);
    for jj = 1:length(lmax_v)-2
        i_single = i(lmax_v(jj):lmax_v(jj+2));
        t_single = t_v(lmax_v(jj):lmax_v(jj+2));
        t_b = t_v(lmax_v(jj));t_e = t_v(lmax_v(jj+2));
        period = t_e-t_b;
        t_norm = 0:period/1000:period;
        v_norm = interp1(t_single-t_single(1),i_single,t_norm);
        v_norm_all(jj,:) = v_norm;
        v_mean = v_mean + v_norm/(length(lmax_v)-2);
        spikeInPeriod = tf_spike(tf_spike>t_b & tf_spike<t_e)-t_b;
        for kk = 1:length(spikeInPeriod)
            idx = find(t_norm>spikeInPeriod(kk),1,'first');
            spike(idx) = spike(idx)+1;
        end
        if ii<4
            fr(jj,:) = instantaneousRate(t_norm,spikeInPeriod,0.015)';
        else
            fr(jj,:) = instantaneousRate(t_norm,spikeInPeriod,0.008)';
        end
    end
    res(ii).cycle = length(lmax_v)-2;
    res(ii).v_norm = v_norm_all;
    res(ii).fr = fr;
    res(ii).t = 0:1:1000;
end
for ii = 1:5
    mn(ii).cycle = 0;
    mn(ii).v = [];
    mn(ii).fr = [];
    mn(ii).t = 0:1:1000;
end
for ii = 1:2
    mn(1).cycle = mn(1).cycle + res(ii).cycle;
    mn(1).v = [mn(1).v ; res(ii).v_norm];
    mn(1).fr = [mn(1).fr ; res(ii).fr];
end
for ii = 3:6
    mn(2).cycle = mn(2).cycle + res(ii).cycle;
    mn(2).v = [mn(2).v ; res(ii).v_norm];
    mn(2).fr = [mn(2).fr ; res(ii).fr];
end
for ii = 7:10
    mn(3).cycle = mn(3).cycle + res(ii).cycle;
    mn(3).v = [mn(3).v ; res(ii).v_norm];
    mn(3).fr = [mn(3).fr ; res(ii).fr];
end
for ii = 11:14
    mn(4).cycle = mn(4).cycle + res(ii).cycle;
    mn(4).v = [mn(4).v ; res(ii).v_norm];
    mn(4).fr = [mn(4).fr ; res(ii).fr];
end
for ii = 15:16
    mn(5).cycle = mn(5).cycle + res(ii).cycle;
    mn(5).v = [mn(5).v ; res(ii).v_norm];
    mn(5).fr = [mn(5).fr ; res(ii).fr];
end
v_norm_total = zeros(5,1001);
for ii = 1:5
    mn(ii).v_norm = mn(ii).v;
    mn(ii).fr_norm = mn(ii).fr;
    v_mean = mean(mn(ii).v);
    fr_mean = mean(mn(ii).fr);
    for jj = 1:height(mn(ii).v)
        mn(ii).v_norm(jj,:) = (mn(ii).v_norm(jj,:) - min(v_mean))./(max(v_mean)-min(v_mean));
        mn(ii).fr_norm(jj,:) = mn(ii).fr(jj,:)./max(fr_mean);
    end
    v_norm_total(ii,:) = (v_mean-min(v_mean))./(max(v_mean-min(v_mean)));
    lmin1 = find(v_mean==min(v_mean(1:floor(length(v_mean)/2))),1,'first');
    lmin2 = find(v_mean==min(v_mean(ceil(length(v_mean)/2):end)),1,'last');
    
    t_norm = 0:1:1000;
    t_norm = t_norm - t_norm(lmin1);
    t_norm = t_norm * (1/t_norm(lmin2));
    
    v_mean = (v_mean-min(v_mean))./(max(v_mean-min(v_mean)));
    fr_mean = fr_mean/max(fr_mean);
    v_norm_std = std(mn(ii).v_norm);
    fr_std = std(mn(ii).fr_norm);
    v_norm_sem = std(mn(ii).v_norm)/sqrt(mn(ii).cycle);
    fr_sem = std(mn(ii).fr_norm)/sqrt(mn(ii).cycle);
    
    xconf = [t_norm,t_norm(end:-1:1)];
    y1conf = [v_mean+v_norm_std,v_mean(end:-1:1)-v_norm_std(end:-1:1)];
    y2conf = [fr_mean+fr_std,fr_mean(end:-1:1)-fr_std(end:-1:1)];
    %         figure;hold on
    fig = figure('Units', 'centimeters', 'Position', [20, 10, 3, 3.5]);hold on
    set(gca, 'Position', [0.3, 0.1, 0.5, 0.4])
    plot(t_norm,v_mean,'lineWidth',0.5,'color',[0,0.45,0.74])
    fill(xconf,y1conf,[0,0.45,0.74],'linestyle','none','facealpha',0.3)
    plot(t_norm,fr_mean,'lineWidth',0.5,'color','k');
    fill(xconf,y2conf,'k','linestyle','none','facealpha',0.3)
    ax = gca;
    set(ax,'color','none')
    set(ax,'box','off')
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    xlim([-0.6 1.4]);
    ylim([-0.4 1.6]);
    hold off
end


v_norm_total_mean = mean(v_norm_total);
v_norm_total_std = std(v_norm_total);
v_norm_total_sem = std(v_norm_total)/sqrt(length(b));
xconf = [t_norm,t_norm(end:-1:1)];
yconf = [v_norm_total_mean+v_norm_total_std,v_norm_total_mean(end:-1:1)-v_norm_total_std(end:-1:1)];
fig = figure('Units', 'centimeters', 'Position', [20, 10, 3, 3]);hold on
set(gca, 'Position', [0.2, 0.2, 0.5, 0.3])
plot(t_norm,v_norm_total_mean,'lineWidth',0.5);
fill(xconf,yconf,[0,0.45,0.74],'linestyle','none','facealpha',0.3)
ax = gca;
set(ax,'color','none')
set(ax,'box','off')
set(ax,'FontName','Arial');
set(ax,'FontSize',7);
set(ax,'TickDir','out');
xlim([-0.4 1.6]);ylim([-0.1 1.1]);

