
ii = 16;
seq = randperm(100);


fig1 = figure('visible','off');hold on
set(gcf,'unit','centimeters','position',[10,5,4,1.5]);
ax = gca;
set(ax,'FontName','Arial','Fontsize',7);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ylim([0 340]);
xlim([-100 100]);
filename = 'sc4.tif';
for k = 1:100
    t = seq(k);
    st = res(ii).isspike_sc(1,t).st-5000;
    for jj = 1:length(st)
        line([ceil(st(jj)/10) ceil(st(jj)/10)],[t-1 t],'color','k','LineWidth',0.5)
    end
    st = res(ii).isspike_sc(4,t).st-5000;
    for jj = 1:length(st)
        line([ceil(st(jj)/10) ceil(st(jj)/10)],[t-1+120 t+120],'color','k','LineWidth',0.5)
    end
    st = res(ii).isspike_sc(7,t).st-5000;
    for jj = 1:length(st)
        line([ceil(st(jj)/10) ceil(st(jj)/10)],[t-1+240 t+240],'color','k','LineWidth',0.5)
    end
end
exportgraphics(fig1,filename,'Resolution',1000)
close(fig1)


fig1 = figure('visible','off');hold on
set(gcf,'unit','centimeters','position',[10,5,4,1.5]);
ax = gca;
set(ax,'FontName','Arial','Fontsize',7);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ylim([0 340]);
xlim([-100 100]);
filename = 'fd4.tif';
for k = 1:100
    t = seq(k);
    st = res(ii).isspike_fd(1,t).st-5000;
    for jj = 1:length(st)
        line([ceil(st(jj)/10) ceil(st(jj)/10)],[t-1 t],'color','k','LineWidth',0.5)
    end
    st = res(ii).isspike_fd(2,t).st-5000;
    for jj = 1:length(st)
        line([ceil(st(jj)/10) ceil(st(jj)/10)],[t-1+120 t+120],'color','k','LineWidth',0.5)
    end
    st = res(ii).isspike_fd(3,t).st-5000;
    for jj = 1:length(st)
        line([ceil(st(jj)/10) ceil(st(jj)/10)],[t-1+240 t+240],'color','k','LineWidth',0.5)
    end
end
exportgraphics(fig1,filename,'Resolution',1000)
close(fig1)
