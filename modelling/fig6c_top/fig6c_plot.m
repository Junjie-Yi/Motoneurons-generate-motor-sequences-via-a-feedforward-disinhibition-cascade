clear
y_lif = zeros(2,9);
x = zeros(2,9);

load('tau_LIF')
k = 2;
for ii = 1:9
    tmp = res(ii).lif_sc_fs(k*3-2,:) - res(ii).lif_sc_fs(k*3-5,:);
    % excluding extreme trails
    x(1,ii) = mean(tmp(tmp>-20&tmp<150));
    y_lif(1,ii) = var(tmp(tmp>-20&tmp<150));
    tmp = res(ii).lif_fd_fs(k,:) - res(ii).lif_fd_fs(k-1,:);
    x(2,ii) = mean(tmp(tmp>-20&tmp<150));
    y_lif(2,ii) = var(tmp(tmp>-20&tmp<150));
end

c = [0,0.45,0.74;
    0.85,0.325,0.098];
figure;hold on
for jj = 1:2
    plot(x(jj,:),y_lif(jj,:),'-o');
    y_up(jj,:) = y_lif(jj,:)*999/chi2inv(0.025,1000);
    y_down(jj,:) = y_lif(jj,:)*999/chi2inv(0.975,1000);

 fill([x(jj,:), fliplr(x(jj,:))], ...
    [ y_lif(jj,:)*999/chi2inv(0.025,1000),fliplr(y_lif(jj,:)*999/chi2inv(0.975,1000))], ...
     c(jj,:),'linestyle','none','FaceAlpha',0.3)
end

