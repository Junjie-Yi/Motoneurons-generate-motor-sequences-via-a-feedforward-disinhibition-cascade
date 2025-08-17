load('fig6d_data.mat');
falseFiring_sc = zeros(50,3);
falseFiring_fd = zeros(50,3);
spontanous_sc = zeros(50,10);
spontanous_fd = zeros(50,4);
for ii = 1:50
    mu_en(ii) = res(ii).mu_en;
    for l = 1:3
        falseFiring_sc(ii,l) = sum(res(ii).lif_sc_fs(l*3+1,:)-res(ii).lif_sc_fs(l*3-2,:)<0)/1000;
    end
    for l = 1:3
        falseFiring_fd(ii,l) = sum(res(ii).lif_fd_fs(l+1,:)-res(ii).lif_fd_fs(l,:)<0)/1000;
    end    
    for l = 1:10
        spontanous_sc(ii,l) = mean(res(ii).spfr_sc(l,:))/0.3;
    end
    for l =1:4
        spontanous_fd(ii,l) = mean(res(ii).spfr_fd(l,:))/0.3;
    end
end
c = [0,0.45,0.74;
    0.85,0.325,0.098];
figure; hold on
plot(mu_en,falseFiring_sc(:,1),'color',c(1,:));
plot(mu_en,falseFiring_fd(:,1),'color',c(2,:));
figure;
imagesc(spontanous_sc(:,end:-1:1)', [0 300])
colorbar
figure;
imagesc(spontanous_fd(:,end:-1:1)', [0 300])