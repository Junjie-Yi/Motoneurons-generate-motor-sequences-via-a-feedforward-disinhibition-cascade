% hh model 
clear
load('tau_HH.mat');
for ii = 1:10
    tau(ii) = res(ii).tau;
end
y_sc = zeros(1,10);
y_fd = zeros(1,10);

for ii = 1:10
    tmp = res(ii).hh_sc_fs(4,:)-res(ii).hh_sc_fs(1,:);
    tmp = tmp(tmp>-20&tmp<10*tau(ii));
    y_sc(ii) = var(tmp);
    x_sc(ii) = mean(tmp);
    bottom_sc(ii) = y_sc(ii) - (length(tmp)-1) * y_sc(ii)/(chi2inv(0.975,length(tmp)-1));
    upper_sc(ii) = (length(tmp)-1) * y_sc(ii)/(chi2inv(0.025,length(tmp)-1)) - y_sc(ii);
    
    tmp = res(ii).hh_fd_fs(2,:)-res(ii).hh_fd_fs(1,:);
    tmp = tmp(tmp>-20&tmp<10*tau(ii));
    y_fd(ii) = var(tmp);
    x_fd(ii) = mean(tmp);
    bottom_fd(ii) = y_fd(ii)-(length(tmp)-1) * y_fd(ii)/(chi2inv(0.975,length(tmp)-1));
    upper_fd(ii) = (length(tmp)-1) * y_fd(ii)/(chi2inv(0.025,length(tmp)-1))-y_fd(ii);
end
% 
c = [0,0.45,0.74;
    0.85,0.325,0.098];
figure;hold on
plot(x_sc,(y_sc));
 fill([x_sc, fliplr(x_sc)], ...
    [ y_sc+upper_sc,fliplr(y_sc-bottom_sc)], ...
     c(1,:),'linestyle','none','FaceAlpha',0.3)
plot(x_fd,(y_fd));
 fill([x_fd, fliplr(x_fd)], ...
    [ y_fd+upper_fd,fliplr(y_fd-bottom_fd)], ...
     c(2,:),'linestyle','none','FaceAlpha',0.3)