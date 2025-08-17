% delta-delay
load('s7a_top.mat');
% LIF model 
delta = 0.1:0.1:3;
y_sc = zeros(1,30);
y_fd = zeros(1,30);

for ii = 1:30
    tmp = res(ii).lif_sc_fs(4,:)-res(ii).lif_sc_fs(1,:);
    tmp = tmp(tmp>-20&tmp<50);
    y_sc(ii) = var(tmp);
    bottom_sc(ii) = y_sc(ii) - (length(tmp)-1) * y_sc(ii)/(chi2inv(0.975,length(tmp)-1));
    upper_sc(ii) = (length(tmp)-1) * y_sc(ii)/(chi2inv(0.025,length(tmp)-1)) - y_sc(ii);
    
    tmp = res(ii).lif_fd_fs(2,:)-res(ii).lif_fd_fs(1,:);
    tmp = tmp(tmp>-20&tmp<50);
    y_fd(ii) = var(tmp);
    bottom_fd(ii) = y_fd(ii)-(length(tmp)-1) * y_fd(ii)/(chi2inv(0.975,length(tmp)-1));
    upper_fd(ii) = (length(tmp)-1) * y_fd(ii)/(chi2inv(0.025,length(tmp)-1))-y_fd(ii);
end
% 
figure;hold on;ylim([0 100])
plot(delta,(y_sc));
plot(delta,(y_fd));