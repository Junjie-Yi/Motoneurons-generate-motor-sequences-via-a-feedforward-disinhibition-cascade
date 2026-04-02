% neuron-delay
clear
load('s7a_bottom.mat');
% LIF model 
for ii = 1:8
    neuron(ii) = 2^(ii-1);
end
y_sc = zeros(1,8);
y_fd = zeros(1,8);

for ii = 1:8
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
figure;hold on
plot(0:7,(y_sc));
plot(0:7,(y_fd));