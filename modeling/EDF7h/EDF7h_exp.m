parpool('local',28);
%%
layer_sc = 10;
layer_fd = 4;


my_ort = zeros(1,61);
my_cas9 = zeros(1,61);
my_mu_inh = ones(1,61)*12;

my_ort(1,1:20) = 0.02:0.02:0.4;
my_cas9(1,21:40) = 0.05:0.05:1;
my_mu_inh(1,41:61) = 12:-0.6:0;

tottrail = 1000;
for i = 1:61
    res(i).hh_fd_fs = zeros(layer_fd,tottrail);
    res(i).ort = my_ort(i);
    res(i).cas9 = my_cas9(i);
    res(i).mu_inh = my_mu_inh(i);
end




%%
for i = 1:61
    if i <= 61
        hh_fd = zeros(layer_fd,tottrail);
        parfor t = 1:tottrail
            if t<=tottrail
                [~,~,~,st_hh_fd] = HH_FD_exp('ort',my_ort(i),'cas9',my_cas9(i),'mu_inh',my_mu_inh(i),'layer',4);

                firstSpike_hh = zeros(layer_fd,2);
                for ii = 1:layer_fd
                    if ii <=layer_fd
                        tmp = find(st_hh_fd(ii,:) == 1);
                        if isempty(tmp(tmp>9000))
                            firstSpike_hh(ii,2) = -9999;
                        else
                            tmp = tmp(tmp>9000);
                            firstSpike_hh(ii,2) = tmp(1)*0.01;
                        end
                    end
                end
                hh_fd(:,t) = firstSpike_hh(:,2);
                
            end
            
        end
        res(i).hh_fd_fs = hh_fd;
    end
end
save('manipulation.mat','res');
