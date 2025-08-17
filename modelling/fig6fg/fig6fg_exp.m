parpool('local',28);
neuron = 1; mu_izh = 1;mu_mn_izh = 1;mu_inh_izh = 9; delta_izh = 4;tau = 5;tau_inh = 5;
mu_lif = 0.1;mu_mn_lif = 0.1;mu_inh_lif = 5.5;delta_lif = 2; mu_ext_lif = 4.8;

%%
layer_sc = 10;
layer_fd = 4;


my_ort = zeros(1,50);
my_cas9 = zeros(1,50);
my_mu_inh = ones(1,50)*5.5;

my_ort(1,1:20) = 0.02:0.02:0.4;
my_cas9(1,21:40) = 0.05:0.05:1;
my_mu_inh(1,41:50) = 5:-0.5:0.5;

tottrail = 1000;
for i = 1:50
    res(i).lif_fd_fs = zeros(layer_fd,tottrail);
    res(i).ort = my_ort(i);
    res(i).cas9 = my_cas9(i);
    res(i).mu_inh = my_mu_inh(i);
end




%%
for i = 1:50
    if i <= 50
        lif_fd = zeros(layer_fd,tottrail);
        parfor t = 1:tottrail
            if t<=tottrail
                [~,~,~,st_lif_fd] = LIF_FD_exp('ort',my_ort(i),'cas9',my_cas9(i),'mu_inh',my_mu_inh(i),'layer',4);

                firstSpike_lif = zeros(layer_fd,2);
                for ii = 1:layer_fd
                    if ii <=layer_fd
                        tmp = find(st_lif_fd(ii,:) == 1);
                        if isempty(tmp(tmp>900))
                            firstSpike_lif(ii,2) = -999;
                        else
                            tmp = tmp(tmp>900);
                            firstSpike_lif(ii,2) = tmp(1)*0.1;
                        end
                    end
                end
                lif_fd(:,t) = firstSpike_lif(:,2);
                
            end
            
        end
        res(i).lif_fd_fs = lif_fd;
    end
end
save('manipulation.mat','res');
