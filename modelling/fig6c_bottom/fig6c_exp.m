
neuron = 1; mu_izh = 1;mu_mn_izh = 1;mu_inh_izh = 9; delta_izh = 4;tau = 5;tau_inh = 5;
mu_lif = 0.1;mu_mn_lif = 0.1;mu_inh_lif = 5.5;delta_lif = 2; mu_ext_lif = 4.8;

%%
mp_sc = zeros(1000,1000);
mp_fd = zeros(1000,1000);



layer_sc = 10;
layer_fd = 4;


for i = 1:1
        tottrail = 1000;
        res(i).lif_fd_fs = zeros(layer_fd,tottrail);
        res(i).lif_sc_fs = zeros(layer_sc,tottrail);
end

for i = 1:1
    if i <= 1

        lif_sc = zeros(layer_sc,tottrail);
        lif_fd = zeros(layer_fd,tottrail);
        
        parfor t = 1:tottrail
            if t<=tottrail
                [sc_exc,~,st_lif_sc] = LIF_SC('layer',10);
                [fd_exc,~,~,st_lif_fd] = LIF_FD('layer',4);
                mp_sc(:,t) = sc_exc(4).V(4501:5500);
                mp_fd(:,t) = fd_exc(2).V(4501:5500);
                
                firstSpike_lif = zeros(layer_sc,2);
                for ii = 1:layer_sc
                    if ii <=layer_sc
                        tmp = find(st_lif_sc(ii,:) == 1);
                        if isempty(tmp(tmp>4000))
                            firstSpike_lif(ii,1) = -9999;
                        else
                            tmp = tmp(tmp>4000);
                            firstSpike_lif(ii,1) = tmp(1)*0.1;
                        end
                    end
                end
                lif_sc(:,t) = firstSpike_lif(:,1);
                firstSpike_lif = zeros(layer_fd,2);
                for ii = 1:layer_fd
                    if ii <=layer_fd
                        tmp = find(st_lif_fd(ii,:) == 1);
                        if isempty(tmp(tmp>4000))
                            firstSpike_lif(ii,2) = -9999;
                        else
                            tmp = tmp(tmp>4000);
                            firstSpike_lif(ii,2) = tmp(1)*0.1;
                        end
                    end
                end
                lif_fd(:,t) = firstSpike_lif(:,2);
                

            end
            
        end
        res(i).lif_fd_fs = lif_fd;
        res(i).lif_sc_fs = lif_sc;
        
    end
end
save('mp_LIF.m','res','mp_sc','mp_fd');