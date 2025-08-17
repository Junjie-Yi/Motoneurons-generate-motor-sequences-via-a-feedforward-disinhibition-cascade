parpool('local',28);
neuron = 1; mu_izh = 1;mu_mn_izh = 1;mu_inh_izh = 9; delta_izh = 4;tau = 5;tau_inh = 5;
mu_lif = 0.1;mu_mn_lif = 0.1;mu_inh_lif = 5.5;delta_lif = 2; mu_ext_lif = 4.8;

%%
layer_sc = 10;
layer_fd = 4;
for i = 1:50
    my_mu_en(i) = 0.1 + 0.02 * (i-1);
end


for i = 1:50
        tottrail = 1000;
        res(i).mu_en = my_mu_en(i);
        res(i).lif_fd_fs = zeros(layer_fd,tottrail);
        res(i).lif_sc_fs = zeros(layer_sc,tottrail);
        res(i).spfr_fd = zeros(layer_fd,tottrail);
        res(i).spfr_sc = zeros(layer_sc,tottrail);     
end

for i = 1:50
    if i <= 50

        lif_sc = zeros(layer_sc,tottrail);
        lif_fd = zeros(layer_fd,tottrail);
        hh_sc = zeros(layer_sc,tottrail);
        hh_fd = zeros(layer_fd,tottrail);
        sp_sc = zeros(layer_sc,tottrail);
        sp_fd = zeros(layer_fd,tottrail);
        
        parfor t = 1:tottrail
            if t<=tottrail
                [sc_exc,~,st_lif_sc] = LIF_SC('layer',10,'mu_en',my_mu_en(i));
                [fd_exc,~,~,st_lif_fd] = LIF_FD('layer',4,'mu_en',my_mu_en(i));
                          
                for ii = 1:layer_sc
                   sp_sc(ii,t) = sum(sum(sc_exc(ii).isspike(:,100/0.1:400/0.1)));
                end
                for ii = 1:layer_fd
                   sp_fd(ii,t) = sum(sum(fd_exc(ii).isspike(:,100/0.1:400/0.1)));
                end
                
                firstSpike_lif = zeros(layer_sc,2);
                for ii = 1:layer_sc
                    if ii <=layer_sc
                        tmp = find(st_lif_sc(ii,:) == 1);
                        isspike_sc(ii,t).st = tmp;
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
                        isspike_fd(ii,t).st = tmp;
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
        
        res(i).spfr_sc = sp_sc;
        res(i).spfr_fd = sp_fd;
        
        res(i).isspike_sc = isspike_sc;
        res(i).isspike_fd = isspike_fd;
    end
end
save('fig6d_data.mat','res');