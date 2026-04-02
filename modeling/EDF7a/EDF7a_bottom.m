parpool('local',10);
%%
layer_sc = 10;
layer_fd = 4;
tottrail = 1000;

for i = 1:8
        res(i).neuron = 2^(i-1);
        myneuron(i) = 2^(i-1);
        res(i).lif_fd_fs = zeros(layer_fd,tottrail);
        res(i).lif_sc_fs = zeros(layer_sc,tottrail);
        res(i).hh_fd_fs = zeros(layer_fd,tottrail);
        res(i).hh_sc_fs = zeros(layer_sc,tottrail);
end
for i = 1:8
    if i <= 8

        lif_sc = zeros(layer_sc,tottrail);
        lif_fd = zeros(layer_fd,tottrail);
        hh_sc = zeros(layer_sc,tottrail);
        hh_fd = zeros(layer_fd,tottrail);
        
        
        parfor t = 1:tottrail
            if t<=tottrail
                [~,~,st_lif_sc] = LIF_SC('layer',10,'delta',2,'neuron',myneuron(i));
                [~,~,~,st_lif_fd] = LIF_FD('layer',4,'delta',2,'neuron',myneuron(i));
                [~,~,st_hh_sc] = HH_SC('layer',10,'delta',10,'neuron',myneuron(i));
                [~,~,~,st_hh_fd] = HH_FD('layer',4,'delta',10,'neuron',myneuron(i));
                
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
                
                % HH
                firstSpike_hh = zeros(layer_sc,2);
                for ii = 1:layer_sc
                    if ii <=layer_sc
                        tmp = find(st_hh_sc(ii,:) == 1);
                        isspike_sc(ii,t).st = tmp;
                        if isempty(tmp(tmp>9000))
                            firstSpike_hh(ii,1) = -9999;
                        else
                            tmp = tmp(tmp>9000);
                            firstSpike_hh(ii,1) = tmp(1)*0.01;
                        end
                    end
                end
                hh_sc(:,t) = firstSpike_hh(:,1);
                firstSpike_hh = zeros(layer_fd,2);
                for ii = 1:layer_fd
                    if ii <=layer_fd
                        tmp = find(st_hh_fd(ii,:) == 1);
                        isspike_fd(ii,t).st = tmp;
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
        res(i).lif_fd_fs = lif_fd;
        res(i).lif_sc_fs = lif_sc;
        res(i).hh_fd_fs = hh_fd;
        res(i).hh_sc_fs = hh_sc;
    end
end
save('./s7a_bottom.mat','res');