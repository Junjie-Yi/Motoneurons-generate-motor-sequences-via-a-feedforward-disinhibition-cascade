
%%
tottrail = 160;
mp_sc = zeros(10000,tottrail);
mp_fd = zeros(10000,tottrail);



layer_sc = 10;
layer_fd = 4;


for i = 1:1
        res(i).hh_fd_fs = zeros(layer_fd,tottrail);
        res(i).hh_sc_fs = zeros(layer_sc,tottrail);
end

for i = 1:1
    if i <= 1

        hh_sc = zeros(layer_sc,tottrail);
        hh_fd = zeros(layer_fd,tottrail);
        
        for t = 1:tottrail
            if t<=tottrail
                [sc_exc,~,st_hh_sc] = HH_SC('layer',10,'neuron',1,'mu_ext',15,'delta',10);
                [fd_exc,~,~,st_hh_fd] = HH_FD('layer',4,'neuron',1,'mu_ext',15,'delta',10);
                mp_sc(:,t) = sc_exc(4).v(1,5001:15000)';
                mp_fd(:,t) = fd_exc(2).v(1,5001:15000)';
                
                firstSpike_hh = zeros(layer_sc,2);
                for ii = 1:layer_sc
                    if ii <=layer_sc
                        tmp = find(st_hh_sc(ii,:) == 1);
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
        res(i).hh_sc_fs = hh_sc;
        
    end
end
save('mp_HH.mat','res','mp_sc','mp_fd');