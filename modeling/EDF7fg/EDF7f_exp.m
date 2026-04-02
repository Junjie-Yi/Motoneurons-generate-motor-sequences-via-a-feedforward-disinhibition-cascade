parpool('local',28);


%%
layer_sc = 10;
layer_fd = 4;
for i = 1:50
    my_mu_en(i) = 0.2 + 0.04 * (i-1);
end

 tottrail = 1000;
for i = 1:50
 
        res(i).mu_en = my_mu_en(i);
        res(i).hh_fd_fs = zeros(layer_fd,tottrail);
        res(i).hh_sc_fs = zeros(layer_sc,tottrail);
        res(i).spfr_fd = zeros(layer_fd,tottrail);
        res(i).spfr_sc = zeros(layer_sc,tottrail);     
end

for i = 1:50
    if i <= 50

        hh_sc = zeros(layer_sc,tottrail);
        hh_fd = zeros(layer_fd,tottrail);
        hh_sc = zeros(layer_sc,tottrail);
        hh_fd = zeros(layer_fd,tottrail);
        sp_sc = zeros(layer_sc,tottrail);
        sp_fd = zeros(layer_fd,tottrail);
        
        parfor t = 1:tottrail
            if t<=tottrail
                [sc_exc,~,st_hh_sc] = HH_SC('layer',10,'mu_en',my_mu_en(i));
                [fd_exc,~,~,st_hh_fd] = HH_FD('layer',4,'mu_en',my_mu_en(i));
                          
                for ii = 1:layer_sc
                   sp_sc(ii,t) = sum(sum(sc_exc(ii).se(:,50/0.01:100/0.01)));
                end
                for ii = 1:layer_fd
                   sp_fd(ii,t) = sum(sum(fd_exc(ii).se(:,50/0.01:100/0.01)));
                end
                
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
        res(i).hh_fd_fs = hh_fd;
        res(i).hh_sc_fs = hh_sc;
        
        res(i).spfr_sc = sp_sc;
        res(i).spfr_fd = sp_fd;
        
        res(i).isspike_sc = isspike_sc;
        res(i).isspike_fd = isspike_fd;
    end
end
save('./figs7f_data.mat','res');