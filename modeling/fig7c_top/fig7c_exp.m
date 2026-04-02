
neuron = 1; mu_izh = 1;mu_mn_izh = 1;mu_inh_izh = 9; delta_izh = 4;tau = 5;tau_inh = 5;
mu_lif = 0.1;mu_mn_lif = 0.1;mu_inh_lif = 5.5;delta_lif = 2; mu_ext_lif = 4.8;

%%
layer_sc = 10;
layer_fd = 4;
for i = 1:20
    my_tau(i) = 2*i;
    my_neuron(i) = 1;
    if i > 10
        my_tau(i) = 5;
        my_neuron(i) = 10*(i-11);
    end
end
my_neuron(11) = 1;

my_ext_lif = [repmat(4.8,1,20),repmat(4.2,1,20),repmat(3.6,1,20),repmat(3.0,1,20),repmat(2.4,1,20)];
my_ext_hh = [repmat(12,1,20),repmat(10,1,20),repmat(8,1,20),repmat(6,1,20),repmat(4,1,20)];
my_tau = repmat(my_tau,1,5);
my_neuron = repmat(my_neuron,1,5);

for i = 1:100
        tottrail = 1000;
        res(i).tau = my_tau(i);
        res(i).neuron = my_neuron(i);
        res(i).ext_lif = my_ext_lif(i);
        res(i).ext_hh = my_ext_hh(i);
        res(i).lif_fd_fs = zeros(layer_fd,tottrail);
        res(i).lif_sc_fs = zeros(layer_sc,tottrail);
        res(i).hh_fd_fs = zeros(layer_fd,tottrail);
        res(i).hh_sc_fs = zeros(layer_sc,tottrail);
end

%%
for i = 1:100
    if i <= 100
        tau = my_tau(i);
        tau_inh = tau;

        lif_sc = zeros(layer_sc,tottrail);
        lif_fd = zeros(layer_fd,tottrail);
        hh_sc = zeros(layer_sc,tottrail);
        hh_fd = zeros(layer_fd,tottrail);
        
        parfor t = 1:tottrail
            if t<=tottrail
                [~,~,st_lif_sc] = LIF_SC('neuron',my_neuron(i),'mu_ext',my_ext_lif(i),'tau',tau,'layer',10);
                [~,~,~,st_lif_fd] = LIF_FD('neuron',my_neuron(i),'mu_ext',my_ext_lif(i),'tau',tau,'layer',4);

                firstSpike_lif = zeros(layer_sc,2);
                for ii = 1:layer_sc
                    if ii <=layer_sc
                        tmp = find(st_lif_sc(ii,:) == 1);
%                         res(i).lif_sc(t).st(ii).st =tmp;
                        if isempty(tmp(tmp>900))
                            firstSpike_lif(ii,1) = -999;
                        else
                            tmp = tmp(tmp>900);
                            firstSpike_lif(ii,1) = tmp(1)*0.1;
                        end
                    end
                end
                lif_sc(:,t) = firstSpike_lif(:,1);
                 firstSpike_lif = zeros(layer_fd,2);
                for ii = 1:layer_fd
                    if ii <=layer_fd
                        tmp = find(st_lif_fd(ii,:) == 1);
%                         res(i).lif_fd(t).st(ii).st =tmp;
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
        res(i).lif_sc_fs = lif_sc;
    end
end
save('20240903_tau.mat','res');
clear res

%%
layer_sc = 10;
layer_fd = 4;
for i = 1:21
    for j = 1:21
    my_mu_en(i*21-21+j) = 1 + 0.02 * (i-1);
    my_mu_mn(i*21-21+j) = 1 + 0.02 * (j-1);
    end
end


for i = 1:21*21
        tottrail = 1000;
        res(i).mu_en = my_mu_en(i);
        res(i).mu_mn = my_mu_mn(i);
        res(i).lif_fd_fs = zeros(layer_fd,tottrail);
        res(i).lif_sc_fs = zeros(layer_sc,tottrail);
        res(i).spfr_fd = zeros(layer_fd,tottrail);
        res(i).spfr_sc = zeros(layer_sc,tottrail);
end

for i = 1:21*21
    if i <= 21*21

        lif_sc = zeros(layer_sc,tottrail);
        lif_fd = zeros(layer_fd,tottrail);
        hh_sc = zeros(layer_sc,tottrail);
        hh_fd = zeros(layer_fd,tottrail);
        sp_sc = zeros(layer_sc,tottrail);
        sp_fd = zeros(layer_fd,tottrail);
        
        
        parfor t = 1:tottrail
            if t<=tottrail
%                 [res(i).exc_lif_sc,res(i).motor_lif_sc,st_lif_sc] = LIF_SC('neuron',my_neuron(i),'mu_ext',my_ext_lif(i),'tau',tau);
%                 [res(i).exc_lif_fd,res(i).motor_lif_fd,res(i).inh_lif_fd,st_lif_fd] = LIF_FD('neuron',my_neuron(i),'mu_ext',my_ext_lif(i),'tau',tau);
%                 [res(i).exc_hh_sc,res(i).motor_hh_sc,st_hh_sc] = HH_SC('neuron',my_neuron(i),'mu_ext',my_ext_hh(i),'tau',tau,'layer',5);
%                 [res(i).exc_hh_fd,res(i).motor_hh_fd,res(i).inh_hh_fd,st_hh_fd] = HH_FD('neuron',my_neuron(i),'mu_ext',my_ext_hh(i),'tau',tau,'layer',5);
                [sc_exc,~,st_lif_sc] = LIF_SC('layer',10,'mu_en',my_mu_en(i),'mu_mn',my_mu_mn(i));
                [fd_exc,~,~,st_lif_fd] = LIF_FD('layer',4,'mu_en',my_mu_en(i),'mu_mn',my_mu_mn(i));
%                 [~,~,st_hh_sc] = HH_SC('neuron',my_neuron(i),'mu_ext',my_ext_hh(i),'tau',tau,'layer',10);
%                 [~,~,~,st_hh_fd] = HH_FD('neuron',my_neuron(i),'mu_ext',my_ext_hh(i),'tau',tau,'layer',4);
                
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
%                         res(i).lif_sc(t).st(ii).st =tmp;
                        if isempty(tmp(tmp>900))
                            firstSpike_lif(ii,1) = -999;
                        else
                            tmp = tmp(tmp>900);
                            firstSpike_lif(ii,1) = tmp(1)*0.1;
                        end
                    end
                end
                lif_sc(:,t) = firstSpike_lif(:,1);
                firstSpike_lif = zeros(layer_fd,2);
                for ii = 1:layer_fd
                    if ii <=layer_fd
                        tmp = find(st_lif_fd(ii,:) == 1);
%                         res(i).lif_fd(t).st(ii).st =tmp;
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
        res(i).lif_sc_fs = lif_sc;
        res(i).spfr_sc = sp_sc;
        res(i).spfr_fd = sp_fd;
    end
end
save('tau_LIF.mat','res');