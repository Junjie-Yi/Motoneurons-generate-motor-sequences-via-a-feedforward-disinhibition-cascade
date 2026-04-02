function [exclayer,inhlayer,motorlayer,readout] = LIF_FD_exp(options)
arguments
    options.neuron (1,1) double = 1
    options.mu_en (1,1) double = 0.1
    options.mu_mn (1,1) double = 0.1
    options.mu_inh (1,1) double = 5.5
    options.delta (1,1) double = 2
    options.tau (1,1) double = 5;
    options.mu_ext (1,1) double = 4.8
    options.layer (1,1) double = 10
    options.ort (1,1) double = 0
    options.cas9 (1,1) double = 0
end
neuron = options.neuron;
mu = options.mu_en;
mu_mn = options.mu_mn;
mu_inh = options.mu_inh;
delta = options.delta;
tau = options.tau;
tau_inh = tau;
mu_ext = options.mu_ext;
totlayer = options.layer;
ort = options.ort;
cas9 = options.cas9;
% LIF model for feedforward
%% basic paramaters
global VL thre rm tref cm Tm Ve Vi dt fi dvdt VL_inh dvdt_inh
dt = 0.1;
VL = -52; thre = -45; rm = 10; tref = 2.2; cm = 2;Tm = cm*rm;
VL_inh = -52;
Ve = 0; Vi = -55;
fi = 15;

dvdt =  @(v,i) (VL-v) / Tm + i /cm ;
dvdt_inh =  @(v,i) (VL_inh-v) / Tm + i /cm ;

tott = 1000;
step = floor(tott/dt);


Iexternal = zeros(1,step);

%% adjustable paramaters
% external inputs
Iexternal(ceil(500/dt):ceil(800/dt)) = mu_ext + 0.05*mu_ext*wgn(ceil(800/dt)-ceil(500/dt)+1,1,0);
kt = 0:dt:40;
kernel = exp(1)*kt/tau^2.*exp(-kt/tau)/neuron/5;
kernel_inh = exp(1)*kt/tau_inh^2.*exp(-kt/tau_inh)/neuron/5*fi;

%% initiation of inputs
for l = 1:1
    inputlayer(1).V = ones(neuron,step) * -52;
    inputlayer(1).Ibg = zeros(neuron,step);
    for n = 1:neuron
        inputlayer(1).Ibg(n,:) = mu + delta*wgn(step,1,2);
    end
    inputlayer(1).isspike = zeros(neuron,step);
    inputlayer(1).totspike = zeros(1,step);
    inputlayer(1).output = zeros(1,step);
    inputlayer(1).Iext = zeros(neuron,step);
    inputlayer(1).ref = zeros(neuron,step);
    Iexternal(ceil(500/dt):ceil(800/dt)) = mu_ext + 0.05*mu_ext*wgn(ceil(800/dt)-ceil(500/dt)+1,1,2);
    for n = 1:neuron
        inputlayer(1).Iext(n,:) = Iexternal;
    end
end
for i = 2:step
    Isyn = 0;
    Ibg = inputlayer(1).Ibg(:,i);
    Iext = inputlayer(1).Iext(:,i);
    Itot = Ibg + Iext + Isyn;
    [inputlayer(1).V(:,i),inputlayer(1).isspike(:,i)] = updateNeuron(inputlayer(1).V(:,i-1),inputlayer(1).ref(:,i),Itot);
    if sum(inputlayer(1).isspike(:,i))>0
        inputlayer(1).ref(inputlayer(1).isspike(:,i)>0,i:i+ceil(tref/dt)-1) = ones(sum(inputlayer(1).isspike(:,i)),ceil(tref/dt));
    end
    
end
inputlayer(1).totspike = sum(inputlayer(1).isspike);
if neuron > 1
    inputlayer(1).output = conv(inputlayer(1).totspike,kernel);
else
    inputlayer(1).output = conv(inputlayer(1).isspike,kernel);
end
inputlayer(1).output = circshift(inputlayer(1).output(1:step),2/dt);

%% initiation of network
% excitatory layer
for l = 1:totlayer
    exclayer(l).V = ones(neuron,step) * -52;
    exclayer(l).Ibg = zeros(neuron,step);
    for n = 1:neuron
        exclayer(l).Ibg(n,:) = mu + delta*wgn(step,1,2);
    end
    exclayer(l).isspike = zeros(neuron,step);
    exclayer(l).totspike = zeros(1,step);
    exclayer(l).output = zeros(1,step);
    exclayer(l).Iext = zeros(neuron,step);
    exclayer(l).ref = zeros(neuron,step);
end
% inhibitory layer
for l = 1:totlayer
    inhlayer(l).V = ones(neuron,step) * VL_inh;
    inhlayer(l).Ibg = zeros(neuron,step);
    for n = 1:neuron
        inhlayer(l).Ibg(n,:) = mu_inh + delta*wgn(step,1,2);
    end
    inhlayer(l).isspike = zeros(neuron,step);
    inhlayer(l).totspike = zeros(1,step);
    inhlayer(l).output = zeros(1,step);
    inhlayer(l).Iext = zeros(neuron,step);
    inhlayer(l).ref = zeros(neuron,step);
end
% motor layer
for l = 1:totlayer
    motorlayer(l).V = ones(1,step) * -52;
    motorlayer(l).Ibg = mu_mn + delta*wgn(step,1,2)';
    motorlayer(l).isspike = zeros(1,step);
    motorlayer(l).output = zeros(1,step);
    motorlayer(l).Iext = zeros(1,step);
    motorlayer(l).ref = zeros(1,step);
end
%% experiment
for l = 1:totlayer
    
    % update exitatory layer
    for i = 2:step
        Isyn = inputlayer(1).output(i)*(Ve-exclayer(1).V(:,i-1));
        if l>1
            Isyn = Isyn + inhlayer(l-1).output(i)*(Vi-exclayer(l).V(:,i-1));
        end
        Ibg = exclayer(l).Ibg(:,i);
        Iext = exclayer(l).Iext(:,i);
        Itot = Ibg + Iext + Isyn;
        [exclayer(l).V(:,i),exclayer(l).isspike(:,i)] = updateNeuron(exclayer(l).V(:,i-1),exclayer(l).ref(:,i),Itot);
        if sum(exclayer(l).isspike(:,i))>0
            exclayer(l).ref(exclayer(l).isspike(:,i)>0,i:i+ceil(tref/dt)-1) = ones(sum(exclayer(l).isspike(:,i)),ceil(tref/dt));
        end
        
    end
    exclayer(l).totspike = sum(exclayer(l).isspike);
    if neuron > 1
        exclayer(l).output = conv(exclayer(l).totspike,kernel);
    else
        exclayer(l).output = conv(exclayer(l).isspike,kernel);
    end
    exclayer(l).output = circshift(exclayer(l).output(1:step),2/dt);
    
    % update MN layer
    for i = 2:step
        Isyn = exclayer(l).output(i)*(Ve-motorlayer(l).V(:,i-1));
        Ibg = motorlayer(l).Ibg(:,i);
        Iext = motorlayer(l).Iext(:,i);
        Itot = Ibg + Iext + Isyn + ort.*(Vi-motorlayer(l).V(:,i-1));
        [motorlayer(l).V(:,i),motorlayer(l).isspike(:,i)] = updateNeuron(motorlayer(l).V(:,i-1),motorlayer(l).ref(:,i),Itot);
        if sum(motorlayer(l).isspike(:,i))>0
            motorlayer(l).ref(motorlayer(l).isspike(:,i)>0,i:i+ceil(tref/dt)-1) = ones(sum(motorlayer(l).isspike(:,i)),ceil(tref/dt));
        end
    end
    motorlayer(l).output = conv(motorlayer(l).isspike,kernel_inh);
    motorlayer(l).output = circshift(motorlayer(l).output(1:step),2/dt)*neuron;
    
    % update inhibitory layer
    for i = 2:step
        Isyn = (1-cas9)*motorlayer(l).output(i)*(Vi-inhlayer(l).V(:,i-1));
        Ibg = inhlayer(l).Ibg(:,i);
        Iext = inhlayer(l).Iext(:,i);
        Itot = Ibg + Iext + Isyn;
        [inhlayer(l).V(:,i),inhlayer(l).isspike(:,i)] = updateNeuron_inh(inhlayer(l).V(:,i-1),inhlayer(l).ref(:,i),Itot);
        if sum(inhlayer(l).isspike(:,i))>0
            inhlayer(l).ref(inhlayer(l).isspike(:,i)>0,i:i+ceil(tref/dt)-1) = ones(sum(inhlayer(l).isspike(:,i)),ceil(tref/dt));
        end
    end
    inhlayer(l).totspike = sum(inhlayer(l).isspike);
    if neuron > 1
        inhlayer(l).output = conv(inhlayer(l).totspike,kernel_inh);
    else
        inhlayer(l).output = conv(inhlayer(l).isspike,kernel_inh);
    end
    inhlayer(l).output = circshift(inhlayer(l).output(1:step),2/dt);
end
readout = zeros(totlayer,step);
for l = 1:totlayer
    readout(l,:) = motorlayer(l).isspike;
end
%% plot
% figure;hold on
% 
% for l = 1:5
%     plot(dt:dt:tott,motorlayer(l).V(1,:))
% end

%% local function
    function [v,ap] = updateNeuron(v,ref,i)
        dv = dvdt(v,i)*dt;
        
        v = v+dv;
        v = v.*(-ref+1) + VL.*ref;
        
        ap = zeros(length(v),1);
        ap(v>thre) = 1;
        v(v>thre) = VL;
        
    end
    function [v,ap] = updateNeuron_inh(v,ref,i)
        dv = dvdt_inh(v,i)*dt;
        
        v = v+dv;
        v = v.*(-ref+1) + VL_inh.*ref;
        
        ap = zeros(length(v),1);
        ap(v>thre) = 1;
        v(v>thre) = VL_inh;
        
    end
end