function [exclayer,motorlayer,readout] = LIF_SC(options)
arguments
    options.neuron (1,1) double = 1
    options.mu_en (1,1) double = 0.1
    options.mu_mn (1,1) double = 0.1
    options.delta (1,1) double = 2
    options.tau (1,1) double = 5;
    options.mu_ext (1,1) double = 4.8
    options.layer (1,1) double = 10
end
neuron = options.neuron;
mu = options.mu_en;
mu_mn = options.mu_mn;
delta = options.delta;
tau = options.tau;
mu_ext = options.mu_ext;
totlayer = options.layer;


% LIF model for synfire
%% basic paramaters
global VL thre rm tref cm Tm Ve Vi dt fi dvdt
dt = 0.1;
VL = -52; thre = -45; rm = 10; tref = 2.2; cm = 2;Tm = cm*rm;
Ve = 0; Vi = -55;
fi = 15;

dvdt =  @(v,i) (VL-v) / Tm + i /cm ;


tott = 1000;
step = floor(tott/dt);


Iexternal = zeros(1,step);

%% adjustable paramaters
% external inputs
Iexternal(ceil(500/dt):ceil(800/dt)) = mu_ext + 0.05*mu_ext*wgn(ceil(800/dt)-ceil(500/dt)+1,1,2);
kt = 0:dt:40;
kernel = exp(1)*kt/tau^2.*exp(-kt/tau)/neuron/5;

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
    exclayer(l).Isyn = zeros(1,step);
end
for l = 1:totlayer
    motorlayer(l).V = ones(1,step) * -52;
    motorlayer(l).Ibg = zeros(1,step);
    motorlayer(l).Ibg(1,:) = mu_mn + delta*wgn(step,1,2);
    motorlayer(l).isspike = zeros(1,step);
    motorlayer(l).totspike = zeros(1,step);
    motorlayer(l).output = zeros(1,step);
    motorlayer(l).Iext = zeros(1,step);
    motorlayer(l).ref = zeros(1,step);
end

%% experiment
for l = 1:totlayer
    for i = 2:step
        if l>1
            Isyn = exclayer(l-1).output(i)*(Ve-exclayer(l).V(:,i-1));
        else
            Isyn = inputlayer(1).output(i)*(Ve-exclayer(l).V(:,i-1));
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
    
    for i = 2:step
        Isyn = exclayer(l).output(i)*(Ve-motorlayer(l).V(:,i-1));
        Ibg = motorlayer(l).Ibg(:,i);
        Iext = motorlayer(l).Iext(:,i);
        Itot = Ibg + Iext + Isyn;
        [motorlayer(l).V(:,i),motorlayer(l).isspike(:,i)] = updateNeuron(motorlayer(l).V(:,i-1),motorlayer(l).ref(:,i),Itot);
        if sum(motorlayer(l).isspike(:,i))>0
            motorlayer(l).ref(motorlayer(l).isspike(:,i)>0,i:i+ceil(tref/dt)-1) = ones(sum(motorlayer(l).isspike(:,i)),ceil(tref/dt));
        end
    end
    
end
readout = zeros(totlayer,step);
for l = 1:totlayer
    readout(l,:) = motorlayer(l).isspike;
end
sp = zeros(1,2);

%% plot

%% local function
    function [v,ap] = updateNeuron(v,ref,i)
        dv = dvdt(v,i)*dt;
        
        v = v+dv;
        v = v.*(-ref+1) + VL.*ref;
        
        ap = zeros(length(v),1);
        ap(v>thre) = 1;
        v(v>thre) = VL;
        
    end
end