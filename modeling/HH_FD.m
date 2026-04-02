function [exclayer,motorlayer,inhlayer,readout] = HH_FD(options)
arguments
    options.neuron (1,1) double = 1
    options.mu_en (1,1) double = 0.5
    options.mu_mn (1,1) double = 0.5
    options.mu_inh (1,1) double = 12
    options.delta (1,1) double = 10
    options.tau (1,1) double = 5;
    options.mu_ext (1,1) double = 15
    options.layer (1,1) double = 4
end
neuron = options.neuron;
mu_en = options.mu_en;
mu_mn = options.mu_mn;
mu_inh = options.mu_inh;
delta = options.delta;
tau = options.tau;
mu_ext = options.mu_ext;
totlayer = options.layer;
% feedforward disinhibition with HH model
% HH model in Moldakarimov et al. 2015 PNAS

%% parameter
gL = 0.05; gNa = 100; gK = 40;dt = 0.01;
VL = -65; VNa = 55; VK = -80; Ve = 0;Vi = -60;
exce = 1.5; gabae = 27; glue = 27;

%%
I = @(v,m,n,h) gL.*(v-VL) + gK.*n.^4.*(v-VK) + gNa.*m.^3.*h.*(v-VNa);

mm = @(am,bm) am./(am+bm);
dmdt = @(m,a,b) 3.*(a.*(1-m) - b.*m);
am = @(v) 0.1.*(v+30)./(1-exp(-0.1*(v+30)));
bm = @(v) 4.*exp(-(v+55)/18);

dndt = @(n,a,b) 3.*(a.*(1-n) - b.*n);
an = @(v) 0.01.*(v+34)./(1-exp(-0.1*(v+34)));
bn = @(v) 0.125.*exp(-(v+44)/80);

dhdt = @(h,a,b) 3.*(a.*(1-h)-b.*h);
ah = @(v) 0.07.*exp(-(v+44)/20);
bh = @(v) 1./(1+exp(-0.01*(v+14)));

%% neuron
dvdt = @(iext,imem,isyn) iext-imem-isyn;
%% synaptic current
dsdt = @(v,s) (exp(-exp((-v+20)/4)/tau)-s)/tau;
kt = 0:dt:40;
kernal = exp(1)*kt/tau^2.*exp(-kt/tau)/5;
%% ajustable parameter

step = 30000;
tott = step*dt;
stimb = 100;stime = 250;


%% initiation inputs
for i = 1:1
    inputlayer(i).v = ones(neuron,step)*VL;
    inputlayer(i).n = ones(neuron,step)*0.18;
    inputlayer(i).m = ones(neuron,step)*0.05;
    inputlayer(i).h = ones(neuron,step)*0.2;
    inputlayer(i).se = zeros(neuron,step);
    inputlayer(1).Iext = zeros(neuron,step);
    inputlayer(i).Ibg = zeros(neuron,step);
    for n = 1:neuron
        inputlayer(i).Ibg(n,:) = mu_en+delta*wgn(1,step,2);
    end
    inputlayer(i).isspike = zeros(neuron,step);
    inputlayer(i).Isyn = zeros(neuron,step);
    inputlayer(i).output = zeros(1,step);
end
for n = 1:neuron
    inputlayer(1).Iext(n,ceil(stimb/dt):ceil(stime/dt)) = mu_ext + 0.05*mu_ext*wgn(1,ceil(stime/dt)-ceil(stimb/dt)+1,2);
end
for i = 2:step
    v = inputlayer(1).v(:,i-1);
    n = inputlayer(1).n(:,i-1);
    m = inputlayer(1).m(:,i-1);
    h = inputlayer(1).h(:,i-1);
    Iext = inputlayer(1).Iext(:,i-1);
    Ibg = inputlayer(1).Ibg(:,i-1);
    Isyn = 0;
    
    dn = dndt(n,an(v),bn(v))*dt;
    dh = dhdt(h,ah(v),bh(v))*dt;
    dv = dvdt(Iext,I(v,m,n,h),Isyn)*dt + Ibg*dt;
    
    inputlayer(1).m(:,i) = mm(am(v),bm(v));
    inputlayer(1).v(:,i) = inputlayer(1).v(:,i-1) + dv;
    inputlayer(1).n(:,i) = inputlayer(1).n(:,i-1) + dn;
    inputlayer(1).h(:,i) = inputlayer(1).h(:,i-1) + dh;
    
end
inputlayer(1).se(:,2:end) = diff(inputlayer(1).v>0,1,2)>0;
if neuron > 1
    inputlayer(1).output = conv(sum(inputlayer(1).se),kernal);
else
    inputlayer(1).output = conv(inputlayer(1).se,kernal);
end
inputlayer(1).output = inputlayer(1).output(1:step);

%% initiation network
% exclayer
for i = 1:totlayer
    exclayer(i).v = ones(neuron,step)*VL;
    exclayer(i).n = ones(neuron,step)*0.18;
    exclayer(i).m = ones(neuron,step)*0.05;
    exclayer(i).h = ones(neuron,step)*0.2;
    exclayer(i).se = zeros(neuron,step);
    exclayer(i).Iext = zeros(neuron,step);
    exclayer(i).Ibg = zeros(neuron,step);
    for n = 1:neuron
        exclayer(i).Ibg(n,:) = mu_en+delta*wgn(1,step,2);
    end
    exclayer(i).isspike = zeros(neuron,step);
    exclayer(i).Isyn = zeros(1,step);
    exclayer(i).output = zeros(1,step);
end
% inhlayer
for i = 1:totlayer
    inhlayer(i).v = ones(neuron,step)*VL;
    inhlayer(i).n = ones(neuron,step)*0.18;
    inhlayer(i).m = ones(neuron,step)*0.05;
    inhlayer(i).h = ones(neuron,step)*0.2;
    inhlayer(i).si = zeros(neuron,step);
    inhlayer(i).Iext = zeros(neuron,step);
    inhlayer(i).Ibg = zeros(neuron,step);
    inhlayer(i).isspike = zeros(neuron,step);
    for n = 1:neuron
        inhlayer(i).Ibg(n,:) = mu_inh+3*delta*wgn(1,step,2);
    end
    inhlayer(i).Isyn = zeros(step,neuron);
    inhlayer(i).output = zeros(step,1);
end
% motorlayer
for i = 1:totlayer
    motorlayer(i).v = ones(1,step)*VL;
    motorlayer(i).n = ones(1,step)*0.18;
    motorlayer(i).m = ones(1,step)*0.05;
    motorlayer(i).h = ones(1,step)*0.2;
    motorlayer(i).si = zeros(1,step);
    motorlayer(i).Iext = zeros(neuron,step);
    motorlayer(i).Ibg = zeros(1,step);
    motorlayer(i).isspike = zeros(1,step);
    for n = 1:1
        motorlayer(i).Ibg(n,:) = mu_mn+delta*wgn(1,step,2);
    end
    motorlayer(i).Isyn = zeros(step,1);
    motorlayer(i).output = zeros(step,1);
end

%% update
for l = 1:totlayer
    % update exclayer
    for i = 2:step
        v = exclayer(l).v(:,i-1);
        n = exclayer(l).n(:,i-1);
        m = exclayer(l).m(:,i-1);
        h = exclayer(l).h(:,i-1);
        Ibg = exclayer(l).Ibg(:,i-1);
        
        Isyn = exce*inputlayer(1).output(i-1).*(v-Ve)/neuron;
        
        if l>1
            Isyn = Isyn + gabae*inhlayer(l-1).output(i-1).*(v-Vi)/neuron;
        end
        dn = dndt(n,an(v),bn(v))*dt;
        dh = dhdt(h,ah(v),bh(v))*dt;
        dv = dvdt(0,I(v,m,n,h),Isyn)*dt + Ibg*dt;
        
        exclayer(l).m(:,i) = mm(am(v),bm(v));
        exclayer(l).v(:,i) = exclayer(l).v(:,i-1) + dv;
        exclayer(l).n(:,i) = exclayer(l).n(:,i-1) + dn;
        exclayer(l).h(:,i) = exclayer(l).h(:,i-1) + dh;
        
    end

    
    exclayer(l).se(:,2:end) = diff(exclayer(l).v>0,1,2)>0;
    if neuron > 1
        exclayer(l).output = conv(sum(exclayer(l).se),kernal);
    else
        exclayer(l).output = conv(exclayer(l).se,kernal);
    end
    exclayer(l).output = exclayer(l).output(1:step);
    
    % update motorlayer
    for i = 2:step
        v = motorlayer(l).v(i-1);
        n = motorlayer(l).n(i-1);
        m = motorlayer(l).m(i-1);
        h = motorlayer(l).h(i-1);
        Iext = motorlayer(l).Iext(i-1);
        Ibg = motorlayer(l).Ibg(i-1);
        
        Isyn = exce*exclayer(l).output(i-1).*(v-Ve)/neuron;
        
        dn = dndt(n,an(v),bn(v))*dt;
        dh = dhdt(h,ah(v),bh(v))*dt;
        dv = dvdt(Iext,I(v,m,n,h),Isyn)*dt + Ibg*dt;
        
        motorlayer(l).m(i) = mm(am(v),bm(v));
        motorlayer(l).v(i) = motorlayer(l).v(i-1) + dv;
        motorlayer(l).n(i) = motorlayer(l).n(i-1) + dn;
        motorlayer(l).h(i) = motorlayer(l).h(i-1) + dh;
    end
    motorlayer(l).si(:,2:end) = diff(motorlayer(l).v>0,1,2)>0;
    motorlayer(l).output = conv(motorlayer(l).si,kernal);
    motorlayer(l).output = motorlayer(l).output(1:step);
    motorlayer(l).output = circshift(motorlayer(l).output,2/dt);
    motorlayer(l).output(1,1:2/dt) = 0;
    
    % update inhlayer
    for i = 2:step
        v = inhlayer(l).v(:,i-1);
        n = inhlayer(l).n(:,i-1);
        m = inhlayer(l).m(:,i-1);
        h = inhlayer(l).h(:,i-1);
        Iext = inhlayer(l).Iext(:,i-1);
        Ibg = inhlayer(l).Ibg(:,i-1);
        
        Isyn =glue*motorlayer(l).output(i-1).*(v-Vi);
        
        dn = dndt(n,an(v),bn(v))*dt;
        dh = dhdt(h,ah(v),bh(v))*dt;
        dv = dvdt(Iext,I(v,m,n,h),Isyn)*dt + Ibg*dt;
        
        inhlayer(l).m(:,i) = mm(am(v),bm(v));
        inhlayer(l).v(:,i) = inhlayer(l).v(:,i-1) + dv;
        inhlayer(l).n(:,i) = inhlayer(l).n(:,i-1) + dn;
        inhlayer(l).h(:,i) = inhlayer(l).h(:,i-1) + dh;
    end
    
    inhlayer(l).si(:,2:end) = diff(inhlayer(l).v>0,1,2)>0;
    if neuron > 1
        inhlayer(l).output = conv(sum(inhlayer(l).si),kernal);
    else
        inhlayer(l).output = conv(inhlayer(l).si,kernal);
    end
    inhlayer(l).output = inhlayer(l).output(1:step);
end
readout = zeros(l,step);
for l = 1:totlayer
    v = motorlayer(l).v;
    idx = find(diff(v>0)>0);
    readout(l,idx) = 1;
end


end