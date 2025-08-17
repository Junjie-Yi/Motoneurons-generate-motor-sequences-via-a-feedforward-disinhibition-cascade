% random delay 8-12v
data = readtable("fig4e_data.xlsx");
data = data{:,1};
dt = 0.01;
t = 0.01:0.01:1000;
cpg = sin(2*pi/150*(t));
thre = sin(180/360*2*pi);
cpg_on = find(cpg>thre);

t0 = rand(1000000,1)*300;
synapticDelay = gamrnd(8,3,1000000,1);

t_8 = t0;
t_12v = t0+synapticDelay;
t_12v_on = zeros(size(t_12v));
for ii = 1:length(t_12v)
    t_12v_on(ii) = t(cpg_on(find(cpg_on>t_12v(ii)/dt,1,'first')));
end

binwidth = 10;
delay_8_12v = t_12v_on-t_8;
xx = 0:1:200;
for ii = 1:length(xx)-1
  yy(ii) = sum(delay_8_12v>=xx(ii)&delay_8_12v<xx(ii+1))/length(delay_8_12v)*binwidth;
end
figure;hold on

histogram(data,0:binwidth:150,'Normalization','probability')
plot(xx(1:150),yy(1:150),'k','linewidth',1.5)