function insR = instantaneousRate(T,spikeLoc,sig,varargin)
% this funciton using Gaussian smoothing method to find the instantaneous
% firing rate
% T          a vector where the rate shell be calculated
% spikeLoc   a vector store the time where a spike happen
% sig        the smooth time window used, in unit of s

gaussFlt = @(t,sigma) 1/(sigma*sqrt(2*pi))*exp(-t.^2/(2*sigma^2));
FRfunc = @(t,T,sigma) sum(gaussFlt(t-T,sigma));

insR = zeros(length(T),1);
for i0 = 1:length(T)
     insR(i0) = sum(FRfunc(T(i0),spikeLoc,sig));
end

end