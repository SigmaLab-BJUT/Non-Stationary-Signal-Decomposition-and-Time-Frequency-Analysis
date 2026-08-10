function IA = estimateIA(IB,s0,wd,lambda)
if nargin<4;lambda = 8*1e-2;end
x = exp(cumsum(IB));
k = max(s0./(x+8*1e-2))*ones(length(s0),1);
%IA = k.*x;
for i = wd+1:length(s0)-wd
    k(i) = max(s0(i-wd:i+wd)./(x(i-wd:i+wd)+lambda));
end
IA = mean(k)*x;