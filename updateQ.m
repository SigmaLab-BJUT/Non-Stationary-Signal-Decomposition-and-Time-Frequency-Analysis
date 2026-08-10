function [q,gamma,snr] = updateQ(s0,s1,s2,p,q,D1,n1,qtrue,wd)
if ~exist('qtrue','var');qtrue = zeros(length(q0),1);end
n = length(s0);
snr = 0;
x1 = 8;
T = [1,floor(linspace(x1,length(s0)-x1+1,n1)),length(s0)];
for i = 1:length(T)
    y = zeros(length(T),1);
    y(i) = 1;
    PHI(:,i)=interp1(T,y,1:length(s0),'spline');
end
tau = 1;
%gamma = 1e-2;
gamma = 0;
q0 = q;%q0 上一步值
%q110 = zeros(length(T),1);
q10 = (PHI'*PHI)\(PHI'*q);% 上一步q的插值系数
%q10 = zeros(length(T),1);
C1 = diag(-s0.*p);
E = -((s2).*p+1/2*(s1).*(D1*p)+s0)*tau;
ii = 0;
q0temp = q0;
while ii<10
    C2 = diag(-2*(s1).*p+1/tau*s0.*p.*q0temp-1/2*s0.*(D1*p));
    W = C1*D1+C2;
    W(1:4,:) = 0;
    W(end-3:end,:)=0;
    E(1:4) = 0;
    E(end-3:end)=0;
    W = W*PHI;
    %W = sparse(W);
    q1 = (W'*W+gamma*(PHI'*PHI))\(W'*E+gamma*(PHI'*PHI)*q10);
    q0temp = PHI*q1;
    ii = ii+1;
end
q = PHI*q1;
q = q/tau;

