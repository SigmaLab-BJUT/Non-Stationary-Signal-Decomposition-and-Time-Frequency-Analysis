function [p,q, snr] = updatePQ(s0,p0,q0,k0,M0,k1,M1,k2,M2,wd,n1,n2,Iter,ptrue,qtrue)
%if nargin<13;beta = 0.1;end
if ~exist('ptrue','var');ptrue = zeros(length(p0),1);end
if ~exist('qtrue','var');qtrue = zeros(length(q0),1);end
s0 = [s0(wd+1:-1:2);s0;s0(end-1:-1:end-wd)];
q0 = [q0(wd+1:-1:2);q0;q0(end-1:-1:end-wd)];
p0 = [p0(wd+1:-1:2);p0;p0(end-1:-1:end-wd)];
[D0,~,~] = sgolayfiltDesign(k0,M0,length(s0));%smoothing
[~,D1,~] = sgolayfiltDesign(k1,M1,length(s0));
[~,~,D2] = sgolayfiltDesign(k2,M2,length(s0));
D1 = sparse(D1);
D2 = sparse(D2);
D0 = sparse(D0);



s1 = D1*s0;
s2 = D2*s0;


[~,D1,~] = sgolayfiltDesign(8,9,length(s0));
D1(1:3,:) = 0;
D1(end-2:end,:) = 0;

%s0 = D0*s0;
q = q0;
p = p0;
i = 0;
eps = 1e-12;
while i<Iter
    [p,~,snr1] = updateP(s0,s1,s2,p,q,D1,n1,ptrue,wd);
    [q,~,snr2] = updateQ(s0,s1,s2,p,q,D1,n2,qtrue,wd);
    if (p0-p)'*(p0-p)<eps*(s0'*s0)
        break;
    end
    i = i+1;
    p0 = p;
    snr(i).p = snr1;
    snr(i).q = snr2;
end
p = p(wd+1:end-wd);
q = q(wd+1:end-wd);

end

