function [v,snr,lambda] = nspWithTrueT(s0,p,q,s1,iter,beta)
if nargin<4;s1 = zeros(size(s0));end;
if ~exist('beta','var');beta = 1;end
if ~exist('iter','var');iter = 100;end
I = eye(length(p));
[~,D1,D2] = sgolayfiltDesign(8,9,length(p));
P = diag(p);
Q = diag(q);
T = P*D2+(1/2*diag(D1*p)-2*P*Q)*D1+(I-P*diag(D1*q)-1/2*diag(D1*p)*Q+P*Q*Q);
T(1:4,:) = 0;
T(end-3:end,:)=0;

D = T'*T;

vv0 = s0;
lambda = 1*ones(10,1);
i = 1;
eps = 1e-2;
verr = zeros(iter,1);
V = zeros(length(s0),iter);
while i < iter
    W = chol(D+lambda(i)*I);
    vv = W\(W'\(lambda(i)*vv0));
    lambda(i+1) = beta*((T*(vv-vv0))'*(T*(vv-vv0))/((vv-vv0)'*(vv-vv0)));
    V(:,i) = vv;

    if i == 1
        verr(i) = norm( vv );
    else
        verr(i) = norm( vv-vv0 );
    end

    snr(i) = snr_compute(s1,[vv0]);
%     fprintf( '%.4f\n', ((vv-vv0)'*(vv-vv0)) / (s0'*s0) );
    if i > 2 && verr(i) < eps*(s0'*s0) && verr(i) > verr(i-1)
        v = V(:,i-2);
        break;
    else
        i = i+1;
        vv0 = vv;
    end
end
