function [u,v,p,q] = extract_one_component( s, wd, n1, n2, beta, eps )

s = s(:);
k0 = 6;M0 = 33;
k1 = 6;M1 = 9;
k2 = 6;M2 = 13;
maxIter = 100;
v = s;
vv = zeros(length(s),maxIter);
for ii = 1:maxIter
    [p,q] = updatePQ(v, zeros(size(v)), zeros(size(v)), k0,M0,k1,M1,k2,M2, wd,n1,n2,20);
    t0 = 100:(length(s)-100);
    t = 1:length(p);
    [p] = polyfitting(t0',p(t0),2,t');
    [V,~,~] = nspWithTrueT(s,p,q,s,5,beta);
    v = V(:,end);
    if mod(ii+1,3)==0
        if M1>k1+2;M1 = M1-2;end
        if M2>k2+4;M2 = M2-2;end
        if n1<300;n1 = n1+10;end
    end
    
    if ii>1 && (vv(:,ii-1)-vv(:,ii))'*(vv(:,ii-1)-vv(:,ii))<eps*(s'*s)
        break;
    end    
    vv(:,ii+1) = v;
%     STFT(v);title(num2str(ii));drawnow;
end
u = s-v;
