%% npr demo
T1 = [505, 495; 495, 505];
T2 = 0.05*ones(2);
T = [T1, zeros(2); zeros(2), T2];
noi = 0.004*(rand(4)-0.5);
Th = T + noi;
% Th = T;
b = [sqrt(2), 0, sqrt(2), 0]';
epsi = 1e-3;

r0 = zeros(4,1);
iters = 300;
lambda = zeros(iters,1);
err = zeros(iters,1);
rerr = zeros(iters,1);
lambda0 = norm(Th);

for i = 1:iters
    r1 = (Th'*Th + lambda0*eye(size(Th))) \ (Th'*b + lambda0*r0);
    lambda(i) = norm( Th*( r1-r0 ) ) / norm( r1-r0 );
    lambda0 = lambda(i);

    rerr(i) = norm(r1-r0);
    err(i) = norm( Th*r1 - b ); 
    rprev = r0;
    if i > 4
        if abs(rerr(i) - mean(rerr(i-3:i-1))) < epsi && rerr(i) > rerr(i-1)
            break;
        end
    end
    r0 = r1;
end
        
figure;
subplot(2,1,1); plot( log10(lambda(1:i)) );
subplot(2,1,2); plot( rerr(1:i) );
   
%% separation demo
n = 1000;
t = linspace(0,1,n)';
dt = t(2)-t(1);
f1 = 80;
f2 = 65; % 50-65 where delta = 0.6; 60 when delta = 0.01-0.1
delta = 0.6;
A1 = 1.0+0.2*cos(2*pi*t);
A2 = 1.0+0.4*cos(3*pi*t);
MaxIter = 10;

A2 = delta*A2;

phi1 = 2*pi*(f1*t+1*cos(4*pi*t));
phi2 = 2*pi*(f2*t+1*cos(4*pi*t));
s1 = A1.*cos(phi1);
s2 = A2.*cos(phi2);

p01 = 1./((dt*2*pi*(f1-1*sin(4*pi*t)*4*pi)).^2);
p02 = 1./((dt*2*pi*(f2-1*sin(4*pi*t)*4*pi)).^2);
q01 = -0.2*sin(2*pi*t)*2*pi*dt./A1;
q02 = -0.2*sin(3*pi*t)*3*pi*dt./A2;

s0 = s1+s2;

gamma = 0;
lambda = 1;
k0 = 6; M0 = 33;
k1 = 6; M1 = 9;
k2 = 6; M2 = 13;
n1 = 20;
n2 = 6;
wd = 0;
p = zeros(size(s0));
q = p;

v0 = s0;
i = 1;
verr = zeros(MaxIter,1);
while i < MaxIter
    [p,q0,~] = updatePQ(v0,zeros(size(p)),zeros(size(q)),k0,M0,k1,M1,k2,M2,wd,n1,n2,2);
    [~,q,~] = updatePQ(v0,zeros(size(p)),zeros(size(q)),k0,M0,k1,M1,k2,M2,160,n1,n2,2);
    [v,~,~] = nspWithTrueT(s0,p,q0);

    if i == 1
        verr(i) = 1;
    else
        verr(i) = norm( v-v0 ) / norm( s0 );
    end
    
    fprintf( 'SNRs:%.4f\t SNRi:%.4f\t SNRa:%.4f\t vErr:%.5f\n', ...
            snr_compute(s1,v(wd+1:end-wd)), ...
            snr_compute(p01,p(wd+1:end-wd)), ...
            snr_compute(q01,q(wd+1:end-wd)), ...
            verr(i) );
    
    if i > 4 && abs( verr(i) - mean(verr(i-3:i-1)) ) < 1e-2
        break;
    else
        i = i + 1;
        v0 = v;
    end    
end
IA = estimateIA(q,v,50);

figure;
subplot(3,1,1); plot( t, s0 ); title('input signal');
subplot(3,1,2); plot( t, v ); hold on; plot( t, IA, 'r' );
title('1st component and its IA');
subplot(3,1,3); plot( t, 1./sqrt(p01), 'r' ); hold on; plot( t, 1./sqrt(p), 'b--' ); 
axis([0 1 0 1]); title('estimated IF');




