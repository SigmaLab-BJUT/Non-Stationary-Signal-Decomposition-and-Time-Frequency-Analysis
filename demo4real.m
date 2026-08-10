% loading real life signal
load('zai.mat')
s00 = s(1:2:end);
s = s00(100:800);
n = length(s);
t = (1:n)';
dt = t(2)-t(1);

% parameters setting
maxiter = 100;
wd = 0;
gamma = 0;
lambda = 1;
n1 = 5;
n2 = 7;

rv = zeros( n, 5 );
rq = zeros( n, 5 );
rp = zeros( n, 5 );
eps = [1e-6,1e-6,1e-6,1e-6,1e-8];
beta = [0.1,0.2,1,1,1];

v = s;
for i = 1:5
    [u,v,rp(:,i),rq(:,i)] = extract_one_component( v, wd, n1, n2, beta(i), eps(i) );
    rv(:,i) = v;
    v = u;
%     STFT(v);title(num2str(ii));drawnow;
end

figure; 
for i = 1:5
    subplot(5,1,i); plot( rv(:,i) );
end
