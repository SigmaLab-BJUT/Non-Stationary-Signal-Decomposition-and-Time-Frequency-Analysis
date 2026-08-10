function [Vs0 Vs Ts Rs] = sst(s,gamma,sigma,ft,bt)
% computes the STFT Vs, the SST Ts and the reassigned STFT Rs of the signal s.
% Uses a Gaussian window with SD sigma. Computes the TF transforms at
% frequencies ft and discrete times bt.
%
% Input:
%   s: signal
%   gamma: threshold
%   sigma: window parameter
%   ft: index of frequencies where to compute the transforms.
%       default: 1:N/2 (N=length(s))
%   bt: time index where to compute the transforms
%       default: 1:N

% checking length of signal
n = length(s);
% n = 512;
nv = log2(n);
if mod(nv,1)~=0
    warning('The signal is not a power of two, truncation to the next power');
    s = s(1:2^floor(nv));
end
n = length(s);
s = s(:);

% Optional parameters
if nargin<5
   ft = 1:n/2;
   bt = 1:n;
end

% Padding
sleft = flipud(s(2:n/2+1));
sright = flipud(s(end-n/2:end-1));
x = [sleft; s; sright];
clear sleft sright;
    
%% STFT and operators tau and omega
nb = length(bt);
neta = length(ft);

t = linspace(-0.5,0.5,n);t=t';

g = exp(-pi/sigma^2*t.^2);
gp = -2*pi/sigma^2*t .* g; % g'

Vs0 = zeros(n,nb);
Vs = zeros(neta,nb);
Ts = zeros(neta,nb);
Rs = zeros(neta,nb);

%% Direct reassignment
df = ft(2)-ft(1);
db = bt(2)-bt(1);
for b=1:nb
    % Computing Vs, omega and tau in t=b
    tmp = fft(x(bt(b):bt(b)+n-1).*g)/n/sigma;
    vtmp = tmp(ft);
    Vs0(:,b) = tmp*n*sigma;
%     if b == floor(nb/2)
%         1;
%     end
    tmp = fft(x(bt(b):bt(b)+n-1).*gp)/n/sigma;
    omegtmp = (ft-1)'-real(tmp(ft)/2/1i/pi./vtmp);
    
    
    
    
    tmp = fft(x(bt(b):bt(b)+n-1).*g)/1/sigma;
    tmp = [tmp(2)-tmp(1); (tmp(3:end)-tmp(1:end-2))/2; tmp(end)-tmp(end-1)];
    tautmp = bt(b)-1+real(tmp(ft)/2/pi/1i ./vtmp);
    
    
    % Storing values
    
    Vs(:,b) = vtmp;
    
    % Reassignment step
	for eta=1:neta
        if abs(vtmp(eta))>gamma
            k = 1+round((omegtmp(eta)-ft(1)+1)/df);
            if k>=1 && k<=neta
                Ts(k,b) = Ts(k,b) + vtmp(eta)*exp(1i*pi*(ft(eta)-1));
                l = 1+round((tautmp(eta)-bt(1)+1)/db);
                if l>=1 && l<=nb
                    Rs(k,l) = Rs(k,l) + abs(vtmp(eta));% 这里对Rs的计算应该有错误，
                                                       % 因为这里只对时间进行压缩了
                end
            end
        end
    end
end
Ts = df*sigma*Ts;
Rs = Rs * df * db;

end