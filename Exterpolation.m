function [y, beta] = Exterpolation(y,wd,k)
X = ones(length(y)-2*wd,k+1);
T = [wd+1:length(y)-wd]';
for i = 2:k+1
    X(:,i) = T.^(i-1);
end
beta = (X'*X)\(X'*y(wd+1:end-wd));

X = ones(length(y),k+1);
T = [1:length(y)]';
for i = 2:k+1
    X(:,i) = T.^(i-1);
end
y = X*beta;