function [D0,D1,D2] = sgolayfiltDesign(k,M,N)
D0 = zeros(N,N);
D1 = D0;
D2 = D0;
for n = (M-1)/2+1:N-(M-1)/2
    [aa,g] = sgolay(k,M);
    D0(n,n-(M-1)/2:n+(M-1)/2) = g(:,1)';
    D1(n,n-(M-1)/2:n+(M-1)/2) = g(:,2)';
    D2(n,n-(M-1)/2:n+(M-1)/2) = 2*g(:,3)';
end
if 1
for n = 2:(M-1)/2
    W = 2*(n-1)+1;
    j = min([k,W-1]);
    [aa,g] = sgolay(j,W);
    D0(n,n-(W-1)/2:n+(W-1)/2) = g(:,1)';
    D1(n,n-(W-1)/2:n+(W-1)/2) = g(:,2)';
    D2(n,n-(W-1)/2:n+(W-1)/2) = 2*g(:,3)';
    
    D0(N-n+1,N+1-n-(W-1)/2:N+1-n+(W-1)/2) = g(:,1)';
    D1(N-n+1,N+1-n-(W-1)/2:N+1-n+(W-1)/2) = g(:,2)';
    D2(N-n+1,N+1-n-(W-1)/2:N+1-n+(W-1)/2) = 2*g(:,3)';
end
D0(1,1) = 1;
D1(1,1:2) = [-1,1];
D2(1,1:3) = [1 -2 1];

D0(N,N) = 1;
D1(N,N-1:N) = [-1,1];
D2(N,N-2:N) = [1 -2 1];
else
    m = (M-1)/2;
    X = zeros(M,k+1);
    Y = zeros(M,1);
    for r = 1:M
        for j = 0:k
            X(r,j+1) = sum([-m:m].^(r+j));
        end
    end
    
    alpha = (X'*X)\(X'*Y);
end