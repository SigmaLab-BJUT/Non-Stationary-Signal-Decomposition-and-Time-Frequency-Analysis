function [x1,y1,z1] = filter_mean(x,y,z,x0,y0,z0,r)
R = sqrt((x-x0).^2+(y-y0).^2+(z-z0).^2);
X = [];
for i = 1:size(x,1)
    for j = 1:size(x,2)
        if R(i,j)<r
            X = [X;[x(i,j),y(i,j),z(i,j)]];
        end
    end
end
X = mean(X);
x1 = X(1);
y1 = X(2);
z1 = X(3);
