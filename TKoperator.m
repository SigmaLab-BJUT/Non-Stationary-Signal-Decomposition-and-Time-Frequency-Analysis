function [IF,IA] = TKoperator(s,wd)
if nargin<2;wd = 0;end
s = [s(wd+1:-1:2);s;s(end-1:-1:end-wd)];
y = s(2:end-1).^2-s(3:end).*s(1:end-2);
y = [y(1);y;y(end)];

if 0
x = diff(s);
x = [x(1);x];
z = x(2:end-1).^2-x(3:end).*x(1:end-2);
z1 = [z(1);z;z(end)];
x = diff(s);
x = [x;x(end)];
z = x(2:end-1).^2-x(3:end).*x(1:end-2);
z2 = [z(1);z;z(end)];
z = z1/2+z2/2;
IF = real(acos(1-z./y/2));
IA = real(sqrt(y./(1-(1-z./y/2).^2)));
else
x = s(3:end)-s(1:end-2);
x = [x(1);x;x(end)];
z = x(2:end-1).^2-x(3:end).*x(1:end-2);
z = [z(1);z;z(end)];

IF = real(asin(sqrt(z./y/4)));
IA = real(2*y./sqrt(z));
end


IF = IF(wd+1:end-wd);
IA = IA(wd+1:end-wd);
% IF = asin(real(sqrt(y)));
% IA = sqrt(b);
