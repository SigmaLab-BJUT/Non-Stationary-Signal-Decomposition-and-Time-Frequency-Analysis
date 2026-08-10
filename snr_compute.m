function y = snr_compute(x0,x)
if size(x0,1)<size(x0,2)
    x0 = conj(x0');
end
if size(x,1)<size(x,2)
    x = conj(x');
end
[r,c] = size(x);
for i = 1:size(x0,2)
x00=x0(:,i);
p_x0 = x00'*x00;
x00 = x00*ones(1,c);
p_noise = diag((x-x00)'*(x-x00))';
y(i,:) = 10*log10(p_x0./p_noise);
end
% if nargout<1
%     disp(num2str(y));
% end