function x = logify(x)
% logify: convert to log10 axis if data is on it already
% useful for avoiding openGL errors from set(gca,'XScale','log')
% 
% 2016-03-30 AZ Created

rsz = 1:size(x,2);

for i = rsz(:)'
   ix = ~isnan(x(:,i));
   if numel(unique(diff(log10(x(ix,i)),1)))==1 || var(diff(log10(x(ix,i)),1))<0.01
      x(:,i) = log10(x(:,i));
   end
end

