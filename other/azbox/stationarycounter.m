function stationarycounter(n,N)
% STATIONARYCOUNTER : replaces text in screen with number n (string)
% 
% ex.
% stationarycounter(2)            --> '2'
% stationarycounter(2,5)          --> '2/5'
% stationarycounter([2 3],[4 5])  --> '8/20'
% 
% 2012-03-12 AZ Created
% 2012-04-09 AZ Added argument for total (N)
% 2014-10-15 AZ Automatically calculates index for 2-level nested loop


if nargin == 2
   if numel(n) == 2
      n = (n(1)-1)*N(2) + n(2);
   end
   
   N = prod(N);
   total = sprintf('/%i',N);
   ndel = ceil(log10(n+eps)) + ceil(log10(N+1+eps)) + 1;
else
   ndel = ceil(log10(n+eps));
   total = [];
end

if n < 2
	fprintf(repmat(' ',[1 ndel]));
end

fprintf([repmat('\b',1,ndel) '%i%s'],n,total);