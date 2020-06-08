function stationarycounter(n,N)
% STATIONARYCOUNTER : replaces text in screen with number n (string).
% 
% Inputs can be vectors of indices/limits of nested loop, where from left to
%    right, they go from outer to inner, e.g.:
% stationarycounter([outern middlen innern],[outerN middleN innerN])
% 
% ex.
% stationarycounter(2)                --> '2'
% stationarycounter(2,5)              --> '2/5'
% stationarycounter([2 3],[4 5])      --> '8/20'
% stationarycounter([2 3 4],[4 5 6])  --> '46/120'
% 
% 2012-03-12 AZ Created
% 2012-04-09 AZ Added argument for total (N)
% 2014-10-15 AZ Automatically calculates index for 2-level nested loop
% 2020-06-08 AZ Generalized for aribtrary-level nested loop

if nargin == 2
   if numel(n) > 1
   	[n,N] = calcIter([0 n],N);
   end
   
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

if n == N, fprintf('\n'); end

return
end


%% helper function
function [n,N] = calcIter(n,N)

	if numel(n) > 1
		[n,N] = calcIter([n(1)+(n(2)-1)*prod(N(end-numel(n)+3:end)) n(3:end)],N);
	else
		n = n + 1;
	end

	N = prod(N);
	return
end
