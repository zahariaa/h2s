function ix = nchoosek_ix(n,k)
% ix = nchoosek_ix(n,k)
% creates indices for nchoosek combinations
% NOTE: only works for nchoose2 right now
% 
% 2018-09-20 AZ Created

if ~exist('k','var') || isempty(k), k = 2; end

N = nchoosek(n,k);

ix = zeros(k,N);

currix = 0;
for i = 2:n
   ix(1,currix+1:currix+n-i+1) = i:n;
   ix(2,currix+1:currix+n-i+1) = i-1;
   currix = currix + n-i+1;
end


