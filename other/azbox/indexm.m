function M = indexm(M,varargin)
% index.m: Returns matrix indexed by those given
%
% e.g., given a 5-D matrix M,
% M = index(M,a1,[],a3,a4);
% is equivalent to
% M = M(a1,:,a3,a4,:);
%
% 2015-04-01 AZ Created

for v = (nargin-1):-1:1
   if ~isempty(varargin{v})
      sz = size(M);
      n  = ndims(M);
      
      perm1 = [v setdiff(1:n,v)];
      perm2 = [perm1(perm1<v)+1 1 perm1(perm1>v)];

      % make v'th dimension first in 2D array
      M = permute(M,perm1);
      M = reshape(M,sz(v),prod(sz(perm1(2:end))));
      
      if     ~any(isnan(varargin{v}(:)));
         % apply index to v'th dimension
         M = M(varargin{v},:);
         % put back v'th dimension in original array
         M = reshape(M,[size(M,1) sz(perm1(2:end))]);
         M = permute(M,perm2);
      elseif ~all(isnan(varargin{v}(:)))
         % if forcing conflicting rebinning, do element-wise
         newM = NaN(size(varargin{v}));
         for i = 1:numel(varargin{v})
            if ~isnan(varargin{v}(i));   newM(i) = M(varargin{v}(i));   end
         end
         M = newM;
      end
   end
end
