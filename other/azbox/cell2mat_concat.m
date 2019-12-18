function M = cell2mat_concat(C,d)
% Minifunction that concatenates cell elements into a vector/matrix
% d = dimension number to concatenate [DEFAULT = 1]
% 
% Code:
% M = cat(d,C{:});
% 
% 2012-11-13 AZ Created
% 2013-02-12 AZ Generalized concatenation dimension

if ~exist('d','var') || isempty(d)
   d = 1;
end

try
   M = cat(d,C{:});
catch
   sz = num2cell(size(C));
   for i = setdiff(1:numel(sz),d)
      sz{i} = ones(sz{i},1);
   end
   M = cellfun(@(x) cat(1,x{:}),mat2cell(C,sz{:}),'UniformOutput',false);
end