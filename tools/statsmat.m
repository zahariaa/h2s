function mat = statsmat(uppert,lowert,diagonal)
% mat = statsmat(uppert,lowert,diagonal)
% Create matrix of stats of interest that do pairwise comparisons (upper and
% lower triangles), and of the values themselves (diagonal)
%
% 2018-09-25 AZ Created

if ~exist('lowert'  ,'var'), lowert = zeros(size(uppert)); end
if ~exist('uppert'  ,'var'), uppert = zeros(size(lowert)); end
if ~exist('diagonal','var')
   n = ceil(sqrt(2*numel(uppert)));
   diagonal = zeros(n,1);
else
   n = numel(diagonal);
end

ix = nchoosek_ix(n);

% Fill in matrix
mat = diag(diagonal);

for i = 1:numel(uppert)
   mat(ix(1,i),ix(2,i)) = lowert(i);
   mat(ix(2,i),ix(1,i)) = uppert(i);
end

