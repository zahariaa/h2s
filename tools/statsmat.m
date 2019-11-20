function mat = statsmat(uppert,lowert,diagonal)
% mat = statsmat(uppert,lowert,<diagonal=zeros>)
% Create matrix of stats of interest that do pairwise comparisons (upper and
% lower triangles), and of the values themselves (diagonal)
%
% 2018-09-25 AZ Created

if ~exist('lowert'  ,'var'), lowert = zeros(size(uppert)); end
if ~exist('uppert'  ,'var'), uppert = zeros(size(lowert)); end
if ~exist('diagonal','var') || isempty(diagonal)
   n = ceil(sqrt(2*numel(uppert)));
   diagonal = ones(n,1);
else
   n = numel(diagonal);
end

% Fill in matrix
mat = diag(diagonal);
mat(triu(true(n), 1)) = uppert;
mat(tril(true(n),-1)) = lowert;

