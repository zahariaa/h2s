function [points,categories] = randnsimplex_of_nballs(n,N,r)
% points = randnsimplex_of_nballs(numDimsOrCenterCoords,Num_points)
%
% 2018-06-11 AZ Created

%% Process inputs
if numel(n)==1, centers = nsimplex(n)';
else            centers = n;
end
[nballs,n] = size(centers);
if ~exist('r','var') || isempty(r)
   r = ones(nballs,1);
end

if numel(N)==1, N = repmat(N,[nballs 1]); end
N = [0 N(:)'];

%% Generate random points
points  = randnball(sum(N),n);

for iball = 1:nballs
   ix = sum(N(1:iball))+1:sum(N(1:iball+1));
   points(ix,:) = r(iball)*points(ix,:) + repmat(centers(iball,:),[N(iball+1) 1]);
end

categories = Categories(N(2:end));

return

