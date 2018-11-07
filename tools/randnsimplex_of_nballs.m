function [points,categories] = randnsimplex_of_nballs(n,N)
% points = randnsimplex_of_nballs(num_of_dimensions,Num_points)
%
% 2018-06-11 AZ Created

if numel(N)==1, N = repmat(N,[n+1 1]); end
N = [0 N(:)'];

points  = randnball(sum(N),n);
centers = nsimplex(n)';

for iball = 1:n+1
   ix = sum(N(1:iball))+1:sum(N(1:iball+1));
   points(ix,:) = points(ix,:) + repmat(centers(iball,:),[N(iball+1) 1]);
end

categories = Categories(N(2:end));

return

