function [F,ix] = getFurthestPoints(X,p)
% [F,ix] = getFurthestPoints(X,p=1)
% Finds the p furthest points in X and returns them in F, where F=X(ix,:)
% 
% 2018-04-18 AZ Created

%% Preliminaries
if ~exist('p','var') || isempty(p),   p = 1;   end
[n,d] = size(X);

%% Find the largest distances
dists = pdist(X,'Euclidean');
maxDistsIX = sortix(dists,'descend');

%% Find the points with the largest distances
pts1=[];pts2=[];i=0;
while numel(unique([pts1 pts2])) < p
   i = i+1;
   % Convert from vectorized triangular matrix to square matrix subscript
   [pts1(i),pts2(i)] = ind2sub([n n], squareformCalcIndex(maxDistsIX(i),numel(dists)));
end

%% Prune points if more are found
ix = unique([pts1(:);pts2(:)]);
ix = ix(1:p);

%% Return the points too
F = X(ix,:);

return

end

