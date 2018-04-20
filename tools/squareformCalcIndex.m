function ixMAT = squareformCalcIndex(ixTRI,N)
% Calculates square matrix index from vectorized triangular half of distance matrix,
% such as output of pdist
% 
% Ex. Find indices of furthest 2 points in X:
% dists = pdist(X,'Euclidean');
% n = numel(dists);
% maxdistIndex = squareformCalcIndex(maxix(dists),n);
% [maxdistPoint1,maxdistPoint2] = ind2sub([n n], maxdistIndex);
% X([maxdistPoint1;maxdistPoint2],:) % the points
% 
% 2018-04-17 AZ Created

m  = ceil(sqrt(N*2));             % Dimensionality of square matrix
IX = tril(reshape(1:m^2,m,m),-1); % Generate indices of square matrix
IX = IX(IX>0);                    % Keep indices only for upper triangular minus diagonal
ixMAT = IX(ixTRI);                % Convert input index to output

return

