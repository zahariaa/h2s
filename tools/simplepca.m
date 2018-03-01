function pointsLow = simplepca(points,dimLow)
% pointsLow = simplepca(points,<dimLow=2>)
% 
% points should be [n X d]
%
% AZ Created 2018-03-01

[n,d] = size(points);

if ~exist('dimLow','var') || isempty(dimLow),   dimLow = 2;   end

if d <= dimLow
   pointsLow = [points zeros(n,dimLow-d)];
   warning('Input dimensionality (%u) <= requested (%u). Zero padded.',d,dimLow);
   return
end

[V,D] = eig(points'*points);

pointsLow = points*V(:,(d-dimLow+1):d);

return

