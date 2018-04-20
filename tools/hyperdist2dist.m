function lo = hyperdist2dist(hi,dim)
% lo = hyperdist2dist(hi,dim=2)
% Visualizes high dimensional vector magnitudes with
% magnitudes with 2 (default) or 3 dimensional magnitudes
% 
% Assumes input hi is a [samples x dimensions] matrix, or
% [samples x dimensions x ... more samples ... ] tensor,
% like sampleSpheres output
% 
% 2018-04-19 AZ Created

% If hi is a tensor, reshape it to a matrix
sz = size(hi);
if numel(sz)>2
   hi = reshape(permute(hi,[1 3:numel(sz) 2]),[],sz(2));
end

if ~exist('dim','var') || isempty(dim),   dim = 2;   end

radii.hi = sqrt(sum(hi         .^2,2));
radii.lo = sqrt(sum(hi(:,1:dim).^2,2));

lo = repmat(radii.hi,[1 dim]).*hi(:,1:dim)./repmat(radii.lo,[1 dim]);

% If hi was a tensor, reshape lo back to lower dimensional tensor
if numel(sz)>2
   lo = permute(reshape(lo,[sz(1) sz(3:end) dim]),[1 numel(sz) 2:(numel(sz)-1)]);
end

return


%% DEBUG/DEMO
d = 2.^(1:9); n = 1000;
figure;
for i = 1:9
   %X = randnball(n,d(i));
   X = randn(    n,d(i));
   H = hyperdist2dist(X);
   subplot(3,3,i);plot(H(:,1),H(:,2),'ko');
   drawCircle;
end

