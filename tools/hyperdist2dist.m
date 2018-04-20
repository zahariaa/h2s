function lo = hyperdist2dist(hi,dim)
% lo = hyperdist2dist(hi,dim=2)
% Visualizes high dimensional vector magnitudes with
% magnitudes with 2 (default) or 3 dimensional magnitudes
% 
% 2018-04-19 AZ Created

if ~exist('dim','var') || isempty(dim),   dim = 2;   end

radii.hi = sqrt(sum(hi         .^2,2));
radii.lo = sqrt(sum(hi(:,1:dim).^2,2));

lo = repmat(radii.hi,[1 dim]).*hi(:,1:dim)./repmat(radii.lo,[1 dim]);

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

