function demo_hyperdist2dist(lodim)
% demo_hyperdist2dist(lodim=2)
% Demo visualizing high dimensional vector magnitudes in lodim (e.g., 2D),
% for Gaussian and Uniform n-ball distributions of varying dimensionality
% 
% 2018-04-19 AZ Created

%% Preliminaries
ds = 2.^(1:9); n = 100;
types = {'Uniform','Gaussian'};
if ~exist('lodim','var') || isempty(lodim),   lodim = 2;
elseif lodim==3
   ds(1) = 3;
end

%% Draw loop
for t = 1:2
   figure;set(gcf,'Name',sprintf('%s%usamples',types{t},n));
   i = 0;
   for d = ds
      i = i+1;
      switch t
         case 1,   X = randnball(n,d);
         case 2,   X = randn(    n,d);
      end
      H = hyperdist2dist(X,lodim);
      subplot(3,3,i);
      switch lodim
         case 2,   plot( H(:,1),H(:,2)       ,'ko');drawCircle;
         case 3,   plot3(H(:,1),H(:,2),H(:,3),'ko');draw3dEllipsoid; lighting gouraud;camlight;
      end
      r = sqrt(2)*exp(gammaln((d+1)/2)-gammaln(d/2));
      if     t==2 && lodim==2,   drawCircle(     [],r); 
      elseif t==2 && lodim==3,   draw3dEllipsoid([],r*eye(3));
      end
      title(sprintf('%u dimensions',d))
   end
end


