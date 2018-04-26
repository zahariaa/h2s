function demo_hyperdist2dist(lodim,JITTER,FAST)
% demo_hyperdist2dist(lodim=2,JITTER=false,FAST=true)
% Demo visualizing high dimensional vector magnitudes in lodim (e.g., 2D),
% for Gaussian and Uniform n-ball distributions of varying dimensionality
% 
% 2018-04-19 AZ Created

%% Preliminaries
if ~exist('JITTER','var') || isempty(JITTER),   JITTER = false;  end
if ~exist('FAST'  ,'var') || isempty(FAST  ),   FAST   = true;   end
if FAST, IMSIZE = [200 200]; 
else     IMSIZE = [1800 1800];
end
if ~exist('lodim','var') || isempty(lodim),   lodim = 2;
elseif lodim==3                               ds(1) = 3;
end
ds = 2.^(1:2:8); n = 100;
types = {'Uniform','Gaussian'};
colors = {[0 0 0],[1 0 0]};
switch lodim
   case 2,  bign = 10000;
   case 3,  bign = n;
end

%% Draw loop
   figure(1001);    set(gcf,'Name',sprintf('%usamples' ,n),'Renderer','painters');
if lodim==2
   figure(1002);clf;set(gcf,'Name',sprintf('%usamplesD',n),'Renderer','painters');
end
for t = 1:2
   i = 0;
   for d = ds
      i = i+1;
      switch t
         case 1,   X = randnball(bign,d);
                   r = 1;
         case 2,   X = randn(    bign,d);
                   r = sqrt(2)*exp(gammaln((d+1)/2)-gammaln(d/2));
      end
      if JITTER,   X = X - repmat(mean(X(1:n,:)),[bign 1]);      end

      H = hyperdist2dist(X,lodim);
      figure(1001);subplot(3,4,i+(t-1)*4);
      switch lodim
         case 2,   plot( H(1:n,1),H(1:n,2)         ,'ko'); drawCircle;
         case 3,   plot3(H(1:n,1),H(1:n,2),H(1:n,3),'ko'); draw3dEllipsoid;
                   lighting gouraud;camlight;
      end
      if     t==2 && lodim==2,   drawCircle(     [],r); 
      elseif t==2 && lodim==3,   draw3dEllipsoid([],r*eye(3));
      end
      title(sprintf('%uD %s',d,types{t}))
      if lodim==2 % 2D "density" histogram of transparent points
         figure(1002); subplot(3,4,i+(t-1)*4);
         scatterDensity(H,IMSIZE,true,0.01); axis equal off;
         title(sprintf('%uD %s',d,types{t}))
      end
      % Histogram of radii
      radii = sqrt(sum(X.^2,2));
      maxbin = max(radii/r)+0.2;
      bins = 0:0.01:maxbin;
      h = hist(radii/r,bins);
      figure(1001); subplot(3,4,i+8); hold on;
      plot([-flip(bins) bins],[flip(h/sum(h)) h/sum(h)],'Color',colors{t})
      xlim([-maxbin maxbin])
   end
end

printFig(1001,[],'eps');
printFig(1002,[],'png',1200);

