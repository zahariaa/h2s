function h2s_mnist_compare(epsilon)
% h2s_mnist:
% runs hypersphere2sphere on MNIST raw images and trained NN activations
% 
% 2018-03-02 AZ Created

if ~nargin, epsilon = 0; end

dimLow = 2;
model  = 'lenet';
%% MNIST, trained model responses, from pyTorch
data = load(sprintf('demo/mnist_%s_eps%0.2f.mat',model,epsilon)); % Trained net on MNIST
label2resp = @(L,chunk) cell2mat_concat(arrayfun(@(i) L(chunk)==i,0:9,'UniformOutput',false))';

chunk = 1:10000;

%% Load raw MNIST (probably not aligned to pyTorch stuff)
data.x{1}   = indexm(loadMNISTImages('data/MNIST/raw/t10k-images-idx3-ubyte'),[],chunk)';
data.labels = {loadMNISTLabels('data/MNIST/raw/t10k-labels-idx1-ubyte')'; data.labels};

%% reorganize data structure
nlayers = numel(fieldnames(data))-2;
for i = 1:nlayers
   eval(sprintf('data.x{%u} = data.x%u;data = rmfield(data,"x%u");',i+1,i,i));
end
%% compute final output
data.x{nlayers+2} = label2resp(maxix(data.x{nlayers+1}(chunk,:),[],2)'-1,chunk);

h2sfile = sprintf('demo/mnist_%s_eps%0.2f_h2s.mat',model,epsilon);
if exist(h2sfile,'file'), load(h2sfile); fprintf('Loaded h2s from %s\n',h2sfile);
else
   %% Construct labels
   % right now this uses the true labels
   for i = 1:2
      cats(i).(model) = Categories(label2resp(data.labels{i},chunk),...
                                   arrayfun(@num2str,0:9,'UniformOutput',false)');
   end
   % % output labels
   % cats.(model) = Categories(data.x{end},...
   %                         arrayfun(@num2str,0:9,'UniformOutput',false)');
   
   %% Estimate hyperspheres, run h2s
   for i = 1:nlayers+2
      stationarycounter(i,nlayers+1)
      hi.(model)(i) = Hypersphere.estimate(double(data.x{i}(chunk,:)),cats(~~(i-1)+1).(model),10000).meanAndMerge;
      lo.(model)(i) = hi.(model)(i).h2s(dimLow);
      if all(lo.(model)(i).radii==0)
         lo.(model)(i).radii = 0.05*ones(1,10);
      end
      hi.(model)(i) = hi.(model)(i).dropBootstraps;
      save(h2sfile,'hi','lo','cats','nlayers')
      fprintf('Saved h2s to %s     \n',h2sfile);
   end
end

%% Compute alternative visualizations for supplementary figure (from H2SandMDS)
[nPoints,nCats] =size(cats(end).(model).vectors);
nPointsPerCat   = sum(cats(end).(model).vectors);
if ~exist('points2D','var')
   for i = 1:nlayers+2
      stationarycounter(i,nlayers+2)
      dists{i} = pdist(data.x{i}(chunk,:),'Euclidean');
      if i>nlayers+1
         dists{i}(~dists{i}) = 1e-3*rand(1,sum(~dists{i}));
      end
      points2D{i}{1} = simplepca(data.x{i}(chunk,:),dimLow);
      points2D{i}{2} = mdscale(dists{i},2,'criterion','metricstress');
      points2D{i}{3} = tsne(single(data.x{i}(chunk,:)));%,[],[],min(size(data.x{i}(chunk,:))));
   end   

%% Draw points from circle/sphere to align other visualizations to
   modelPointsAlign = NaN(0,dimLow);
   % if all vectors are less than nPoints, add more proportionally
   if sum(nPointsPerCat)<nPoints
      nPointsPerCat = nPointsPerCat + floor((nPoints-sum(nPointsPerCat))*nPointsPerCat/sum(nPointsPerCat));
      catI = 1;
      while sum(nPointsPerCat)<nPoints
         nPointsPerCat(catI) = nPointsPerCat(catI) + 1;
         catI=catI+1;
      end
   end
   i=nlayers+1;
   for catI = 1:nCats
      modelPointsAlign = [modelPointsAlign;
                          randsphere(nPointsPerCat(catI),dimLow,lo.(model)(i).radii(catI)) ...
                          + repmat(lo.(model)(i).centers(catI,:),[nPointsPerCat(catI) 1])];
   end

%% Align PCA, MDS, and t-SNE relative to H2S
   for i = 1:nlayers+2
      if size(data.x{i}(chunk,:),2) > 1
         for j = 1:numel(points2D{i})
            if numel(points2D{i}{j})
               [~,points2D{i}{j}] = procrustes(modelPointsAlign,points2D{i}{j});
            end
         end
      end   
   end
   save(h2sfile,'hi','lo','cats','nlayers','points2D')
   fprintf('Saved h2s to %s     \n',h2sfile);
end
keyboard

%% PLOT
fh = newfigure([3 nlayers+2],[],[],sprintf( 'combo_%s_eps%0.2f_h2s%u',model,epsilon,dimLow));
fh = newfigure([1 nlayers+2],fh,[],sprintf('stats2_%s_eps%0.2f_h2s%u',model,epsilon,dimLow));
for i = 1:nlayers+2
   if i==1,            titlestr = 'Raw images';
   elseif i<nlayers+2, titlestr = sprintf('Layer %u',i-1);
   else                titlestr = 'Final readout';
   end
   lo.(model)(i).show(fh.a(1).h(i));                               title(titlestr);
   lo.(model)(i).showValues(fh.a(1).h(i+nlayers+2),hi.(model)(i));
   hi.(model)(i).showSig(fh.a(1).h(i+2*(nlayers+2)));
   hi.(model)(i).showSig(fh.a(2).h(i),hi.sec(i));                  title(titlestr);
end

%% Supplementary figure
ms=5;
chunk = 1:1000;
fh = newfigure([4 nlayers+2],fh,[],sprintf('compare_%s_eps%0.2f_h2s%u',model,epsilon,dimLow));
for i = 1:nlayers+2
   if i==1,            titlestr = 'Raw images';
   elseif i<nlayers+2, titlestr = sprintf('Layer %u',i-1);
   else                titlestr = 'Final readout';
   end
   lo.(model)(i).show(fh.a(end).h(i));   title(titlestr);
   for j = 1:3
      axtivate(fh.a(end).h(i+j*(nlayers+2)));
      for catI = 1:nCats
         if numel(points2D{i}{j})
             plot(points2D{i}{j}(logical(cats(~~(i-1)+1).(model).vectors(chunk,catI)),1), ...
                  points2D{i}{j}(logical(cats(~~(i-1)+1).(model).vectors(chunk,catI)),2),...
                 'o', 'MarkerFaceColor',cats(~~(i-1)+1).(model).colors(catI,:), ...
                 'MarkerEdgeColor','none', 'MarkerSize',ms);
          end
      end
      axis equal square off
   end
end

if dimLow==2,   printFig(fh.f,[],'eps')
else            printFig(fh.f(2:end),[],'eps'); printFig(fh.f(1),[],'png',300)
end

